#!/usr/bin/env python3
"""Build HA panel taxon list: mandatory flu+outgroup, regional global, regional flu_costa window, american anchor.

``--out-include-list`` lists taxa excluded from Augur subsampling pools (always kept in
``context_taxa``). Panel trim/QC uses the global ``TRIM_MAX_DIVERGENCE`` only; flu are not exempt.
"""
import argparse
import csv
import os
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", required=True)
    p.add_argument("--ha-fasta", required=True)
    p.add_argument("--out-context-taxa", required=True)
    p.add_argument("--out-include-list", required=True)
    p.add_argument("--out-costa-audit", required=True)
    p.add_argument("--outgroup", required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--regional-context-max", type=int, default=200)
    p.add_argument("--american-anchor-max", type=int, default=100)
    p.add_argument("--regional-costa-adjacent-max", type=int, default=80)
    p.add_argument("--costa-window-padding-months", type=int, default=2)
    return p.parse_args()


def load_ha_ids(path: str) -> set[str]:
    ids: set[str] = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip())
    return ids


def parse_iso_date(value: str) -> datetime | None:
    value = (value or "").strip()
    if not value:
        return None
    try:
        if len(value) >= 10:
            return datetime.strptime(value[:10], "%Y-%m-%d")
        if len(value) >= 7:
            return datetime.strptime(value[:7], "%Y-%m")
        if len(value) >= 4:
            return datetime.strptime(value[:4], "%Y")
    except ValueError:
        return None
    return None


def add_months(dt: datetime, months: int) -> datetime:
    month_index = dt.month - 1 + months
    year = dt.year + month_index // 12
    month = month_index % 12 + 1
    return datetime(year, month, 1)


def write_lines(path: str, lines: list[str]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


def read_strain_list(path: str) -> list[str]:
    strains = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                strains.append(line)
    return strains


def run_augur_filter(
    metadata_csv: str,
    query: str,
    max_sequences: int,
    seed: int,
    output_strains: str,
) -> None:
    cmd = [
        "augur",
        "filter",
        "--metadata",
        metadata_csv,
        "--metadata-id-columns",
        "file_name",
        "--query",
        query,
        "--group-by",
        "country",
        "year_month",
        "--subsample-max-sequences",
        str(max_sequences),
        "--subsample-seed",
        str(seed),
        "--output-strains",
        output_strains,
    ]
    print("Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def compute_costa_window(rows: list[dict], ha_ids: set[str], padding_months: int) -> tuple[str, str, list[str]]:
    costa_dates: list[datetime] = []
    costa_ha_ids: list[str] = []
    for row in rows:
        if row.get("expected_role") != "flu_costa":
            continue
        name = (row.get("file_name") or "").strip()
        if name not in ha_ids:
            continue
        costa_ha_ids.append(name)
        dt = parse_iso_date(row.get("collection_date", ""))
        if dt:
            costa_dates.append(dt)

    if not costa_dates:
        raise ValueError("No flu_costa samples with HA and parseable collection_date found.")

    min_dt = min(costa_dates)
    max_dt = max(costa_dates)
    window_start_dt = add_months(datetime(min_dt.year, min_dt.month, 1), -padding_months)
    window_end_month_start = add_months(datetime(max_dt.year, max_dt.month, 1), padding_months)
    window_end_dt = add_months(window_end_month_start, 1) - timedelta(days=1)

    start_iso = window_start_dt.strftime("%Y-%m-%d")
    end_iso = window_end_dt.strftime("%Y-%m-%d")
    return start_iso, end_iso, costa_ha_ids


def main() -> None:
    args = parse_args()
    ha_ids = load_ha_ids(args.ha_fasta)

    with open(args.metadata, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    flu_ha = [
        (row.get("file_name") or "").strip()
        for row in rows
        if (row.get("expected_role") or "").startswith("flu_")
        and (row.get("file_name") or "").strip() in ha_ids
    ]
    protected = sorted(set([args.outgroup] + flu_ha))
    write_lines(args.out_include_list, protected)
    print(f"Mandatory panel strains (outgroup + flu with HA): {len(protected)}", file=sys.stderr)

    window_start, window_end, costa_ha_ids = compute_costa_window(
        rows, ha_ids, args.costa_window_padding_months
    )
    print(f"flu_costa HA window: {window_start} .. {window_end} ({len(costa_ha_ids)} costa taxa)", file=sys.stderr)

    audit_path = Path(args.out_costa_audit)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    with open(audit_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["flu_costa_ha_count", len(costa_ha_ids)])
        writer.writerow(["window_start", window_start])
        writer.writerow(["window_end", window_end])
        writer.writerow(["padding_months", args.costa_window_padding_months])
        writer.writerow(["regional_context_max", args.regional_context_max])
        writer.writerow(["regional_costa_adjacent_max", args.regional_costa_adjacent_max])
        writer.writerow(["american_anchor_max", args.american_anchor_max])

    # Metadata for augur: add year_month, exclude protected
    tmp_meta = str(audit_path.parent / "metadata_augur.tmp.csv")

    def year_month(d: str) -> str:
        d = (d or "").strip()
        if not d:
            return "unknown"
        parts = d.split("-")
        if len(parts) >= 2:
            return f"{parts[0]}-{parts[1]}"
        return parts[0]

    fieldnames = list(rows[0].keys()) + ["year_month"]
    protected_set = set(protected)
    with open(tmp_meta, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            if (row.get("file_name") or "").strip() in protected_set:
                continue
            row = dict(row)
            row["year_month"] = year_month(row.get("collection_date", ""))
            writer.writerow(row)

    work_dir = audit_path.parent
    rc_global = str(work_dir / "strains_rc_global.tmp")
    rc_costa = str(work_dir / "strains_rc_costa.tmp")
    aa_out = str(work_dir / "strains_aa.tmp")

    run_augur_filter(
        tmp_meta,
        "expected_role == 'regional_context'",
        args.regional_context_max,
        args.seed,
        rc_global,
    )

    costa_query = (
        f"expected_role == 'regional_context' and "
        f"collection_date >= '{window_start}' and collection_date <= '{window_end}'"
    )
    run_augur_filter(
        tmp_meta,
        costa_query,
        args.regional_costa_adjacent_max,
        args.seed + 1,
        rc_costa,
    )

    run_augur_filter(
        tmp_meta,
        "expected_role == 'american_anchor'",
        args.american_anchor_max,
        args.seed + 2,
        aa_out,
    )

    merged: list[str] = []
    seen: set[str] = set()
    for path in [args.out_include_list, rc_global, rc_costa, aa_out]:
        for strain in read_strain_list(path):
            if strain not in seen:
                seen.add(strain)
                merged.append(strain)

    write_lines(args.out_context_taxa, merged)
    print(f"Final context_taxa: {len(merged)} strains", file=sys.stderr)
    print(
        f"  protected={len(protected)}, rc_global={len(read_strain_list(rc_global))}, "
        f"rc_costa={len(read_strain_list(rc_costa))}, aa={len(read_strain_list(aa_out))}",
        file=sys.stderr,
    )

    for tmp in [tmp_meta, rc_global, rc_costa, aa_out]:
        try:
            os.remove(tmp)
        except OSError:
            pass


if __name__ == "__main__":
    main()
