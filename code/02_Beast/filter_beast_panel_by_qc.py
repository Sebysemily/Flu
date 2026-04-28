#!/usr/bin/env python3
import argparse
import csv
import os


def read_tsv_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(value, default=0.0):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def is_true_flag(value):
    text = str(value).strip().lower()
    return text in {"1", "true", "yes"}


def main():
    parser = argparse.ArgumentParser(
        description="Filter BEAST panel taxa using conservative QC-based exclusions"
    )
    parser.add_argument("--panel-taxa", required=True)
    parser.add_argument("--source-qc-metrics", required=True)
    parser.add_argument("--filtered-panel-out", required=True)
    parser.add_argument("--exclusions-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--high-n-threshold", type=float, default=0.10)
    args = parser.parse_args()

    panel_rows = read_tsv_rows(args.panel_taxa)
    qc_rows = read_tsv_rows(args.source_qc_metrics)
    qc_by_taxon = {
        str(row.get("taxon", "")).strip(): row
        for row in qc_rows
        if str(row.get("taxon", "")).strip()
    }

    filtered_rows = []
    excluded_rows = []

    for row in panel_rows:
        taxon = str(row.get("taxon", "")).strip()
        qc = qc_by_taxon.get(taxon, {})
        high_n_flag = is_true_flag(qc.get("high_raw_n_fraction", "0"))
        raw_total_n_fraction = as_float(qc.get("raw_total_n_fraction", 0.0))
        outlier_score = as_float(qc.get("outlier_score", 0.0))
        flags = str(qc.get("flags", "")).strip()

        excluded = high_n_flag and raw_total_n_fraction >= args.high_n_threshold
        if excluded:
            excluded_rows.append(
                {
                    "taxon": taxon,
                    "excluded": "true",
                    "reason_code": "high_raw_n_fraction_extreme",
                    "reason_detail": (
                        f"high_raw_n_fraction=1 and raw_total_n_fraction={raw_total_n_fraction:.6f} "
                        f">= {args.high_n_threshold:.2f}"
                    ),
                    "raw_total_n_fraction": f"{raw_total_n_fraction:.6f}",
                    "outlier_score": f"{outlier_score:.6f}",
                    "flags": flags,
                }
            )
        else:
            filtered_rows.append(row)

    write_tsv(
        args.filtered_panel_out,
        filtered_rows,
        ["taxon", "role", "lineage", "distance_to_seed"],
    )
    write_tsv(
        args.exclusions_out,
        excluded_rows,
        [
            "taxon",
            "excluded",
            "reason_code",
            "reason_detail",
            "raw_total_n_fraction",
            "outlier_score",
            "flags",
        ],
    )
    write_tsv(
        args.summary_out,
        [
            {"metric": "panel_taxa_input", "value": len(panel_rows)},
            {"metric": "panel_taxa_excluded", "value": len(excluded_rows)},
            {"metric": "panel_taxa_kept", "value": len(filtered_rows)},
            {"metric": "high_n_threshold", "value": f"{args.high_n_threshold:.2f}"},
        ],
        ["metric", "value"],
    )


if __name__ == "__main__":
    main()
