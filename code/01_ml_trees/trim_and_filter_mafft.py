#!/usr/bin/env python3
"""Trim alignment to CDS ORF and drop sequences with excess gaps/ambiguities."""
import argparse
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qc.flu_role_utils import load_role_map, write_discarded_rows

DEFAULT_MAX_DIVERGENCE = 0.10


def read_fasta(path: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")


def load_protected_ids(path: str | None) -> set[str]:
    if not path:
        return set()
    with open(path, encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def gap_n_fraction(seq: str, valid_len: int) -> float:
    if valid_len <= 0:
        return 0.0
    bad = sum(1 for base in seq if base in "-N")
    return bad / valid_len


def get_consensus_and_trim_limits(records: list[tuple[str, str]]) -> tuple[int, int, np.ndarray]:
    mat = np.array([list(s.upper()) for _, s in records])
    aln_len = mat.shape[1]
    consensus = np.array(
        [Counter(mat[:, i]).most_common(1)[0][0] for i in range(aln_len)]
    )

    ungapped_cons: list[str] = []
    col_map: list[int] = []
    for i, char in enumerate(consensus):
        if char not in ("-", "N"):
            ungapped_cons.append(char)
            col_map.append(i)

    ungapped_str = "".join(ungapped_cons)
    best_start, best_end, max_len = -1, -1, -1
    for i in range(len(ungapped_str) - 2):
        if ungapped_str[i : i + 3] != "ATG":
            continue
        for j in range(i + 3, len(ungapped_str) - 2, 3):
            if ungapped_str[j : j + 3] in ("TAA", "TAG", "TGA"):
                orf_len = j + 3 - i
                if orf_len > max_len:
                    max_len, best_start, best_end = orf_len, i, j + 3
                break

    if best_start == -1:
        print("Warning: No clear ORF found. Falling back to whole alignment.")
        return 0, aln_len, consensus

    return col_map[best_start], col_map[best_end - 1] + 1, consensus


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument(
        "--max-divergence",
        type=float,
        default=DEFAULT_MAX_DIVERGENCE,
        help="Max fraction of gaps/N in trimmed CDS (default 0.10)",
    )
    parser.add_argument(
        "--protect-ids",
        default=None,
        help="Taxon IDs kept regardless of gap/N rate (one per line)",
    )
    parser.add_argument(
        "--qc-audit",
        default=None,
        help="Optional TSV of protected taxa kept above --max-divergence",
    )
    parser.add_argument(
        "--discarded-csv",
        default=None,
        help="CSV of sequences dropped by gap/N threshold",
    )
    parser.add_argument(
        "--role-metadata",
        default=None,
        help="Metadata CSV with file_name and expected_role (e.g. H5N1_context.csv)",
    )
    parser.add_argument(
        "--filter-step",
        default="trim_gap_n",
        help="Label written to discarded CSV (e.g. trim_gap_n_HA, trim_gap_n_panel_HA)",
    )
    parser.add_argument(
        "--discarded-csv-only",
        action="store_true",
        help="Only write discarded CSVs; do not read/write --output alignment",
    )
    args = parser.parse_args()
    if not args.discarded_csv_only and not args.output:
        parser.error("--output is required unless --discarded-csv-only is set")
    if args.discarded_csv_only and not args.discarded_csv:
        parser.error("--discarded-csv-only requires --discarded-csv")
    return args


DISCARD_FIELDS = [
    "taxon",
    "expected_role",
    "gap_n_fraction",
    "max_divergence",
    "filter_step",
    "discard_reason",
]


def evaluate_trim(
    records: list[tuple[str, str]],
    max_divergence: float,
    protected: set[str],
) -> tuple[list[tuple[str, str]], list[dict], list[tuple[str, float]]]:
    start_col, end_col, consensus = get_consensus_and_trim_limits(records)
    trimmed_cons = consensus[start_col:end_col]
    valid_len = len(trimmed_cons)
    print(f"Trimming columns {start_col}:{end_col} (CDS length {valid_len})")

    kept: list[tuple[str, str]] = []
    dropped_rows: list[dict] = []
    protected_kept: list[tuple[str, float]] = []

    for header, seq in records:
        trimmed = seq[start_col:end_col].upper()
        div = gap_n_fraction(trimmed, valid_len)
        if header in protected or div <= max_divergence:
            kept.append((header, trimmed))
            if header in protected and div > max_divergence:
                protected_kept.append((header, div))
        else:
            dropped_rows.append(
                {
                    "taxon": header,
                    "gap_n_fraction": f"{div:.4f}",
                    "max_divergence": max_divergence,
                    "discard_reason": f"gaps/N > {max_divergence:.0%} in trimmed CDS",
                }
            )

    return kept, dropped_rows, protected_kept


def main() -> None:
    args = parse_args()
    protected = load_protected_ids(args.protect_ids)
    role_map = load_role_map(args.role_metadata)

    records = read_fasta(args.input)
    if not records:
        print("Empty input.")
        if not args.discarded_csv_only and args.output:
            write_fasta(args.output, [])
        return

    kept, dropped_rows, protected_kept = evaluate_trim(
        records, args.max_divergence, protected
    )

    for row in dropped_rows:
        row["expected_role"] = role_map.get(row["taxon"], "")
        row["filter_step"] = args.filter_step

    flu_dropped = [r for r in dropped_rows if r["expected_role"].startswith("flu_")]
    print(
        f"Kept {len(kept)} / {len(records)} sequences; dropped {len(dropped_rows)} "
        f"(>{args.max_divergence:.0%} gaps/N); flu_* dropped: {len(flu_dropped)}."
    )
    if protected:
        print(f"Protected set size: {len(protected)}; kept above threshold: {len(protected_kept)}")

    if args.discarded_csv:
        write_discarded_rows(args.discarded_csv, dropped_rows, DISCARD_FIELDS)
        print(f"Written discarded CSV: {args.discarded_csv} ({len(dropped_rows)} rows)")

    if args.discarded_csv_only:
        return

    write_fasta(args.output, kept)

    if args.qc_audit:
        import csv

        audit_dir = os.path.dirname(args.qc_audit)
        if audit_dir:
            os.makedirs(audit_dir, exist_ok=True)
        with open(args.qc_audit, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["taxon", "gap_n_fraction", "note"])
            for taxon, div in protected_kept:
                writer.writerow([taxon, f"{div:.4f}", "protected_kept"])


if __name__ == "__main__":
    main()
