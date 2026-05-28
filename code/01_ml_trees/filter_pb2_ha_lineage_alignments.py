#!/usr/bin/env python3
"""Subset PB2/HA QC alignments to taxa with both segments (flu_filtrado + context)."""

from __future__ import annotations

import argparse
import csv
import os
import sys

import pandas as pd

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

from build_inputs.process_raw_to_segments import read_fasta, wrap_seq  # noqa: E402


def load_alignment_ids(path: str) -> set[str]:
    return {header.split()[0] for header, _seq in read_fasta(path)}


def write_subset_fasta(in_path: str, out_path: str, keep_taxa: set[str]) -> int:
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    kept = 0
    with open(out_path, "w", encoding="utf-8") as handle:
        for header, seq in read_fasta(in_path):
            taxon = header.split()[0]
            if taxon not in keep_taxa:
                continue
            handle.write(f">{header}\n")
            handle.write(wrap_seq(seq) + "\n")
            kept += 1
    return kept


def dual_segment_ecuador_epis(filtrado_csv: str) -> set[str]:
    df = pd.read_csv(filtrado_csv, dtype=str)
    epis = set()
    for _, row in df.iterrows():
        epi_raw = row.get("EPI_ISL")
        if pd.isna(epi_raw):
            continue
        epi = str(epi_raw).strip()
        pb2 = str(row.get("PB2", "")).strip().upper()
        ha = str(row.get("HA", "")).strip().upper()
        if pb2 == "SI" and ha == "SI":
            epis.add(epi)
    return epis


def select_lineage_taxa(
    pb2_ids: set[str],
    ha_ids: set[str],
    filtrado_csv: str,
) -> set[str]:
    in_both = pb2_ids & ha_ids
    df = pd.read_csv(filtrado_csv, dtype=str)
    ecuador_epis = set(df["EPI_ISL"].dropna().astype(str).str.strip())
    dual_ecuador = dual_segment_ecuador_epis(filtrado_csv)
    context_taxa = in_both - ecuador_epis
    ecuador_taxa = in_both & dual_ecuador
    return context_taxa | ecuador_taxa


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter PB2/HA alignments for lineage analysis (dual-segment Ecuador + context)."
    )
    parser.add_argument("--pb2-alignment", required=True)
    parser.add_argument("--ha-alignment", required=True)
    parser.add_argument("--flu-filtrado", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--audit-csv", required=True)
    args = parser.parse_args()

    pb2_ids = load_alignment_ids(args.pb2_alignment)
    ha_ids = load_alignment_ids(args.ha_alignment)
    keep_taxa = select_lineage_taxa(pb2_ids, ha_ids, args.flu_filtrado)

    pb2_out = os.path.join(args.output_dir, "H5N1_PB2.mafft")
    ha_out = os.path.join(args.output_dir, "H5N1_HA.mafft")
    n_pb2 = write_subset_fasta(args.pb2_alignment, pb2_out, keep_taxa)
    n_ha = write_subset_fasta(args.ha_alignment, ha_out, keep_taxa)

    dual_ecuador = dual_segment_ecuador_epis(args.flu_filtrado)
    df = pd.read_csv(args.flu_filtrado, dtype=str)
    ecuador_epis = set(df["EPI_ISL"].dropna().astype(str).str.strip())
    context_kept = keep_taxa - ecuador_epis
    ecuador_kept = keep_taxa & dual_ecuador

    audit_dir = os.path.dirname(args.audit_csv)
    if audit_dir:
        os.makedirs(audit_dir, exist_ok=True)
    with open(args.audit_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerow({"metric": "pb2_alignment_taxa", "value": len(pb2_ids)})
        writer.writerow({"metric": "ha_alignment_taxa", "value": len(ha_ids)})
        writer.writerow({"metric": "intersection_taxa", "value": len(pb2_ids & ha_ids)})
        writer.writerow({"metric": "ecuador_dual_pb2_ha_filtrado", "value": len(dual_ecuador)})
        writer.writerow({"metric": "kept_total", "value": len(keep_taxa)})
        writer.writerow({"metric": "kept_ecuador", "value": len(ecuador_kept)})
        writer.writerow({"metric": "kept_context", "value": len(context_kept)})
        writer.writerow({"metric": "kept_pb2_records", "value": n_pb2})
        writer.writerow({"metric": "kept_ha_records", "value": n_ha})

    print(
        f"Wrote lineage alignments to {args.output_dir}: "
        f"PB2={n_pb2}, HA={n_ha} taxa (Ecuador={len(ecuador_kept)}, context={len(context_kept)})"
    )


if __name__ == "__main__":
    main()
