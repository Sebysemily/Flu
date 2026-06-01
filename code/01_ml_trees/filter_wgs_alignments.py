#!/usr/bin/env python3
"""Subset all 8 QC alignments to taxa that are fully complete (WGS)."""

import argparse
import csv
import os
import sys

import pandas as pd

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

from build_inputs.process_raw_to_segments import read_fasta, wrap_seq

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

def complete_wgs_ecuador_epis(filtrado_csv: str) -> set[str]:
    df = pd.read_csv(filtrado_csv, dtype=str)
    epis = set()
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    for _, row in df.iterrows():
        epi_raw = row.get("EPI_ISL")
        if pd.isna(epi_raw):
            continue
        epi = str(epi_raw).strip()
        
        is_complete = True
        for seg in segments:
            if str(row.get(seg, "")).strip().upper() != "SI":
                is_complete = False
                break
        if is_complete:
            epis.add(epi)
    return epis

def main() -> None:
    parser = argparse.ArgumentParser(description="Filter alignments to keep complete WGS Ecuador + all context.")
    parser.add_argument("--alignments", nargs="+", required=True)
    parser.add_argument("--flu-filtrado", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--audit-csv", required=True)
    args = parser.parse_args()

    wgs_ecuador = complete_wgs_ecuador_epis(args.flu_filtrado)
    df = pd.read_csv(args.flu_filtrado, dtype=str)
    all_ecuador_epis = set(df["EPI_ISL"].dropna().astype(str).str.strip())
    
    all_ids = set()
    for aln in args.alignments:
        all_ids.update(load_alignment_ids(aln))
        
    context_taxa = all_ids - all_ecuador_epis
    
    keep_taxa = context_taxa | wgs_ecuador
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    kept_counts = {}
    for aln in args.alignments:
        seg = os.path.basename(aln).split(".")[0]
        out_path = os.path.join(args.output_dir, os.path.basename(aln))
        kept_counts[seg] = write_subset_fasta(aln, out_path, keep_taxa)
        
    audit_dir = os.path.dirname(args.audit_csv)
    if audit_dir:
        os.makedirs(audit_dir, exist_ok=True)
    with open(args.audit_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerow({"metric": "total_input_taxa", "value": len(all_ids)})
        writer.writerow({"metric": "ecuador_wgs_filtrado", "value": len(wgs_ecuador)})
        writer.writerow({"metric": "kept_total", "value": len(keep_taxa)})
        for seg, count in kept_counts.items():
            writer.writerow({"metric": f"kept_{seg}_records", "value": count})

    print(f"Wrote WGS alignments to {args.output_dir}: kept {len(keep_taxa)} max taxa (Ecuador WGS={len(wgs_ecuador)})")

if __name__ == "__main__":
    main()
