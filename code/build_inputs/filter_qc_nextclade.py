#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qc.flu_role_utils import load_role_map, write_discarded_rows

NEXTCLADE_DISCARD_FIELDS = [
    "seqName",
    "expected_role",
    "qc_action",
    "discard_reason",
    "filter_step",
]

def read_fasta(path):
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)

def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as out_fh:
        for header, seq in records:
            out_fh.write(f">{header}\n")
            width = 80
            wrapped = "\n".join(seq[i : i + width] for i in range(0, len(seq), width))
            out_fh.write(wrapped + "\n")

def load_core_ids(metadata_path):
    import csv
    import re
    core_ids = set()
    if not metadata_path or not os.path.exists(metadata_path):
        return core_ids
    with open(metadata_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("EPI_ISL") or row.get("file_name") or list(row.values())[0]).strip()
            role = row.get("expected_role")
            if role is not None:
                if name and role.strip() and re.match(r"^flu_.*", role.strip(), re.IGNORECASE):
                    core_ids.add(name)
            else:
                # If expected_role is not present (e.g. flu_filtrado.csv), assume all are core ids
                if name:
                    core_ids.add(name)
    return core_ids

def main():
    parser = argparse.ArgumentParser(description="Filter segment alignments using Nextclade QC report.")
    parser.add_argument("--input-alignment", help="Path to single input alignment file")
    parser.add_argument("--output-alignment", help="Path to single output alignment file")
    parser.add_argument("--input-dir", help="Directory containing input alignments")
    parser.add_argument("--output-dir", help="Directory to write filtered alignments")
    parser.add_argument("--report", required=True, help="Path to Nextclade report TSV or CSV")
    parser.add_argument("--metadata", default=None, help="Path to metadata CSV to protect core flu sequences")
    parser.add_argument("--discarded-csv", default=None, help="Path to write CSV listing discarded sequences and reasons")
    parser.add_argument(
        "--role-metadata",
        default=None,
        help="Metadata with file_name and expected_role (e.g. H5N1_context.csv)",
    )
    parser.add_argument(
        "--filter-step",
        default="nextclade_ha",
        help="Label for filter_step column in discarded CSV",
    )
    parser.add_argument("--skip-filter", action="store_true", help="Do not filter out any sequences from alignment, but still generate reports")
    parser.add_argument(
        "--discarded-csv-only",
        action="store_true",
        help="Only write discarded/flu-discarded CSVs; do not read or write alignment FASTA",
    )
    args = parser.parse_args()


    if args.discarded_csv_only:
        if not args.discarded_csv:
            print("Error: --discarded-csv-only requires --discarded-csv.")
            sys.exit(1)
    else:
        if (args.input_alignment or args.output_alignment) and not (
            args.input_alignment and args.output_alignment
        ):
            print("Error: Both --input-alignment and --output-alignment must be specified if one is.")
            sys.exit(1)
        if not args.input_alignment and not args.input_dir:
            print("Error: Must specify either --input-alignment or --input-dir.")
            sys.exit(1)

    # 1. Load Nextclade report and core ids
    if not os.path.exists(args.report):
        print(f"Error: Nextclade report {args.report} does not exist.")
        sys.exit(1)

    sep = None
    if args.report.endswith(".csv"):
        # Detect if nextclade used semicolon or comma
        with open(args.report, "r", encoding="utf-8") as f:
            first_line = f.readline()
        if ";" in first_line:
            sep = ";"
        elif "," in first_line:
            sep = ","
    
    if sep is None:
        sep = "," if args.report.endswith(".csv") else "\t"

    df = pd.read_csv(args.report, sep=sep)
    core_ids = load_core_ids(args.metadata)
    role_map = load_role_map(args.role_metadata)
    
    # 2. Identify failed sequences (where totalFrameShifts > 0, qc.stopCodons.totalStopCodons > 0, or qc.overallStatus == 'bad')
    discard_set = set()
    failed_rows = []
    
    frameshift_col = "totalFrameShifts" if "totalFrameShifts" in df.columns else None
    stopcodon_col = "qc.stopCodons.totalStopCodons" if "qc.stopCodons.totalStopCodons" in df.columns else None
    overall_status_col = "qc.overallStatus" if "qc.overallStatus" in df.columns else None
            
    if not frameshift_col:
        print("Warning: totalFrameShifts column not found in Nextclade report. Skipping frameshift check.")
    if not stopcodon_col:
        print("Warning: totalStopCodons column not found in Nextclade report. Skipping stop codon check.")
    if not overall_status_col:
        print("Warning: qc.overallStatus column not found in Nextclade report. Skipping overall status check.")
        
    for idx, row in df.iterrows():
        seq_name = str(row["seqName"]).strip()
        failed = False
        reasons = []
        
        if frameshift_col:
            val = row[frameshift_col]
            try:
                if pd.notna(val) and float(val) > 0:
                    failed = True
                    reasons.append(f"frameshifts ({int(float(val))})")
            except ValueError:
                pass
                
        if stopcodon_col:
            val = row[stopcodon_col]
            try:
                if pd.notna(val) and float(val) > 0:
                    failed = True
                    reasons.append(f"stop codons ({int(float(val))})")
            except ValueError:
                pass
                
        if overall_status_col:
            val = str(row[overall_status_col]).strip().lower()
            if val in ("bad", "mediocre"):
                failed = True
                reasons.append(f"{val} overall status")
                
        if failed:
            is_local_core = (seq_name in core_ids)
            reason_str = " & ".join(reasons)
            
            row_dict = row.to_dict()
            if is_local_core:
                row_dict["qc_action"] = "KEPT"
                row_dict["discard_reason"] = f"Local Core - Kept despite: {reason_str}"
                print(f"Sequence '{seq_name}' kept (Local Core) despite QC issues: {reason_str}")
            elif args.skip_filter:
                row_dict["qc_action"] = "KEPT"
                row_dict["discard_reason"] = f"Filter Skipped - Kept despite: {reason_str}"
                print(f"Sequence '{seq_name}' kept (Filter Skipped) despite QC issues: {reason_str}")
            else:
                discard_set.add(seq_name)
                row_dict["qc_action"] = "DISCARDED"
                row_dict["discard_reason"] = reason_str
                print(f"Sequence '{seq_name}' discarded due to: {reason_str}")
            failed_rows.append(row_dict)

    summary_rows = []
    for row_dict in failed_rows:
        seq_name = str(row_dict.get("seqName", "")).strip()
        role = role_map.get(seq_name, "")
        summary_rows.append(
            {
                "seqName": seq_name,
                "expected_role": role,
                "qc_action": row_dict.get("qc_action", ""),
                "discard_reason": row_dict.get("discard_reason", ""),
                "filter_step": args.filter_step,
            }
        )

    print(f"Total sequences discarded based on Nextclade QC: {len(discard_set)}")
    if discard_set:
        print(f"Discarded IDs: {sorted(list(discard_set))}")

    discarded_only = [r for r in summary_rows if r.get("qc_action") == "DISCARDED"]
    if args.discarded_csv:
        write_discarded_rows(args.discarded_csv, discarded_only, NEXTCLADE_DISCARD_FIELDS)
        print(f"Written QC report CSV: {args.discarded_csv} ({len(discarded_only)} discarded)")
    if args.discarded_csv_only:
        return

    # 3. Filter alignment files
    if args.input_alignment and args.output_alignment:
        in_path = args.input_alignment
        out_path = args.output_alignment
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        
        filtered_records = []
        total_count = 0
        kept_count = 0
        
        for header, seq in read_fasta(in_path):
            total_count += 1
            clean_header = header.strip()
            seq_id = clean_header.split()[0]
            if clean_header in discard_set or seq_id in discard_set:
                continue
                
            filtered_records.append((header, seq))
            kept_count += 1
            
        write_fasta(out_path, filtered_records)
        print(f"Filtered {in_path}: kept {kept_count}/{total_count} sequences (discarded {total_count - kept_count})")
    
    elif args.input_dir and args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        for file_name in os.listdir(args.input_dir):
            if not file_name.endswith(".mafft"):
                continue
                
            in_path = os.path.join(args.input_dir, file_name)
            out_path = os.path.join(args.output_dir, file_name)
            
            filtered_records = []
            total_count = 0
            kept_count = 0
            
            for header, seq in read_fasta(in_path):
                total_count += 1
                clean_header = header.strip()
                seq_id = clean_header.split()[0]
                if clean_header in discard_set or seq_id in discard_set:
                    continue
                    
                filtered_records.append((header, seq))
                kept_count += 1
                
            write_fasta(out_path, filtered_records)
            print(f"Filtered {file_name}: kept {kept_count}/{total_count} sequences (discarded {total_count - kept_count})")

if __name__ == "__main__":
    main()
