#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

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

def main():
    parser = argparse.ArgumentParser(description="Filter segment alignments using Nextclade HA QC report.")
    parser.add_argument("--input-dir", required=True, help="Directory containing input alignments")
    parser.add_argument("--report", required=True, help="Path to Nextclade HA report TSV")
    parser.add_argument("--output-dir", required=True, help="Directory to write filtered alignments")
    parser.add_argument("--discarded-csv", default=None, help="Path to write CSV listing discarded sequences and reasons")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Load Nextclade report
    if not os.path.exists(args.report):
        print(f"Error: Nextclade report {args.report} does not exist.")
        sys.exit(1)

    df = pd.read_csv(args.report, sep="\t")
    
    # 2. Identify failed sequences (where totalFrameShifts > 0, qc.stopCodons.totalStopCodons > 0, or qc.overallStatus == 'bad')
    discard_set = set()
    discarded_records = []
    
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
            is_local_core = seq_name.lower().startswith("flu")
            reason_str = " & ".join(reasons)
            
            if is_local_core:
                discarded_records.append({
                    "sequence": seq_name,
                    "status": "KEPT",
                    "reason": f"Local Core - Kept despite: {reason_str}"
                })
                print(f"Sequence '{seq_name}' kept (Local Core) despite QC issues: {reason_str}")
            else:
                discard_set.add(seq_name)
                discarded_records.append({
                    "sequence": seq_name,
                    "status": "DISCARDED",
                    "reason": reason_str
                })
                print(f"Sequence '{seq_name}' discarded due to: {reason_str}")

    print(f"Total sequences discarded based on HA Nextclade QC: {len(discard_set)}")
    if discard_set:
        print(f"Discarded IDs: {sorted(list(discard_set))}")

    # Write discarded CSV if requested
    if args.discarded_csv:
        os.makedirs(os.path.dirname(args.discarded_csv), exist_ok=True)
        # Sort records by sequence ID
        discarded_records = sorted(discarded_records, key=lambda x: (x["status"], x["sequence"]))
        pd.DataFrame(discarded_records).to_csv(args.discarded_csv, index=False)
        print(f"Written QC report CSV: {args.discarded_csv}")

    # 3. Filter each fasta alignment file in the input directory
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
            if clean_header in discard_set:
                continue
                
            filtered_records.append((header, seq))
            kept_count += 1
            
        write_fasta(out_path, filtered_records)
        print(f"Filtered {file_name}: kept {kept_count}/{total_count} sequences (discarded {total_count - kept_count})")

if __name__ == "__main__":
    main()
