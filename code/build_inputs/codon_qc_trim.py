#!/usr/bin/env python3
import argparse
import os
import sys

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

import csv

def wrap_seq(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))

def main():
    parser = argparse.ArgumentParser(description="Alignment-first codon-aware QC and reference-based trimming.")
    parser.add_argument("--alignment", required=True, help="Path to input aligned FASTA")
    parser.add_argument("--segment", required=True, help="Segment name (e.g. HA, PB2)")
    parser.add_argument("--output", required=True, help="Path to output trimmed/QC FASTA")
    parser.add_argument("--min-coverage", type=float, default=0.8, help="Minimum sequence coverage relative to reference length")
    parser.add_argument("--qc-summary", help="Path to output QC metrics CSV file")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # 1. Read all sequences from alignment
    aligned_records = []
    for header, seq in read_fasta(args.alignment):
        aligned_records.append((header, seq))

    if not aligned_records:
        print(f"Error: Alignment file {args.alignment} is empty.")
        sys.exit(1)

    # 2. Find reference sequence (first Ecuador core starting with 'Flu-')
    ref_header = None
    ref_seq = None
    for header, seq in aligned_records:
        if header.startswith("Flu-"):
            ref_header = header
            ref_seq = seq
            break

    # Fallback to first sequence if no Flu- sequence found
    if ref_seq is None:
        ref_header, ref_seq = aligned_records[0]
        print(f"Warning: No Ecuador core (Flu-*) sequence found in {args.alignment}. Using first sequence as reference: {ref_header}")

    # 3. Find CDS start and end columns based on reference non-gap characters
    # Find first non-gap character
    cds_start_idx = None
    for i, char in enumerate(ref_seq):
        if char != "-":
            cds_start_idx = i
            break

    # Find last non-gap character
    cds_end_idx = None
    for i in range(len(ref_seq) - 1, -1, -1):
        if ref_seq[i] != "-":
            cds_end_idx = i
            break

    if cds_start_idx is None or cds_end_idx is None:
        print(f"Error: Reference sequence {ref_header} contains only gaps.")
        sys.exit(1)

    ref_cds_length = sum(1 for char in ref_seq[cds_start_idx:cds_end_idx + 1] if char != "-")
    print(f"Segment: {args.segment}")
    print(f"Reference: {ref_header}")
    print(f"CDS Alignment Columns: {cds_start_idx} to {cds_end_idx} (Length: {cds_end_idx - cds_start_idx + 1})")
    print(f"Reference CDS Nucleotides (excluding gaps): {ref_cds_length} bp")

    # 4. Trim and filter sequences
    trimmed_records = []
    fragmented_count = 0
    qc_records = []

    for header, seq in aligned_records:
        # Extract sub-sequence corresponding to reference CDS columns
        trimmed_seq = seq[cds_start_idx : cds_end_idx + 1]
        
        # Count non-gap nucleotides in this range
        n_nucleotides = sum(1 for char in trimmed_seq if char != "-")
        
        # Check coverage
        coverage = n_nucleotides / ref_cds_length if ref_cds_length > 0 else 0
        passed = coverage >= args.min_coverage
        
        qc_records.append({
            "sequence": header,
            "segment": args.segment,
            "coverage": round(coverage, 4),
            "status": "PASS" if passed else "FAIL"
        })

        if not passed:
            fragmented_count += 1
            continue  # Discard fragmented sequence

        # Keep the trimmed sequence WITH its internal gaps intact to preserve frame
        trimmed_records.append((header, trimmed_seq))

    # 5. Write outputs
    with open(args.output, "w") as out_fh:
        for header, seq in trimmed_records:
            out_fh.write(f">{header}\n")
            out_fh.write(wrap_seq(seq) + "\n")

    # 6. Write CSV QC metrics
    if args.qc_summary:
        os.makedirs(os.path.dirname(args.qc_summary), exist_ok=True)
        with open(args.qc_summary, "w", newline="", encoding="utf-8") as csv_fh:
            writer = csv.DictWriter(csv_fh, fieldnames=["sequence", "segment", "coverage", "status"])
            writer.writeheader()
            writer.writerows(qc_records)
        print(f"Written QC metrics CSV: {args.qc_summary}")

    print(f"Written trimmed alignment: {args.output}")
    print(f"Total input: {len(aligned_records)}")
    print(f"Total written: {len(trimmed_records)}")
    print(f"Total fragmented discarded (coverage < {args.min_coverage}): {fragmented_count}")

if __name__ == "__main__":
    main()
