#!/usr/bin/env python3
import argparse
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

def main():
    parser = argparse.ArgumentParser(description="Extract a reference sequence by ID from a FASTA file.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--output-fasta", required=True)
    args = parser.parse_args()

    ref_header = None
    ref_seq = None

    # Search for sequence containing reference-id
    for header, seq in read_fasta(args.input_fasta):
        if args.reference_id in header:
            ref_header = header
            ref_seq = seq
            break

    # Fallback to first sequence starting with 'Flu-' if not found
    if ref_seq is None:
        for header, seq in read_fasta(args.input_fasta):
            if header.startswith("Flu-"):
                ref_header = header
                ref_seq = seq
                break

    # Ultimate fallback to first sequence in file
    if ref_seq is None:
        for header, seq in read_fasta(args.input_fasta):
            ref_header = header
            ref_seq = seq
            break

    if ref_seq is None:
        print(f"Error: No sequences found in {args.input_fasta}")
        sys.exit(1)

    print(f"Extracted reference: {ref_header}")

    # Write output
    with open(args.output_fasta, "w", encoding="utf-8") as out_fh:
        out_fh.write(f">{ref_header}\n")
        # wrap sequence
        width = 80
        wrapped_seq = "\n".join(ref_seq[i : i + width] for i in range(0, len(ref_seq), width))
        out_fh.write(wrapped_seq + "\n")

if __name__ == "__main__":
    main()
