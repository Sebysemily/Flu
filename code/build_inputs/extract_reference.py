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
    parser.add_argument("--output-gff", required=False, help="Path to write GFF3 genome annotation.")
    parser.add_argument("--output-json", required=False, help="Path to write pathogen.json configuration.")
    parser.add_argument("--segment", required=False, help="Segment name (e.g. PB2, HA, etc.).")
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

    # Write output FASTA
    with open(args.output_fasta, "w", encoding="utf-8") as out_fh:
        out_fh.write(f">{ref_header}\n")
        # wrap sequence
        width = 80
        wrapped_seq = "\n".join(ref_seq[i : i + width] for i in range(0, len(ref_seq), width))
        out_fh.write(wrapped_seq + "\n")

    # Extract clean accession/ID (first word of header)
    ref_header_id = ref_header.split()[0] if ref_header else "reference"

    # Write GFF3 and pathogen.json if paths are provided
    if args.output_gff and args.output_json and args.segment:
        # Find first in-frame stop codon
        cds_end = len(ref_seq)
        for i in range(0, len(ref_seq) - 2, 3):
            codon = ref_seq[i:i+3].upper()
            if codon in ["TAA", "TAG", "TGA"]:
                cds_end = i + 3
                break

        # Write GFF3
        with open(args.output_gff, "w", encoding="utf-8") as gff_fh:
            gff_fh.write("##gff-version 3\n")
            gff_fh.write(f"##sequence-region {ref_header_id} 1 {len(ref_seq)}\n")
            gff_fh.write(f"{ref_header_id}\tcustom\tCDS\t1\t{cds_end}\t.\t+\t0\tgene_name={args.segment};gene={args.segment}\n")

        # Write pathogen.json
        import json
        pat = {
            "$schema": "https://raw.githubusercontent.com/nextstrain/nextclade/refs/heads/release/packages/nextclade-schemas/input-pathogen-json.schema.json",
            "alignmentParams": {
                "excessBandwidth": 9,
                "terminalBandwidth": 100,
                "allowedMismatches": 4,
                "gapAlignmentSide": "right",
                "minSeedCover": 0.1
            },
            "attributes": {
                "name": f"Influenza A H5N1 {args.segment}",
                "reference accession": ref_header_id,
                "segment": args.segment
            },
            "defaultCds": args.segment,
            "qc": {
                "frameShifts": { "enabled": True },
                "missingData": { "enabled": False },
                "mixedSites": { "enabled": True },
                "privateMutations": { "enabled": False },
                "stopCodons": { "enabled": True }
            },
            "schemaVersion": "3.0.0"
        }
        with open(args.output_json, "w", encoding="utf-8") as json_fh:
            json.dump(pat, json_fh, indent=2)

if __name__ == "__main__":
    main()
