#!/usr/bin/env python3
import argparse
import csv
import os

from Bio import SeqIO


def load_taxa(path):
    taxa = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxon = str(row.get("taxon", "")).strip()
            if taxon:
                taxa.append(taxon)
    return taxa


def main():
    parser = argparse.ArgumentParser(description="Subset an alignment FASTA to a taxon panel")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--taxa", required=True)
    parser.add_argument("--out-alignment", required=True)
    parser.add_argument("--audit", default=None)
    args = parser.parse_args()

    taxa = set(load_taxa(args.taxa))
    records = list(SeqIO.parse(args.alignment, "fasta"))
    kept = [record for record in records if record.id in taxa]

    out_dir = os.path.dirname(args.out_alignment)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    SeqIO.write(kept, args.out_alignment, "fasta")

    if args.audit:
        audit_dir = os.path.dirname(args.audit)
        if audit_dir:
            os.makedirs(audit_dir, exist_ok=True)
        aln_ids = {record.id for record in records}
        missing_in_alignment = sorted(taxa - aln_ids)
        with open(args.audit, "w", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write(f"requested_taxa\t{len(taxa)}\n")
            handle.write(f"kept_alignment\t{len(kept)}\n")
            handle.write(f"missing_in_alignment\t{';'.join(missing_in_alignment)}\n")


if __name__ == "__main__":
    main()
