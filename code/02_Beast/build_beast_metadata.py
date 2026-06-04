#!/usr/bin/env python3
import argparse
import csv
import os
import sys

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build BEAST metadata TSV from H5N1_context.csv"
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata/H5N1_context.csv")
    parser.add_argument("--out", required=True, help="Output metadata TSV")
    parser.add_argument("--aln", required=False, help="Path to FASTA alignment to filter taxa")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    valid_taxa = None
    if args.aln:
        valid_taxa = set()
        with open(args.aln, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    valid_taxa.add(line[1:].strip())

    seen = set()
    rows = []

    with open(args.metadata, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            file_name = (row.get("file_name") or "").strip()
            date_val = (row.get("collection_date") or "").strip()
            country = (row.get("country") or "").strip()
            expected_role = (row.get("expected_role") or "").strip()

            if file_name and file_name not in seen:
                if valid_taxa is not None and file_name not in valid_taxa:
                    continue

                # Location mapping logic
                if country == "Ecuador":
                    if expected_role == "flu_costa":
                        location = "Ecuador_Coastal"
                    elif expected_role == "flu_sierra":
                        location = "Ecuador_Andine"
                    elif expected_role == "flu_amazonia":
                        location = "Ecuador_Amazon"
                    else:
                        location = "Ecuador"
                else:
                    location = country
                
                rows.append((file_name, date_val, location))
                seen.add(file_name)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("Taxon\tDate\tLocation\n")
        for file_name, date_val, location in sorted(rows):
            handle.write(f"{file_name}\t{date_val}\t{location}\n")

    print(f"Wrote {len(rows)} entries to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
