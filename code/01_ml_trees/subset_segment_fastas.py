#!/usr/bin/env python3
import argparse
import csv
import os
import re
from typing import Dict, Set
from Bio import SeqIO

def load_context_taxa(path: str) -> Set[str]:
    """Load context taxa IDs from text file."""
    taxa = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                taxa.add(line)
    return taxa

def load_core_taxa_from_metadata(path: str) -> Set[str]:
    """Load core taxa (flu_*) from metadata."""
    core_taxa = set()
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("file_name") or "").strip()
            role = (row.get("expected_role") or "").strip()
            if name and role and re.match(r"^flu_.*", role, re.IGNORECASE):
                core_taxa.add(name)
    return core_taxa

def main():
    parser = argparse.ArgumentParser(description="Subset segment FASTA using context list and metadata core flu samples.")
    parser.add_argument("--unaligned-fasta", required=True, help="Input segment FASTA")
    parser.add_argument("--context-list", required=True, help="List of selected context taxa")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV")
    parser.add_argument("--out-fasta", required=True, help="Output pruned FASTA")
    args = parser.parse_args()

    # Load context and core taxa
    context_taxa = load_context_taxa(args.context_list)
    core_taxa = load_core_taxa_from_metadata(args.metadata)
    
    # Combined set of allowed taxa
    allowed_taxa = context_taxa.union(core_taxa)

    print(f"Loaded {len(context_taxa)} context taxa.")
    print(f"Loaded {len(core_taxa)} core flu taxa.")
    print(f"Total allowed taxa pool: {len(allowed_taxa)}.")

    # Read original FASTA, filter and strip gaps (if any)
    print(f"Filtering {args.unaligned_fasta}...")
    records = list(SeqIO.parse(args.unaligned_fasta, "fasta"))
    kept_records = []
    
    for r in records:
        if r.id in allowed_taxa:
            unaligned_seq = str(r.seq).replace("-", "").replace("~", "")
            r.seq = r.seq.__class__(unaligned_seq)
            kept_records.append(r)

    out_dir = os.path.dirname(args.out_fasta)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    SeqIO.write(kept_records, args.out_fasta, "fasta")
    print(f"Wrote {len(kept_records)} sequences to {args.out_fasta}.")

if __name__ == "__main__":
    main()
