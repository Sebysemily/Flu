#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
from copy import deepcopy
from io import StringIO

from Bio import Phylo, SeqIO


def load_taxa(path: str):
    taxa = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxa.append(row["taxon"])
    return taxa


def read_tree(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        return Phylo.read(StringIO(handle.read().strip()), "newick")


def main():
    parser = argparse.ArgumentParser(description="Subset alignment and prune tree to panel taxa")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--taxa", required=True)
    parser.add_argument("--out-alignment", required=True)
    parser.add_argument("--out-tree", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()

    taxa = set(load_taxa(args.taxa))

    # Subset alignment
    records = list(SeqIO.parse(args.alignment, "fasta"))
    kept = [r for r in records if r.id in taxa]
    out_dir = os.path.dirname(args.out_alignment)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    tmp_subset = args.out_alignment + ".subset_input.fasta"
    SeqIO.write(kept, tmp_subset, "fasta")

    # Realign panel subset to avoid slicing bias from full concat alignment.
    mafft_cmd = ["mafft", "--auto", "--thread", "1", tmp_subset]
    with open(args.out_alignment, "w", encoding="utf-8") as handle:
        proc = subprocess.run(mafft_cmd, stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
    if proc.returncode != 0:
        raise SystemExit(f"MAFFT failed ({proc.returncode}): {proc.stderr}")
    try:
        os.remove(tmp_subset)
    except OSError:
        pass

    # Prune tree
    tree = read_tree(args.tree)
    pruned = deepcopy(tree)
    for tip in list(pruned.get_terminals()):
        if tip.name not in taxa:
            try:
                pruned.prune(tip)
            except ValueError:
                # Ignore edge-case if node already pruned through parent collapse
                pass
    Phylo.write(pruned, args.out_tree, "newick")

    tree_tips = {t.name for t in tree.get_terminals() if t.name}
    aln_ids = {r.id for r in records}
    missing_in_tree = sorted(taxa - tree_tips)
    missing_in_alignment = sorted(taxa - aln_ids)

    with open(args.audit, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"requested_taxa\t{len(taxa)}\n")
        handle.write(f"kept_alignment\t{len(kept)}\n")
        handle.write(f"kept_tree\t{len(pruned.get_terminals())}\n")
        handle.write(f"missing_in_tree\t{';'.join(missing_in_tree)}\n")
        handle.write(f"missing_in_alignment\t{';'.join(missing_in_alignment)}\n")


if __name__ == "__main__":
    main()
