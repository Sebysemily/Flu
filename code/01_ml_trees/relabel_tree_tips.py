#!/usr/bin/env python3
"""
relabel_tree_tips.py
--------------------
Renames tip labels in a Newick tree from their full verbose form to the
short stable ID stored in metadata/main_panel.csv.

The metadata CSV must have (at minimum) two columns:
  label  – the full tip label as it appears in the input Newick
  id     – the short replacement label

Any tip not found in the metadata is kept as-is (with a warning).

Usage:
    python relabel_tree_tips.py \\
        --tree   results/phylogeny/iq-tree/HA/HA_final.treefile \\
        --metadata metadata/main_panel.csv \\
        --output results/phylogeny/iq-tree/HA/HA_viz.treefile
"""

import argparse
import csv
import os
import re
import sys


def load_label_map(metadata_path: str) -> dict:
    """Return {full_label: short_id} from the metadata CSV."""
    mapping = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            mapping[row["label"]] = row["id"]
    return mapping


def relabel_newick(newick_str: str, label_map: dict) -> tuple[str, int, int]:
    """
    Replace every tip label in the Newick string with its mapped short ID.

    Returns (relabeled_newick, n_replaced, n_missing).
    Tips with no mapping entry are kept unchanged and counted as missing.
    """
    replaced = 0
    missing = 0

    def replace_tip(m):
        nonlocal replaced, missing
        full = m.group(1)
        short = label_map.get(full)
        if short is None:
            missing += 1
            return m.group(0)          # keep original
        replaced += 1
        # Re-attach the branch-length / separator that followed
        return short + m.group(2)

    # Newick tips: appear after '(' or ',' and are followed by ':' or ')' or ','
    # Pattern captures: (tip_label)(trailing_separator)
    # We use a negative-lookbehind to avoid matching branch lengths (pure digits/dots)
    pattern = re.compile(r"(?<=[(,])([^():,;\s][^():,;]*)(:)")
    result = pattern.sub(replace_tip, newick_str)
    return result, replaced, missing


def main():
    parser = argparse.ArgumentParser(
        description="Relabel Newick tip labels using a metadata CSV."
    )
    parser.add_argument("--tree",     required=True, help="Input Newick tree file.")
    parser.add_argument("--metadata", required=True, help="metadata/main_panel.csv path.")
    parser.add_argument("--output",   required=True, help="Output Newick file path.")
    args = parser.parse_args()

    if not os.path.exists(args.tree):
        sys.exit(f"ERROR: tree not found: {args.tree}")
    if not os.path.exists(args.metadata):
        sys.exit(f"ERROR: metadata not found: {args.metadata}")

    label_map = load_label_map(args.metadata)

    with open(args.tree) as fh:
        newick = fh.read()

    relabeled, n_ok, n_miss = relabel_newick(newick, label_map)

    if n_miss:
        print(f"WARNING: {n_miss} tip(s) had no metadata entry and were kept as-is.",
              file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(relabeled)

    print(f"Relabeled {n_ok} tips → {args.output}  ({n_miss} kept unchanged)")


if __name__ == "__main__":
    main()
