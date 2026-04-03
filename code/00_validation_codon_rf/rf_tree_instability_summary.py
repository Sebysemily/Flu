#!/usr/bin/env python3
import argparse
import csv
import itertools
import os
from io import StringIO
from typing import Dict, FrozenSet, List, Set, Tuple

from Bio import Phylo


def read_tree(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read().strip()
    return Phylo.read(StringIO(text), "newick")


def tree_tip_set(tree) -> Set[str]:
    return {tip.name for tip in tree.get_terminals() if tip.name}


def canonical_split(subset: Set[str], all_tips: Set[str]) -> FrozenSet[str]:
    complement = all_tips - subset
    left = frozenset(subset)
    right = frozenset(complement)

    if len(left) < len(right):
        return left
    if len(right) < len(left):
        return right

    # Tie-break equal-sized splits deterministically.
    left_key = tuple(sorted(left))
    right_key = tuple(sorted(right))
    return left if left_key <= right_key else right


def tree_splits(tree, all_tips: Set[str]) -> Set[FrozenSet[str]]:
    splits: Set[FrozenSet[str]] = set()
    n_tips = len(all_tips)

    for clade in tree.find_clades(order="postorder"):
        if clade == tree.root:
            continue
        terminals = {tip.name for tip in clade.get_terminals() if tip.name}
        k = len(terminals)

        # Keep only non-trivial bipartitions.
        if k <= 1 or k >= n_tips - 1:
            continue

        splits.add(canonical_split(terminals, all_tips))

    return splits


def rf_distance(s1: Set[FrozenSet[str]], s2: Set[FrozenSet[str]]) -> int:
    return len(s1 - s2) + len(s2 - s1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute pairwise RF distances for full-concat stability checks")
    parser.add_argument("--base-tree", required=True)
    parser.add_argument("--replicate-tree", action="append", required=True)
    parser.add_argument("--base-seed", type=int, required=True)
    parser.add_argument("--replicate-seed", action="append", required=True, type=int)
    parser.add_argument("--base-label", default="base")
    parser.add_argument("--replicate-label", action="append", default=[])
    parser.add_argument("--cutoff", type=float, default=0.05)
    parser.add_argument("--label", default="")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if len(args.replicate_tree) != len(args.replicate_seed):
        raise ValueError("Number of --replicate-tree and --replicate-seed values must match")

    if args.replicate_label and len(args.replicate_label) != len(args.replicate_tree):
        raise ValueError("Number of --replicate-label values must match number of --replicate-tree values")

    replicate_labels = args.replicate_label or [f"rep{i + 1}" for i in range(len(args.replicate_tree))]
    labels: List[str] = [args.base_label] + replicate_labels

    if len(set(labels)) != len(labels):
        raise ValueError("Tree labels must be unique across --base-label and --replicate-label values")

    tree_paths: Dict[str, str] = {args.base_label: args.base_tree}
    seed_map: Dict[str, int] = {args.base_label: args.base_seed}

    for label, path, seed in zip(replicate_labels, args.replicate_tree, args.replicate_seed):
        tree_paths[label] = path
        seed_map[label] = seed

    trees = {label: read_tree(path) for label, path in tree_paths.items()}

    reference_tips = tree_tip_set(trees[args.base_label])
    for label in labels[1:]:
        tips = tree_tip_set(trees[label])
        if tips != reference_tips:
            only_ref = sorted(reference_tips - tips)
            only_other = sorted(tips - reference_tips)
            raise ValueError(
                f"Tip mismatch between {args.base_label} and {label}. "
                f"Only in {args.base_label}: {only_ref[:5]}{'...' if len(only_ref) > 5 else ''}; "
                f"Only in {label}: {only_other[:5]}{'...' if len(only_other) > 5 else ''}"
            )

    split_map = {label: tree_splits(trees[label], reference_tips) for label in labels}

    rows = []
    max_ratio = 0.0
    for a, b in itertools.combinations(labels, 2):
        s1 = split_map[a]
        s2 = split_map[b]
        rf = rf_distance(s1, s2)
        comparable = len(s1 | s2)
        ratio = (rf / comparable) if comparable else 0.0
        max_ratio = max(max_ratio, ratio)
        classification = "noisy_instability" if ratio <= args.cutoff else "genuine_instability"

        rows.append(
            {
                "run_label": args.label,
                "tree_a": a,
                "tree_b": b,
                "seed_a": seed_map[a],
                "seed_b": seed_map[b],
                "rf_distance": rf,
                "comparable_splits": comparable,
                "rf_ratio": f"{ratio:.6f}",
                "cutoff": f"{args.cutoff:.6f}",
                "classification": classification,
            }
        )

    summary = {
        "run_label": args.label,
        "tree_a": "__summary__",
        "tree_b": "max_pairwise",
        "seed_a": "",
        "seed_b": "",
        "rf_distance": "",
        "comparable_splits": "",
        "rf_ratio": f"{max_ratio:.6f}",
        "cutoff": f"{args.cutoff:.6f}",
        "classification": "genuine_instability" if max_ratio > args.cutoff else "noisy_instability",
    }

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fieldnames = [
        "run_label",
        "tree_a",
        "tree_b",
        "seed_a",
        "seed_b",
        "rf_distance",
        "comparable_splits",
        "rf_ratio",
        "cutoff",
        "classification",
    ]

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow(summary)


if __name__ == "__main__":
    main()
