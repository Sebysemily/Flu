#!/usr/bin/env python3
import argparse
import os
import subprocess
import pandas as pd


def merge_trees(tree_paths, output_path):
    print(f"Merging {len(tree_paths)} trees into {output_path}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as out_f:
        for p in tree_paths:
            if not os.path.exists(p):
                raise FileNotFoundError(f"Tree file not found: {p}")
            with open(p, "r") as in_f:
                tree_str = in_f.read().strip()
                # Ensure it ends with a semicolon
                if tree_str and not tree_str.endswith(";"):
                    tree_str += ";"
                out_f.write(tree_str + "\n")


def parse_rfdist(rfdist_file, labels, rf_max):
    if not os.path.exists(rfdist_file):
        raise FileNotFoundError(f"rfdist file not found: {rfdist_file}")

    print(f"Parsing RF distances from {rfdist_file}...")
    with open(rfdist_file, "r") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"Empty rfdist file: {rfdist_file}")

    # First line contains size, e.g. "9 9"
    size_parts = lines[0].strip().split()
    n_trees = int(size_parts[0])

    if n_trees != len(labels):
        raise ValueError(f"Number of trees in rfdist ({n_trees}) does not match number of labels ({len(labels)})")

    matrix = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        # parts[0] is 'TreeX', parts[1:] are the distance values
        dists = [float(x) / rf_max for x in parts[1:]]  # Scale RF distances to 0-1
        matrix.append(dists)

    df = pd.DataFrame(matrix, index=labels, columns=labels)
    return df


def main():
    parser = argparse.ArgumentParser(description="Calculate all-to-all RF distance matrix between segment and concat trees.")
    parser.add_argument("--segments", required=True, help="Comma-separated list of segment names (e.g. PB2,PB1...)")
    parser.add_argument("--trees", required=True, nargs="+", help="Paths to final trees (in the order of segments)")
    parser.add_argument("--concat-tree", required=True, help="Path to final concatenated tree")
    parser.add_argument("--work-dir", required=True, help="Directory to store intermediate outputs")
    parser.add_argument("--output", required=True, help="Path to output compiled CSV matrix")
    args = parser.parse_args()

    segments = args.segments.split(",")
    labels = segments + ["concat"]

    # 1. Merge all trees into a single trees file
    trees_file_path = os.path.join(args.work_dir, "all_final_trees.trees")
    all_tree_paths = args.trees + [args.concat_tree]
    merge_trees(all_tree_paths, trees_file_path)

    # 2. Run iqtree -rf_all
    prefix = os.path.join(args.work_dir, "rf_all_matrix")
    
    # We remove previous files to avoid issues
    rfdist_file = f"{prefix}.rfdist"
    if os.path.exists(rfdist_file):
        os.remove(rfdist_file)

    cmd = [
        "iqtree",
        "-rf_all", trees_file_path,
        "-pre", prefix
    ]
    
    print(f"Running IQ-TREE: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # 3. Parse rfdist file and normalize
    from Bio import Phylo
    first_tree = Phylo.read(all_tree_paths[0], "newick")
    N = len(first_tree.get_terminals())
    rf_max = 2 * (N - 3)
    print(f"Number of leaves N = {N}, RF max = {rf_max}")

    df = parse_rfdist(rfdist_file, labels, rf_max)

    # 4. Save labeled CSV matrix
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output)
    print(f"Success! RF distance matrix saved to {args.output}")


if __name__ == "__main__":
    main()
