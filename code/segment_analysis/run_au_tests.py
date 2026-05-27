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


def parse_au_pvalues(iqtree_file):
    if not os.path.exists(iqtree_file):
        raise FileNotFoundError(f"IQ-TREE output file not found: {iqtree_file}")

    print(f"Parsing AU p-values from {iqtree_file}...")
    p_au_values = []
    in_table = False

    with open(iqtree_file, "r") as f:
        for line in f:
            line = line.strip()
            if "Tree" in line and "p-AU" in line:
                in_table = True
                continue
            if in_table:
                if not line or line.startswith("---") or line.startswith("="):
                    if p_au_values:
                        break  # End of table
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[0].isdigit():
                    # The last element might be a sign like '-' or '+' or '*'
                    val_str = parts[-1]
                    try:
                        p_au_val = float(val_str)
                    except ValueError:
                        p_au_val = float(parts[-2])
                    p_au_values.append(p_au_val)

    return p_au_values


def main():
    parser = argparse.ArgumentParser(description="Run AU topology tests on a matrix of alignments and trees.")
    parser.add_argument("--segments", required=True, help="Comma-separated list of segment names (e.g. PB2,PB1...)")
    parser.add_argument("--trees", required=True, nargs="+", help="Paths to final trees (in the order of segments)")
    parser.add_argument("--concat-tree", required=True, help="Path to final concatenated tree")
    parser.add_argument("--alignments-dir", required=True, help="Directory containing segment alignments")
    parser.add_argument("--concat-alignment", required=True, help="Path to final concatenated alignment")
    parser.add_argument("--concat-partitions", required=True, help="Path to final concatenated partitions")
    parser.add_argument("--work-dir", required=True, help="Directory to store intermediate outputs")
    parser.add_argument("--output", required=True, help="Path to output compiled CSV matrix")
    args = parser.parse_args()

    segments = args.segments.split(",")
    all_names = segments + ["concat"]
    
    # 1. Merge all trees into a single candidate topologies file
    candidate_trees_path = os.path.join(args.work_dir, "candidate_topologies.trees")
    all_tree_paths = args.trees + [args.concat_tree]
    merge_trees(all_tree_paths, candidate_trees_path)

    # 2. Run AU test for each alignment
    matrix_data = {}
    
    for i, aln_name in enumerate(all_names):
        print(f"\n--- Evaluating Alignment: {aln_name} ---")
        
        # Setup files for this alignment
        if aln_name == "concat":
            alignment_path = args.concat_alignment
            partition_path = args.concat_partitions
        else:
            alignment_path = os.path.join(args.alignments_dir, f"H5N1_{aln_name}.fasta")
            # If it's a codon segment, it has partitions
            partition_path = os.path.join(args.alignments_dir, f"H5N1_{aln_name}.partitions")
            if not os.path.exists(partition_path):
                partition_path = None
        
        if not os.path.exists(alignment_path):
            raise FileNotFoundError(f"Alignment file not found: {alignment_path}")
            
        prefix = os.path.join(args.work_dir, f"eval_{aln_name}")
        
        # Build command
        # Using -nt 2 to be fast but lightweight
        # -zb 10000 specifies the number of RELL bootstrap replicates required for the AU test
        cmd = [
            "iqtree",
            "-s", alignment_path,
            "-z", candidate_trees_path,
            "-au",
            "-zb", "10000",
            "-n", "0",
            "-pre", prefix,
            "-nt", "2"
        ]
        if partition_path and os.path.exists(partition_path):
            cmd.extend(["-spp", partition_path])
            
        print(f"Running IQ-TREE: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        # Parse results
        iqtree_out_path = f"{prefix}.iqtree"
        p_values = parse_au_pvalues(iqtree_out_path)
        
        if len(p_values) != len(all_names):
            raise ValueError(f"Expected {len(all_names)} p-values from {iqtree_out_path}, but parsed {len(p_values)}")
            
        matrix_data[aln_name] = p_values

    # 3. Create DataFrame and export
    # Rows = alignments, Columns = topologies
    df = pd.DataFrame.from_dict(matrix_data, orient='index', columns=all_names)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output)
    print(f"\nSuccess! AU test matrix saved to {args.output}")


if __name__ == "__main__":
    main()
