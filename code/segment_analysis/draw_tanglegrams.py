#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
from Bio import Phylo


def prune_to_common(tree1, tree2):
    leaves1 = {l.name for l in tree1.get_terminals()}
    leaves2 = {l.name for l in tree2.get_terminals()}
    common = leaves1.intersection(leaves2)

    # Prune tree1
    for l in list(tree1.get_terminals()):
        if l.name not in common:
            tree1.prune(l)

    # Prune tree2
    for l in list(tree2.get_terminals()):
        if l.name not in common:
            tree2.prune(l)

    return common


def optimize_tree_layout(tree, reference_leaf_order):
    """
    Applies the barycenter heuristic to rotate tree nodes in-place
    to match the reference leaf order as closely as possible.
    """
    memo = {}

    def get_descendant_indices(node):
        if node in memo:
            return memo[node]
        if node.is_terminal():
            val = [reference_leaf_order.get(node.name, 0)]
            memo[node] = val
            return val
        indices = []
        for child in node.clades:
            indices.extend(get_descendant_indices(child))
        memo[node] = indices
        return indices

    def sort_clades(node):
        if node.is_terminal():
            return
        for child in node.clades:
            sort_clades(child)
        
        # Sort children based on the mean index of their descendant leaves in Tree A
        node.clades.sort(key=lambda c: sum(get_descendant_indices(c)) / len(get_descendant_indices(c)))

    sort_clades(tree.root)


def build_parent_dict(tree):
    parent_dict = {}
    for p in tree.find_clades():
        for c in p.clades:
            parent_dict[c] = p
    return parent_dict


def get_coordinates(tree, side='left'):
    leaves = tree.get_terminals()
    y_coords = {leaf.name: i for i, leaf in enumerate(leaves)}
    
    # Compute depths (distances from root)
    depths = tree.depths()
    max_depth = max(depths.values()) if depths.values() else 1.0
    if max_depth == 0:
        max_depth = 1.0

    node_coords = {}
    parent_dict = build_parent_dict(tree)

    def calc_coords(node):
        if node.is_terminal():
            y = y_coords[node.name]
            x = 0.8 if side == 'left' else 2.2
            node_coords[node] = (x, y)
            return x, y

        child_coords = [calc_coords(child) for child in node.clades]
        y = sum(cy for cx, cy in child_coords) / len(child_coords)
        
        depth = depths[node]
        if side == 'left':
            # scale X from 0.0 (root) to 0.8 (leaves)
            x = (depth / max_depth) * 0.8
        else:
            # scale X from 3.0 (root) to 2.2 (leaves)
            x = 3.0 - (depth / max_depth) * 0.8
            
        node_coords[node] = (x, y)
        return x, y

    calc_coords(tree.root)
    return node_coords, y_coords


def draw_tree_branches(ax, tree, coords):
    parent_dict = build_parent_dict(tree)
    for child, parent in parent_dict.items():
        x_p, y_p = coords[parent]
        x_c, y_c = coords[child]
        
        # Draw vertical line at parent x
        ax.plot([x_p, x_p], [y_p, y_c], color='black', linewidth=1.2)
        # Draw horizontal line to child
        ax.plot([x_p, x_c], [y_c, y_c], color='black', linewidth=1.2)


def main():
    parser = argparse.ArgumentParser(description="Generate a beautiful tanglegram between two trees.")
    parser.add_argument("--tree1", required=True, help="Path to Left Newick Tree")
    parser.add_argument("--tree2", required=True, help="Path to Right Newick Tree")
    parser.add_argument("--label1", default="Tree 1", help="Label for Left Tree")
    parser.add_argument("--label2", default="Tree 2", help="Label for Right Tree")
    parser.add_argument("--output", required=True, help="Path to output PNG image")
    args = parser.parse_args()

    # 1. Load trees
    t1 = Phylo.read(args.tree1, "newick")
    t2 = Phylo.read(args.tree2, "newick")

    # 2. Prune to common taxa
    common_taxa = prune_to_common(t1, t2)
    print(f"Drawing tanglegram for {len(common_taxa)} common taxa.")

    if not common_taxa:
        raise ValueError("No common taxa between the two trees. Cannot draw tanglegram.")

    # 3. Optimize layout of Tree 2 to match Tree 1
    leaves1_names = [l.name for l in t1.get_terminals()]
    ref_order = {name: idx for idx, name in enumerate(leaves1_names)}
    optimize_tree_layout(t2, ref_order)

    # 4. Generate coordinates
    coords1, y_coords1 = get_coordinates(t1, side='left')
    coords2, y_coords2 = get_coordinates(t2, side='right')

    # 5. Set up plot
    num_taxa = len(common_taxa)
    fig_height = max(6, num_taxa * 0.25)  # Scale height based on number of taxa
    fig, ax = plt.subplots(figsize=(14, fig_height))
    
    # Hide axes
    ax.axis('off')
    ax.set_xlim(-0.1, 3.1)
    ax.set_ylim(-1.0, num_taxa + 1.0)

    # 6. Draw tree branches
    draw_tree_branches(ax, t1, coords1)
    draw_tree_branches(ax, t2, coords2)

    # 7. Draw labels and connection lines
    # Define connector boundaries
    x_conn_left = 1.35
    x_conn_right = 1.65
    
    for taxon in common_taxa:
        y1 = y_coords1[taxon]
        y2 = y_coords2[taxon]
        
        # Draw leaf labels
        # Left label (left-aligned at X=0.82)
        # We replace underscores/pipes with cleaner text if desired, or leave as is
        label_text = taxon.split("|")[0] if "|" in taxon else taxon
        ax.text(0.82, y1, label_text, ha='left', va='center', fontsize=8, color='#333333', family='monospace')
        # Right label (right-aligned at X=2.18)
        ax.text(2.18, y2, label_text, ha='right', va='center', fontsize=8, color='#333333', family='monospace')
        
        # Connection lines
        # If they match exactly (no crossing), draw a light gray line
        # If they cross, draw a vibrant tomato line to highlight reassortment
        is_crossing = abs(y1 - y2) > 0.01
        if is_crossing:
            color = '#ff4d4d'  # tomato
            alpha = 0.5
            linewidth = 1.0
            style = '-'
        else:
            color = '#cccccc'  # light gray
            alpha = 0.4
            linewidth = 0.8
            style = '--'
            
        ax.plot([x_conn_left, x_conn_right], [y1, y2], color=color, alpha=alpha, linewidth=linewidth, linestyle=style)

    # 8. Add Titles
    ax.text(0.4, num_taxa, args.label1, ha='center', va='bottom', fontsize=12, fontweight='bold', color='#1a1a1a')
    ax.text(2.6, num_taxa, args.label2, ha='center', va='bottom', fontsize=12, fontweight='bold', color='#1a1a1a')
    ax.text(1.5, num_taxa + 0.5, f"Tanglegram: {args.label1} vs {args.label2}", ha='center', va='bottom', fontsize=14, fontweight='bold', color='#111111')

    # Save output
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Success! Tanglegram saved to {args.output}")


if __name__ == "__main__":
    main()
