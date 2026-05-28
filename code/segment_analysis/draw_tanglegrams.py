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


def collapse_maximal_context_clades(tree, core_prefix="Flu-"):
    """
    Finds and collapses maximal clades in-place that contain only context sequences.
    Returns a dictionary mapping the new leaf name to the set of original leaf names.
    """
    collapsed_groups = {}
    group_counter = 0

    def is_context_only(node):
        return all(not leaf.name.lower().startswith(core_prefix.lower()) for leaf in node.get_terminals())

    clades_to_collapse = []

    def find_maximal_context(node):
        if is_context_only(node):
            clades_to_collapse.append(node)
            return
        if node.is_terminal():
            return
        for child in node.clades:
            find_maximal_context(child)

    find_maximal_context(tree.root)

    # Collapse each identified clade
    for node in clades_to_collapse:
        terminals = node.get_terminals()
        leaf_names = {l.name for l in terminals}
        count = len(leaf_names)

        if count > 1:
            group_counter += 1
            new_name = f"Regional Context Group {group_counter} (N={count})"
            node.name = new_name
            node.clades = []  # Collapse descendants
            collapsed_groups[new_name] = leaf_names
        else:
            # If it's a single leaf, keep it as is
            leaf_name = list(leaf_names)[0]
            collapsed_groups[leaf_name] = leaf_names

    # Add uncollapsed core leaves to the mapping
    for leaf in tree.get_terminals():
        if leaf.name not in collapsed_groups:
            collapsed_groups[leaf.name] = {leaf.name}

    return collapsed_groups


def optimize_tree_layout(tree, reference_leaf_order, collapsed_members):
    """
    Applies the barycenter heuristic to rotate tree nodes in-place
    to match the reference leaf order as closely as possible.
    """
    memo = {}

    def get_descendant_indices(node):
        if node in memo:
            return memo[node]
        if node.is_terminal():
            # For collapsed groups, find the average index of its members in the reference order
            members = collapsed_members.get(node.name, {node.name})
            indices = [reference_leaf_order.get(m, 0) for m in members if m in reference_leaf_order]
            val = indices if indices else [0]
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
    
    depths = tree.depths()
    max_depth = max(depths.values()) if depths.values() else 1.0
    if max_depth == 0:
        max_depth = 1.0

    node_coords = {}

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
            x = (depth / max_depth) * 0.8
        else:
            x = 3.0 - (depth / max_depth) * 0.8
            
        node_coords[node] = (x, y)
        return x, y

    calc_coords(tree.root)
    return node_coords, y_coords


def load_itol_colors(style_file):
    color_map = {}
    if not style_file or not os.path.exists(style_file):
        return color_map
    with open(style_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("TREE_COLORS") or line.startswith("SEPARATOR") or line.startswith("DATA") or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) >= 3:
                name = parts[0].strip("'\"")
                color = parts[2].strip()
                color_map[name] = color
    return color_map


def get_node_color(node, color_map, collapsed_members, default_color='#555555'):
    """
    Determines the color for a node based on its descendants.
    """
    terminals = node.get_terminals()
    
    # Expand collapsed leaves to original names
    all_original_leaves = []
    for t in terminals:
        members = collapsed_members.get(t.name, {t.name})
        all_original_leaves.extend(members)
        
    colors = [color_map.get(leaf_name) for leaf_name in all_original_leaves if leaf_name in color_map]
    
    if not colors:
        return default_color
        
    # If all leaves share a single color, use it
    first_color = colors[0]
    if all(c == first_color for c in colors):
        return first_color
        
    return default_color


def draw_tree_branches(ax, tree, coords, color_map, collapsed_members):
    parent_dict = build_parent_dict(tree)
    for child, parent in parent_dict.items():
        x_p, y_p = coords[parent]
        x_c, y_c = coords[child]
        
        # Determine branch color
        branch_color = get_node_color(child, color_map, collapsed_members)
        
        # Draw vertical line at parent x
        ax.plot([x_p, x_p], [y_p, y_c], color=branch_color, linewidth=1.2)
        # Draw horizontal line to child
        ax.plot([x_p, x_c], [y_c, y_c], color=branch_color, linewidth=1.2)


def main():
    parser = argparse.ArgumentParser(description="Generate a beautiful tanglegram between two trees.")
    parser.add_argument("--tree1", required=True, help="Path to Left Newick Tree")
    parser.add_argument("--tree2", required=True, help="Path to Right Newick Tree")
    parser.add_argument("--label1", default="Tree 1", help="Label for Left Tree")
    parser.add_argument("--label2", default="Tree 2", help="Label for Right Tree")
    parser.add_argument("--itol-colors", help="Path to iTOL tree_colors.txt template file")
    parser.add_argument("--collapse-context", action="store_true", help="Collapse maximal context-only clades")
    parser.add_argument("--hide-non-core-labels", action="store_true", help="Hide labels for non-core sequences (except collapsed group titles)")
    parser.add_argument("--output", required=True, help="Path to output PNG image")
    args = parser.parse_args()

    # 1. Load trees
    t1 = Phylo.read(args.tree1, "newick")
    t2 = Phylo.read(args.tree2, "newick")

    # 2. Prune to common taxa
    common_taxa = prune_to_common(t1, t2)
    print(f"Loaded {len(common_taxa)} common taxa.")

    if not common_taxa:
        raise ValueError("No common taxa between the two trees. Cannot draw tanglegram.")

    # Load colors
    color_map = load_itol_colors(args.itol_colors)

    # 3. Collapse context clades if requested
    collapsed_members1 = {}
    collapsed_members2 = {}
    if args.collapse_context:
        collapsed_members1 = collapse_maximal_context_clades(t1)
        collapsed_members2 = collapse_maximal_context_clades(t2)
        print(f"Collapsed Tree 1 leaves count: {len(t1.get_terminals())}")
        print(f"Collapsed Tree 2 leaves count: {len(t2.get_terminals())}")
    else:
        # Default identity mapping
        collapsed_members1 = {l.name: {l.name} for l in t1.get_terminals()}
        collapsed_members2 = {l.name: {l.name} for l in t2.get_terminals()}

    # 4. Optimize layout of Tree 2 to match Tree 1
    leaves1_names = [l.name for l in t1.get_terminals()]
    # Reference order is based on original uncollapsed leaf order to maintain detailed alignment
    original_leaves_order = {name: idx for idx, name in enumerate(leaves1_names)}
    optimize_tree_layout(t2, original_leaves_order, collapsed_members2)

    # 5. Generate coordinates
    coords1, y_coords1 = get_coordinates(t1, side='left')
    coords2, y_coords2 = get_coordinates(t2, side='right')

    # 6. Set up plot
    num_taxa1 = len(t1.get_terminals())
    num_taxa2 = len(t2.get_terminals())
    max_taxa = max(num_taxa1, num_taxa2)
    
    fig_height = max(6, max_taxa * 0.25)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    
    ax.axis('off')
    ax.set_xlim(-0.1, 3.1)
    ax.set_ylim(-1.0, max_taxa + 1.0)

    # 7. Draw tree branches
    draw_tree_branches(ax, t1, coords1, color_map, collapsed_members1)
    draw_tree_branches(ax, t2, coords2, color_map, collapsed_members2)

    # 8. Draw labels and connection lines
    x_conn_left = 1.35
    x_conn_right = 1.65

    # A: Draw leaf nodes (dots and text)
    # Left Tree
    for leaf in t1.get_terminals():
        y = y_coords1[leaf.name]
        is_core = leaf.name.lower().startswith("flu-")
        color = get_node_color(leaf, color_map, collapsed_members1)
        
        # Plot point
        ax.plot(0.8, y, 'o', color=color, markersize=4)
        
        # Plot text
        show_label = is_core or not args.hide_non_core_labels or leaf.name.startswith("Regional Context Group")
        if show_label:
            label_text = leaf.name.split("|")[0] if "|" in leaf.name else leaf.name
            ax.text(0.82, y, label_text, ha='left', va='center', fontsize=8, color='#333333', family='monospace')

    # Right Tree
    for leaf in t2.get_terminals():
        y = y_coords2[leaf.name]
        is_core = leaf.name.lower().startswith("flu-")
        color = get_node_color(leaf, color_map, collapsed_members2)
        
        # Plot point
        ax.plot(2.2, y, 'o', color=color, markersize=4)
        
        # Plot text
        show_label = is_core or not args.hide_non_core_labels or leaf.name.startswith("Regional Context Group")
        if show_label:
            label_text = leaf.name.split("|")[0] if "|" in leaf.name else leaf.name
            ax.text(2.18, y, label_text, ha='right', va='center', fontsize=8, color='#333333', family='monospace')

    # B: Draw Connection Lines
    # We draw lines between any pair of collapsed leaves that share original leaves.
    for leaf1, members1 in collapsed_members1.items():
        y1 = y_coords1[leaf1]
        for leaf2, members2 in collapsed_members2.items():
            y2 = y_coords2[leaf2]
            
            shared = members1.intersection(members2)
            if not shared:
                continue
                
            # Line properties
            is_crossing = abs(y1 - y2) > 0.01
            
            # Determine color from shared leaves
            shared_list = list(shared)
            first_shared = shared_list[0]
            shared_color = color_map.get(first_shared, '#cccccc')
            
            if is_crossing:
                # Highlight crossings with a bit more alpha
                color = '#ff4d4d' if not args.itol_colors else shared_color
                alpha = 0.6
                linewidth = 1.2
                style = '-'
            else:
                color = '#cccccc' if not args.itol_colors else shared_color
                alpha = 0.3
                linewidth = 0.8
                style = '--'
                
            ax.plot([x_conn_left, x_conn_right], [y1, y2], color=color, alpha=alpha, linewidth=linewidth, linestyle=style)

    # 9. Add Titles
    ax.text(0.4, max_taxa, args.label1, ha='center', va='bottom', fontsize=12, fontweight='bold', color='#1a1a1a')
    ax.text(2.6, max_taxa, args.label2, ha='center', va='bottom', fontsize=12, fontweight='bold', color='#1a1a1a')
    ax.text(1.5, max_taxa + 0.5, f"Tanglegram: {args.label1} vs {args.label2}", ha='center', va='bottom', fontsize=14, fontweight='bold', color='#111111')

    # Save output
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Success! Tanglegram saved to {args.output}")


if __name__ == "__main__":
    main()
