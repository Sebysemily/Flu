import sys

tree_count = 0
with open('results/beast/GSS/ucln_bay/H5N1_HA_panel_postQC.trees', 'r') as f:
    for line in f:
        if line.startswith('tree '):
            tree_count += 1
            
print(f"Total trees: {tree_count}")
