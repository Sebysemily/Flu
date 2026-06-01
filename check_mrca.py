import sys
from Bio import Phylo
import csv

tree_path = "results/phylogeny/iq-tree/HA/lineage/H5N1_HA_fast.treefile"
meta_path = "metadata/H5N1_context.csv"

# Load metadata
meta = {}
with open(meta_path, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row.get("file_name", "").strip()
        role = row.get("expected_role", "").strip()
        if name:
            meta[name] = role

try:
    tree = Phylo.read(tree_path, "newick")
except Exception as e:
    print(f"Error reading tree: {e}")
    sys.exit(1)

# Find flu tips
flu_tips = []
for tip in tree.get_terminals():
    if tip.name in meta and meta[tip.name].startswith("flu_"):
        flu_tips.append(tip)

print(f"Total flu_ tips in tree: {len(flu_tips)}")

if len(flu_tips) < 2:
    print("Not enough flu_ tips to find MRCA.")
    sys.exit(0)

mrca = tree.common_ancestor(flu_tips)

# Find all descendants of this MRCA
descendants = mrca.get_terminals()
print(f"Total descendant tips in MRCA clade: {len(descendants)}")

# Breakdown of descendants by expected role
roles_count = {}
for tip in descendants:
    role = meta.get(tip.name, "unknown")
    roles_count[role] = roles_count.get(role, 0) + 1

print("\nBreakdown of MRCA descendants by role:")
for r, c in sorted(roles_count.items()):
    print(f"{r}: {c}")

