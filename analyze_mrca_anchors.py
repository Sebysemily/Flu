import sys
import re
from Bio import Phylo
import csv

tree_path = "results/phylogeny/iq-tree/HA/lineage/H5N1_HA_fast.treefile"
meta_path = "metadata/H5N1_context.csv"

meta = {}
with open(meta_path, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row.get("file_name", "").strip()
        if name:
            meta[name] = {
                "role": row.get("expected_role", "").strip(),
                "country": row.get("country", "").strip(),
                "date": row.get("collection_date", "").strip()
            }

tree = Phylo.read(tree_path, "newick")

flu_tips = []
for tip in tree.get_terminals():
    if tip.name in meta and meta[tip.name]["role"].startswith("flu_"):
        flu_tips.append(tip)

mrca = tree.common_ancestor(flu_tips)
descendants = mrca.get_terminals()

american_anchors = [t.name for t in descendants if meta.get(t.name, {}).get("role") == "american_anchor"]

months = {}
countries = {}

for tip in american_anchors:
    m = meta.get(tip, {})
    # Parse month
    match = re.search(r"(\d{4})-(\d{2})", m.get("date", ""))
    month = f"{match.group(1)}-{match.group(2)}" if match else "UNKNOWN"
    months[month] = months.get(month, 0) + 1
    
    country = m.get("country", "Unknown")
    countries[country] = countries.get(country, 0) + 1

print("American Anchors in MRCA Clade by Month:")
for m in sorted(months.keys()):
    print(f"{m}: {months[m]}")

print("\nAmerican Anchors in MRCA Clade by Country:")
for c, count in sorted(countries.items()):
    print(f"{c}: {count}")

