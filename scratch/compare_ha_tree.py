#!/usr/bin/env python3
"""compare_ha_tree.py

Utility script to compare the set of EPI_ISL identifiers found in the HA
Newick tree file with those listed in the metadata CSV
metadata/H5N1_context.csv.

Usage::
    python compare_ha_tree.py <ha_tree_file> <metadata_csv>

The script prints two sections:
  * IDs present in the HA tree but missing from the CSV
  * IDs present in the CSV but missing from the HA tree

Both sets are written to ``missing_from_csv.txt`` and ``missing_from_tree.txt``
for downstream inspection.
"""

import sys
import re
import csv
from pathlib import Path

def extract_ids_from_tree(tree_path: Path) -> set:
    """Extract all EPI_ISL identifiers from a Newick tree file.
    The function searches for the pattern ``EPI_ISL_\d+`` anywhere in the file.
    """
    pattern = re.compile(r"EPI_ISL_\d+")
    ids = set()
    with tree_path.open("r", encoding="utf-8") as f:
        for line in f:
            ids.update(pattern.findall(line))
    return ids

def extract_ids_from_csv(csv_path: Path) -> set:
    """Extract the ``file_name`` column (which holds the EPI_ISL IDs) from the CSV."""
    ids = set()
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            eid = row.get("file_name", "").strip()
            if eid:
                ids.add(eid)
    return ids

def main(tree_file: str, csv_file: str) -> None:
    tree_path = Path(tree_file)
    csv_path = Path(csv_file)

    if not tree_path.is_file():
        sys.stderr.write(f"[ERROR] HA tree file not found: {tree_path}\n")
        sys.exit(1)
    if not csv_path.is_file():
        sys.stderr.write(f"[ERROR] CSV file not found: {csv_path}\n")
        sys.exit(1)

    tree_ids = extract_ids_from_tree(tree_path)
    csv_ids = extract_ids_from_csv(csv_path)

    missing_from_csv = sorted(tree_ids - csv_ids)
    missing_from_tree = sorted(csv_ids - tree_ids)

    print("=== EPI_ISL IDs in HA tree but NOT in metadata CSV ===")
    for eid in missing_from_csv:
        print(eid)
    print(f"Total missing from CSV: {len(missing_from_csv)}\n")

    print("=== EPI_ISL IDs in metadata CSV but NOT in HA tree ===")
    for eid in missing_from_tree:
        print(eid)
    print(f"Total missing from HA tree: {len(missing_from_tree)}\n")

    # Save results for convenience
    Path("missing_from_csv.txt").write_text("\n".join(missing_from_csv))
    Path("missing_from_tree.txt").write_text("\n".join(missing_from_tree))
    print("Results saved to missing_from_csv.txt and missing_from_tree.txt")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python compare_ha_tree.py <ha_tree_file> <metadata_csv>\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
