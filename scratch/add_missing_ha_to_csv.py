#!/usr/bin/env python3
"""add_missing_ha_to_csv.py

Extracts EPI_ISL identifiers from a HA Newick tree file and ensures they are
present in ``metadata/H5N1_context.csv``. If an ID is missing, a placeholder row
is appended.

Placeholder fields (adjust as needed):
  - host            = "unknown"
  - collection_date = "" (empty)
  - country         = "unknown"
  - expected_role   = "american_anchor"
  - genotype        = "A1"
  - all segment columns (PB2, PB1, PA, HA, NP, NA, MP, NS) = ""

Usage::
    python add_missing_ha_to_csv.py <ha_tree_file> <metadata_csv>

The script is idempotent – running it again will not duplicate rows.
"""

import csv
import re
import sys
from pathlib import Path

# Expected CSV columns (taken from the existing H5N1_context.csv header)
CSV_COLUMNS = [
    "file_name",
    "host",
    "collection_date",
    "country",
    "expected_role",
    "PB2",
    "PB1",
    "PA",
    "HA",
    "NP",
    "NA",
    "MP",
    "NS",
    "genotype",
    "PB2_lineage",
    "PB1_lineage",
    "PA_lineage",
    "HA_lineage",
    "NP_lineage",
    "NA_lineage",
    "MP_lineage",
    "NS_lineage",
]

EPI_REGEX = re.compile(r"EPI_ISL_\d+")

def extract_ids_from_tree(tree_path: Path) -> set:
    ids = set()
    with tree_path.open("r", encoding="utf-8") as f:
        for line in f:
            ids.update(EPI_REGEX.findall(line))
    return ids

def load_csv_ids(csv_path: Path) -> set:
    ids = set()
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.add(row["file_name"].strip())
    return ids

def append_missing_rows(csv_path: Path, missing_ids: set):
    if not missing_ids:
        print("No missing IDs – CSV already up‑to‑date.")
        return
    with csv_path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        for eid in sorted(missing_ids):
            row = {col: "" for col in CSV_COLUMNS}
            row.update({
                "file_name": eid,
                "host": "unknown",
                "collection_date": "",
                "country": "unknown",
                "expected_role": "american_anchor",
                "HA": eid,
                "genotype": "A1",
            })
            writer.writerow(row)
    print(f"Appended {len(missing_ids)} rows to {csv_path}")

def main(tree_file: str, csv_file: str):
    tree_path = Path(tree_file)
    csv_path = Path(csv_file)
    if not tree_path.is_file():
        sys.stderr.write(f"[ERROR] Tree file not found: {tree_path}\n")
        sys.exit(1)
    if not csv_path.is_file():
        sys.stderr.write(f"[ERROR] CSV file not found: {csv_path}\n")
        sys.exit(1)

    tree_ids = extract_ids_from_tree(tree_path)
    csv_ids = load_csv_ids(csv_path)
    missing = tree_ids - csv_ids
    print(f"Total IDs in tree: {len(tree_ids)}")
    print(f"IDs already in CSV: {len(csv_ids)}")
    print(f"Missing from CSV: {len(missing)}")
    append_missing_rows(csv_path, missing)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python add_missing_ha_to_csv.py <ha_tree_file> <metadata_csv>\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
