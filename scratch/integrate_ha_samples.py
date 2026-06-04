#!/usr/bin/env python3
"""integrate_ha_samples.py

Add HA‑only samples from ``data/context/HA_only_context.fasta`` to the
main context FASTA (``data/context/gisaid_epiflu_sequence.fasta``) and to the
metadata CSV (``metadata/H5N1_context.csv``).

The HA‑only FASTA header format is:
    >A/Host/Country/Accession|HA|YYYY-MM-DD|EPI_ISL_XXXXXX
We extract the fields and, if the EPI_ISL is not already present in the CSV,
append a new row with placeholder values for the missing columns.

Usage::
    python integrate_ha_samples.py

The script is idempotent – running it multiple times will not duplicate
entries.
"""

import csv
import os
from pathlib import Path
import re

# Paths (relative to the project root)
HA_FASTA = Path("data/context/HA_only_context.fasta")
MAIN_FASTA = Path("data/context/gisaid_epiflu_sequence.fasta")
CSV_FILE = Path("metadata/H5N1_context.csv")

# Expected CSV columns (from existing file header)
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

# Helper to parse a header line from the HA‑only FASTA
HEADER_RE = re.compile(r"^>([^|]+)\|HA\|([^|]+)\|([^|]+)$")

def parse_ha_header(header: str):
    """Return a dict with keys: file_name, host, country, collection_date, HA.
    Example header: >A/Avian/Chile/0114|HA|2023-03-21|EPI_ISL_20444471
    """
    m = HEADER_RE.match(header.strip())
    if not m:
        raise ValueError(f"Unexpected HA header format: {header}")
    # m.group(1) = A/Avian/Chile/0114
    parts = m.group(1).split('/')
    # Expected parts: ["A", host, country, accession]
    if len(parts) < 4:
        raise ValueError(f"Header does not contain enough parts: {header}")
    host = parts[1]
    country = parts[2]
    accession = parts[3]
    collection_date = m.group(2)
    epi_isl = m.group(3)
    return {
        "file_name": epi_isl,
        "host": host,
        "collection_date": collection_date,
        "country": country,
        "HA": epi_isl,  # in the original CSV the HA column holds the ID
    }

def load_existing_csv(csv_path: Path):
    existing = {}
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            existing[row["file_name"]] = row
    return existing

def load_existing_fasta_ids(fasta_path: Path):
    ids = set()
    with fasta_path.open("r") as f:
        for line in f:
            if line.startswith(">"):
                # The ID is the last field after the last pipe
                eid = line.strip().split("|")[-1]
                ids.add(eid)
    return ids

def main():
    if not HA_FASTA.is_file():
        raise FileNotFoundError(f"HA‑only FASTA not found: {HA_FASTA}")
    if not MAIN_FASTA.is_file():
        raise FileNotFoundError(f"Main context FASTA not found: {MAIN_FASTA}")
    if not CSV_FILE.is_file():
        raise FileNotFoundError(f"Metadata CSV not found: {CSV_FILE}")

    # Load current CSV content
    csv_data = load_existing_csv(CSV_FILE)
    existing_fasta_ids = load_existing_fasta_ids(MAIN_FASTA)

    # Prepare to collect new FASTA entries and CSV rows
    new_fasta_entries = []
    new_csv_rows = []

    with HA_FASTA.open("r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        header = lines[i].strip()
        seq = lines[i+1].strip() if i+1 < len(lines) else ""
        i += 2
        if not header.startswith(">"):
            continue
        info = parse_ha_header(header)
        eid = info["file_name"]
        # Append to FASTA if missing
        if eid not in existing_fasta_ids:
            new_fasta_entries.append(f"{header}\n{seq}\n")
        # Append to CSV if missing
        if eid not in csv_data:
            # Build a full row with placeholders for missing columns
            row = {col: "" for col in CSV_COLUMNS}
            row.update({
                "file_name": info["file_name"],
                "host": info["host"],
                "collection_date": info["collection_date"],
                "country": info["country"],
                "expected_role": "american_anchor",  # default placeholder
                "HA": info["HA"],
                "genotype": "A1",  # generic placeholder – adjust if known
            })
            new_csv_rows.append(row)

    # Write additions to FASTA (append mode)
    if new_fasta_entries:
        with MAIN_FASTA.open("a") as out_f:
            out_f.writelines(new_fasta_entries)
        print(f"Appended {len(new_fasta_entries)} sequences to {MAIN_FASTA}")
    else:
        print("No new FASTA entries needed.")

    # Write additions to CSV (append mode, preserving header order)
    if new_csv_rows:
        with CSV_FILE.open("a", newline="") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=CSV_COLUMNS)
            for row in new_csv_rows:
                writer.writerow(row)
        print(f"Appended {len(new_csv_rows)} rows to {CSV_FILE}")
    else:
        print("No new CSV rows needed.")

if __name__ == "__main__":
    main()
