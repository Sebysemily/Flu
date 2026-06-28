#!/usr/bin/env python3
import os
import re
import csv
import subprocess
import pandas as pd

def parse_header(header):
    # Example: A/turkey/Minnesota/22-012123-001/2022|EPI_ISL_16171254|2022-04-19|NP
    parts = header.split('|')
    strain_name = parts[0]
    epi_isl = ""
    date = ""
    for p in parts[1:]:
        if p.startswith("EPI_ISL_"):
            epi_isl = p
        elif re.match(r"^\d{4}(-\d{2}(-\d{2})?)?$", p):
            date = p
    return strain_name, epi_isl, date

def get_local_epis(csv_path):
    epis = set()
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path, dtype=str)
        if "EPI_ISL" in df.columns:
            epis = set(df["EPI_ISL"].dropna().str.strip())
    return epis

def filter_american_anchor(fasta_path, local_epis):
    # Create rep fasta and metadata
    rep_fasta = "data/context/temp_rep.fasta"
    meta_tsv = "data/context/temp_meta.tsv"
    
    records = {}
    with open(fasta_path, "r") as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    _, epi_isl, date = parse_header(header[1:])
                    if epi_isl and epi_isl not in local_epis:
                        if epi_isl not in records:
                            records[epi_isl] = {"date": date, "seq": "".join(seq), "full_headers": []}
                        records[epi_isl]["full_headers"].append((header, "".join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            _, epi_isl, date = parse_header(header[1:])
            if epi_isl and epi_isl not in local_epis:
                if epi_isl not in records:
                    records[epi_isl] = {"date": date, "seq": "".join(seq), "full_headers": []}
                records[epi_isl]["full_headers"].append((header, "".join(seq)))

    with open(rep_fasta, "w") as f_fa, open(meta_tsv, "w", newline="") as f_meta:
        writer = csv.writer(f_meta, delimiter="\t")
        writer.writerow(["strain", "date", "year_month"])
        for epi_isl, data in records.items():
            f_fa.write(f">{epi_isl}\n{data['seq']}\n")
            date_val = data["date"]
            if not date_val:
                date_val = "2022-01-01" # Dummy to avoid augur dropping it
            elif len(date_val) == 4:
                date_val += "-01-01"
            elif len(date_val) == 7:
                date_val += "-01"
            
            # create year_month
            parts = date_val.split("-")
            ym = f"{parts[0]}-{parts[1]}" if len(parts) >= 2 else "2022-01"
            
            writer.writerow([epi_isl, date_val, ym])
            
    # Run augur filter
    out_strains = "data/context/temp_strains.txt"
    cmd = [
        "augur", "filter",
        "--sequences", rep_fasta,
        "--metadata", meta_tsv,
        "--group-by", "year_month",
        "--sequences-per-group", "100",
        "--subsample-seed", "42",
        "--output-strains", out_strains
    ]
    subprocess.run(cmd, check=True)
    
    keep_epis = set()
    with open(out_strains, "r") as f:
        for line in f:
            keep_epis.add(line.strip())
            
    final_records = []
    for epi_isl in keep_epis:
        if epi_isl in records:
            final_records.extend(records[epi_isl]["full_headers"])
            
    # cleanup
    os.remove(rep_fasta)
    os.remove(meta_tsv)
    os.remove(out_strains)
    
    return final_records

def process_regional(full_fasta, ha_fasta, local_epis):
    full_epis = set()
    full_records = []
    
    # Process Full
    with open(full_fasta, "r") as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    _, epi_isl, _ = parse_header(header[1:])
                    if epi_isl and epi_isl not in local_epis:
                        full_epis.add(epi_isl)
                        full_records.append((header, "".join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            _, epi_isl, _ = parse_header(header[1:])
            if epi_isl and epi_isl not in local_epis:
                full_epis.add(epi_isl)
                full_records.append((header, "".join(seq)))
                
    # Process HA only
    ha_records = []
    with open(ha_fasta, "r") as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    _, epi_isl, _ = parse_header(header[1:])
                    if epi_isl and epi_isl not in local_epis and epi_isl not in full_epis:
                        ha_records.append((header, "".join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            _, epi_isl, _ = parse_header(header[1:])
            if epi_isl and epi_isl not in local_epis and epi_isl not in full_epis:
                ha_records.append((header, "".join(seq)))
                
    return full_records, ha_records

if __name__ == "__main__":
    local_epis = get_local_epis("metadata/flu_filtrado.csv")
    print(f"Loaded {len(local_epis)} local EPI_ISLs to exclude.")
    
    print("Filtering American Anchor...")
    anchor_records = filter_american_anchor("data/context/american_anchor_12-31-21_09-30-22.fasta", local_epis)
    
    print("Processing Regional Contexts...")
    full_records, ha_records = process_regional(
        "data/context/regional_context_full.fasta",
        "data/context/HA_only_regional_context.fasta",
        local_epis
    )
    
    out_path = "data/context/gisaid_epiflu_sequence.fasta"
    print(f"Writing to {out_path}...")
    with open(out_path, "w") as f:
        for header, seq in anchor_records:
            f.write(f"{header}\n{seq}\n")
        for header, seq in full_records:
            f.write(f"{header}\n{seq}\n")
        for header, seq in ha_records:
            f.write(f"{header}\n{seq}\n")
            
    print(f"Done! Anchor: {len(anchor_records)} segments, Regional Full: {len(full_records)} segments, Regional HA: {len(ha_records)} segments.")
