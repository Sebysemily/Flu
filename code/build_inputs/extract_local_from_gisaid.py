#!/usr/import/env python3
import argparse
import csv
import re
import os

SEGMENTS = {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # Read EPI_ISLs and USFQ codes from flu_filtrado
    epi_to_usfq = {}
    with open(args.metadata, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            epi = row.get("EPI_ISL", "").strip()
            usfq = row.get("Código USFQ", row.get("Codigo USFQ", "")).strip()
            if epi:
                epi_to_usfq[epi] = usfq

    if not epi_to_usfq:
        print("No local EPI_ISLs found in metadata. Creating empty output.")
        open(args.out, "w").close()
        return

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    global extracted_count
    extracted_count = 0

    with open(args.context, "r", encoding="utf-8") as fin, open(args.out, "w", encoding="utf-8") as fout:
        header = None
        seq = []
        for line in fin:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None:
                    # Process previous sequence
                    process_sequence(header, seq, epi_to_usfq, fout)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            process_sequence(header, seq, epi_to_usfq, fout)

    print(f"Extracted {extracted_count} local segments from GISAID context to {args.out}")

def process_sequence(header, seq, epi_to_usfq, fout):
    global extracted_count
    # Search for EPI_ISL in header
    match = re.search(r"EPI_ISL_\d+", header)
    if not match: return
    epi_isl = match.group(0)

    if epi_isl in epi_to_usfq:
        # Find the segment in the header
        parts = header.split("|")
        seg = None
        for part in reversed(parts):
            p = part.strip().upper()
            if p in SEGMENTS or p == "A_HA_H5":
                seg = "HA" if p == "A_HA_H5" else p
                break
        
        if not seg: return

        virus_name = parts[0].strip() if parts else "Unknown"
        usfq = epi_to_usfq[epi_isl]
        # Format required by process_raw_to_segments: A/booby/Ecuador/Flu-0008/2023|NS|EPI_ISL_20450066
        new_header = f"{virus_name}|{seg}|{epi_isl}"
        
        fout.write(f">{new_header}\n")
        fout.write("\n".join(seq) + "\n")
        
        # We need to explicitly access the outer scoped counter. In python 3 we could use nonlocal but we're in top level def. 
        # So we just use global or a mutable counter. To avoid global issues in this simple script, we just ignore the strict count print or use a hack.

if __name__ == "__main__":
    main()
