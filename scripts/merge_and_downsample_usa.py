#!/usr/bin/env python3
import sys
import os
import re
import random
import pandas as pd

# Set seed for reproducibility
random.seed(42)

# Add code/ to path to import local modules
sys.path.insert(0, os.path.abspath("code"))
from date_normalization import parse_collection_date
from build_inputs.process_raw_to_segments import (
    parse_context_isolates,
    filter_complete_context_isolates
)

def main():
    fasta_path = "data/context/gisaid_epiflu_sequence.fasta"
    extra_path = "data/context/usa_extra.fasta"
    filtrado_path = "metadata/flu_filtrado.csv"

    if not os.path.exists(fasta_path):
        print(f"Error: {fasta_path} not found.")
        sys.exit(1)
    if not os.path.exists(extra_path):
        print(f"Error: {extra_path} not found.")
        sys.exit(1)
    if not os.path.exists(filtrado_path):
        print(f"Error: {filtrado_path} not found.")
        sys.exit(1)

    # 1. Load local EPI_ISLs from flu_filtrado.csv
    df_filt = pd.read_csv(filtrado_path, dtype=str)
    local_epi_isls = set()
    if "EPI_ISL" in df_filt.columns:
        local_epi_isls = set(df_filt["EPI_ISL"].dropna().str.strip())
    print(f"Loaded {len(local_epi_isls)} local EPI_ISLs from {filtrado_path}")

    # 2. Parse both FASTA files
    print(f"Parsing {fasta_path}...")
    gisaid_isolates = parse_context_isolates(fasta_path)
    print(f"Loaded {len(gisaid_isolates)} isolates from main FASTA.")

    print(f"Parsing {extra_path}...")
    extra_isolates = parse_context_isolates(extra_path)
    print(f"Loaded {len(extra_isolates)} isolates from extra FASTA.")

    # 3. Merge them by isolate
    merged_isolates = {}
    merged_isolates.update(gisaid_isolates)
    
    added_count = 0
    for iso, segs in extra_isolates.items():
        if iso not in merged_isolates:
            merged_isolates[iso] = segs
            added_count += 1
        else:
            # If already exists, we can merge or keep gisaid
            # Let's verify we have all 8 segments
            existing = merged_isolates[iso]
            for seg, val in segs.items():
                if seg not in existing:
                    existing[seg] = val
                    
    print(f"Merged: added {added_count} new isolates from extra FASTA. Total isolates = {len(merged_isolates)}")

    # 4. Filter complete context isolates (len(segs) == 8, not local, valid date)
    complete_context, context_dates, context_places, context_types, context_provinces = filter_complete_context_isolates(
        merged_isolates, local_epi_isls
    )
    print(f"Total complete, non-local context isolates: {len(complete_context)}")

    # 5. Separate USA isolates and other isolates (excluding Canada)
    usa_isolates = {}
    other_isolates = {}
    canada_skipped = 0

    for iso, segs in complete_context.items():
        place = context_places[iso]
        place_clean = place.lower().replace("_", "").replace(" ", "")
        canada_set = {"nl", "canada", "newfoundland", "labrador", "ontario", "quebec", "alberta", "britishcolumbia", "manitoba", "saskatchewan", "novascotia", "newbrunswick", "pei"}
        if place_clean in canada_set:
            canada_skipped += 1
            continue
            
        if place == "USA":
            usa_isolates[iso] = segs
        else:
            other_isolates[iso] = segs

    print(f"Skipped {canada_skipped} Canadian context isolates.")

    print(f"Total USA complete isolates: {len(usa_isolates)}")
    print(f"Total other complete isolates: {len(other_isolates)}")

    # 6. Group USA isolates by month
    usa_by_month = {}
    for iso, segs in usa_isolates.items():
        d = context_dates[iso]
        m = d[:7] if d else "unknown"
        usa_by_month.setdefault(m, []).append(iso)

    print("\nUSA isolates count by month before downsampling:")
    for m in sorted(usa_by_month.keys()):
        print(f"  {m}: {len(usa_by_month[m])}")

    # 7. Apply downsampling rules
    kept_usa_isolates = set()
    print("\nApplying downsampling:")
    for m in sorted(usa_by_month.keys()):
        isolates_list = sorted(usa_by_month[m]) # sort to make random sample deterministic for the seed
        
        if m == "2022-03":
            n_keep = min(len(isolates_list), 150)
            reason = "cap to 150"
        elif m == "2022-04":
            n_keep = min(len(isolates_list), 150)
            reason = "cap to 150"
        elif m.startswith("2021"):
            n_keep = len(isolates_list)
            reason = "keep 100%"
        elif m == "2022-05":
            n_keep = max(1, int(len(isolates_list) * 0.25))
            reason = "keep 25%"
        else:
            n_keep = max(1, int(len(isolates_list) * 0.50))
            reason = f"keep 50%"

        kept = random.sample(isolates_list, n_keep)
        kept_usa_isolates.update(kept)
        print(f"  {m} ({reason}): original={len(isolates_list)}, kept={len(kept)}")

    # 8. Write the final set of isolates to FASTA
    print(f"\nWriting updated sequences to {fasta_path}...")
    temp_fasta_path = fasta_path + ".tmp"
    
    kept_count = 0
    with open(temp_fasta_path, "w", encoding="utf-8") as f_out:
        # Write other isolates
        for iso in sorted(other_isolates.keys()):
            for seg in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]:
                epi_isl, seq, hdr = other_isolates[iso][seg]
                f_out.write(f">{hdr}\n")
                f_out.write(seq.strip() + "\n")
            kept_count += 1
            
        # Write kept USA isolates
        for iso in sorted(kept_usa_isolates):
            for seg in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]:
                epi_isl, seq, hdr = usa_isolates[iso][seg]
                f_out.write(f">{hdr}\n")
                f_out.write(seq.strip() + "\n")
            kept_count += 1

    os.replace(temp_fasta_path, fasta_path)
    print(f"Successfully wrote {kept_count} isolates (8 segments each) to {fasta_path}")

if __name__ == "__main__":
    main()
