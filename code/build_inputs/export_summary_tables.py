#!/usr/bin/env python3
import pandas as pd
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out_eco", required=True)
    args = parser.parse_args()

    meta = pd.read_csv(args.metadata)
    
    # 1. Ecological Groups
    eco_map = {}
    for _, row in meta.iterrows():
        ht = str(row.get('host_type', ''))
        h = str(row.get('host', '')).lower()
        if ht and ht != 'nan' and h and h != 'nan':
            if ht not in eco_map:
                eco_map[ht] = set()
            eco_map[ht].add(h)
            
    with open(args.out_eco, "w") as f:
        f.write('"Ecological_Group","Original_Species_from_GISAID"\n')
        for ht in sorted(eco_map.keys()):
            species_str = ", ".join(sorted(list(eco_map[ht])))
            f.write(f'"{ht}","{species_str}"\n')

if __name__ == "__main__":
    main()
