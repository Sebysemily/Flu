#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Parse GenoFLU-multi results and add {segment}_lineage columns to metadata file."
    )
    parser.add_argument("--genoflu-results", required=True, help="Path to metadata/genoflu_results.tsv")
    parser.add_argument("--metadata-csv", required=True, help="Path to metadata/flu_filtrado.csv")
    args = parser.parse_args()

    if not os.path.exists(args.genoflu_results):
        print(f"ERROR: GenoFLU results file not found: {args.genoflu_results}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.metadata_csv):
        print(f"ERROR: Metadata CSV file not found: {args.metadata_csv}", file=sys.stderr)
        sys.exit(1)

    # 1. Read and parse genoflu results
    df_genoflu = pd.read_csv(args.genoflu_results, sep='\t', dtype=str)
    
    # Identify the column containing the genotype list
    list_col = next((c for c in df_genoflu.columns if c.startswith("Genotype List Used")), None)
    if not list_col:
        print("ERROR: Could not find 'Genotype List Used' column in GenoFLU results.", file=sys.stderr)
        sys.exit(1)

    strain_to_lineages = {}
    for _, row in df_genoflu.iterrows():
        strain = str(row["Strain"]).strip()
        val = row[list_col]
        lineages = {}
        if pd.notna(val) and val:
            parts = str(val).split(",")
            for part in parts:
                part = part.strip()
                if ":" in part:
                    seg, lineage = part.split(":", 1)
                    lineages[seg.strip().upper()] = lineage.strip()
        strain_to_lineages[strain] = lineages

    # 2. Read metadata CSV
    df_meta = pd.read_csv(args.metadata_csv, dtype=str)
    
    # Find key column: 'EPI_ISL' or the second column
    if "EPI_ISL" in df_meta.columns:
        key_col = "EPI_ISL"
    elif len(df_meta.columns) >= 2:
        key_col = df_meta.columns[1]
    else:
        key_col = df_meta.columns[0]
        print(f"WARNING: Metadata CSV has only one column. Using {key_col} as key.", file=sys.stderr)

    print(f"Using column '{key_col}' as key matching to GenoFLU results.")

    # 3. Add segment lineage columns
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    for seg in segments:
        col_name = f"{seg}_lineage"
        if col_name not in df_meta.columns:
            df_meta[col_name] = ""

    # Populate columns
    for idx, row in df_meta.iterrows():
        key_val = str(row[key_col]).strip()
        if key_val in strain_to_lineages:
            lineages = strain_to_lineages[key_val]
            for seg in segments:
                col_name = f"{seg}_lineage"
                df_meta.at[idx, col_name] = lineages.get(seg, "")

    # 4. Save updated metadata CSV
    df_meta.to_csv(args.metadata_csv, index=False)
    print(f"Successfully updated metadata file: {args.metadata_csv}")

if __name__ == "__main__":
    main()
