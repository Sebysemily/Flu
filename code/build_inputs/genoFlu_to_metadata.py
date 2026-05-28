#!/usr/bin/env python3
import argparse
import os
import sys

import pandas as pd

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

from build_inputs.genoflu_parse import SEGMENTS, lineages_by_epi_isl, parse_genoflu_results  # noqa: E402


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

    strain_to_lineages = parse_genoflu_results(args.genoflu_results)
    if not strain_to_lineages:
        print("INFO: GenoFLU results contain no rows or missing genotype column; nothing to add.")
        sys.exit(0)

    by_epi = lineages_by_epi_isl(strain_to_lineages)

    df_meta = pd.read_csv(args.metadata_csv, dtype=str)

    if "EPI_ISL" in df_meta.columns:
        key_col = "EPI_ISL"
    elif len(df_meta.columns) >= 2:
        key_col = df_meta.columns[1]
    else:
        key_col = df_meta.columns[0]
        print(f"WARNING: Metadata CSV has only one column. Using {key_col} as key.", file=sys.stderr)

    print(f"Using column '{key_col}' as key matching to GenoFLU results.")

    for seg in SEGMENTS:
        col_name = f"{seg}_lineage"
        if col_name not in df_meta.columns:
            df_meta[col_name] = ""

    for idx, row in df_meta.iterrows():
        key_val = str(row[key_col]).strip()
        if key_val in by_epi:
            lineages = by_epi[key_val]
            for seg in SEGMENTS:
                col_name = f"{seg}_lineage"
                df_meta.at[idx, col_name] = lineages.get(seg, "")

    df_meta.to_csv(args.metadata_csv, index=False)
    print(f"Successfully updated metadata file: {args.metadata_csv}")


if __name__ == "__main__":
    main()
