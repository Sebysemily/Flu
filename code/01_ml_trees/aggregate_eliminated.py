import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    dfs = []
    for f in args.inputs:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"Warning: could not read {f}: {e}")

    if not dfs:
        # Create empty with known columns
        df_out = pd.DataFrame(columns=["taxon", "country", "gap_n_fraction", "max_divergence", "filter_step", "discard_reason", "qc_action"])
    else:
        df_out = pd.concat(dfs, ignore_index=True)

    df_out.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
