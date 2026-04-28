#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

def main() -> None:
    parser = argparse.ArgumentParser(description="Run TreeTime clock")
    parser.add_argument("--treetime-exe", default="treetime")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--aln", required=True)
    parser.add_argument("--dates", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--done", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    cmd = [
        args.treetime_exe, "clock",
        "--tree", args.tree,
        "--aln", args.aln,
        "--dates", args.dates,
        "--outdir", args.outdir,
    ]

    with open(args.log, "w", encoding="utf-8") as log_handle:
        proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        raise SystemExit(proc.returncode)

    outliers_path = os.path.join(args.outdir, "outliers.tsv")
    if not os.path.exists(outliers_path):
        with open(outliers_path, "w", encoding="utf-8") as handle:
            handle.write("\tgiven_date\tapparent_date\tresidual\n")

    with open(args.done, "w", encoding="utf-8") as handle:
        handle.write("done\n")


if __name__ == "__main__":
    main()
