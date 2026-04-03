#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import sys
from typing import Dict, Optional

import pandas as pd
from Bio import SeqIO

from date_normalization import parse_collection_date, pick_ecuador_date


YEAR_RE = re.compile(r"(?:^|/)(19\d{2}|20\d{2})(?:$|/)")
ACCESSION_RE = re.compile(r"([A-Z]{1,2}\d{5,8}\.\d+)")


def load_ecuador_dates(flu_filtrado_csv: str, ecuador_date_source: str) -> Dict[str, str]:
    """Load Ecuador sample dates from flu_filtrado.csv using the configured Ecuador date source."""
    dates = {}
    try:
        df = pd.read_csv(flu_filtrado_csv)
        for _, row in df.iterrows():
            code = str(row.get("Código USFQ", "")).strip()
            date_value = pick_ecuador_date(row, ecuador_date_source)
            if code and code != "nan" and date_value:
                dates[code] = date_value
    except Exception as e:
        print(f"Warning: could not load Ecuador dates: {e}", file=sys.stderr)
    return dates


def load_context_dates(context_metadata_tsv: str) -> Dict[str, str]:
    """Load context sequence dates from final_metadata (accession → collection_date)."""
    dates = {}
    try:
        df = pd.read_csv(context_metadata_tsv, sep="\t")
        for _, row in df.iterrows():
            acc = str(row.get("accession", "")).strip()
            date_str = str(row.get("collection_date", "")).strip()
            if acc and acc != "nan" and date_str and date_str != "nan":
                parsed_date = parse_collection_date(date_str)
                if parsed_date:
                    dates[acc] = parsed_date
    except Exception as e:
        print(f"Warning: could not load context dates: {e}", file=sys.stderr)
    return dates

def extract_base_id(label: str) -> str:
    """Extract base Flu ID from header (e.g., 'Flu-0316' from 'Flu-0316/Cotopaxi/2023')."""
    if not label.startswith("Flu-"):
        return ""
    parts = label.split("/")
    return parts[0]


def extract_accession(label: str) -> Optional[str]:
    """Extract GenBank accession from header (e.g., 'PP779064.1' from '0155-N_PP779064.1__regional_context/Brazil/2024')."""
    match = ACCESSION_RE.search(label)
    if match:
        return match.group(1)
    return None


def build_dates_tsv_from_alignment(
    aln_path: str,
    out_path: str,
    ecuador_dates: Dict[str, str],
    context_dates: Dict[str, str],
) -> None:
    """Build dates TSV with hierarchical lookup: Ecuador exact dates > context metadata dates > year fallback."""
    rows = []
    seen = set()
    
    for rec in SeqIO.parse(aln_path, "fasta"):
        name = rec.id
        if name in seen:
            continue
        seen.add(name)
        
        date_value = None
        
        # Try Ecuador lookup by Flu-XXXX base ID
        if name.startswith("Flu-"):
            base_id = extract_base_id(name)
            if base_id in ecuador_dates:
                date_value = ecuador_dates[base_id]
        
        # Try context lookup by accession
        if not date_value:
            accession = extract_accession(name)
            if accession and accession in context_dates:
                date_value = context_dates[accession]
        
        # Fallback: extract year from header
        if not date_value:
            m = YEAR_RE.search(name)
            if m:
                date_value = f"{m.group(1)}-07-01"
        
        # If still no date, use generous fallback (don't fail yet)
        if not date_value:
            print(f"Warning: no date found for {name}, using 2023-07-01 as placeholder", file=sys.stderr)
            date_value = "2023-07-01"
        
        rows.append((name, date_value))
    
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("name\tdate\n")
        for name, date in rows:
            handle.write(f"{name}\t{date}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run TreeTime clock with date lookup from multiple sources")
    parser.add_argument("--treetime-exe", default="treetime")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--aln", required=True)
    parser.add_argument("--flu-filtrado-csv", default=None)
    parser.add_argument("--ecuador-date-source", default="reception")
    parser.add_argument("--context-metadata", default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--done", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    log_dir = os.path.dirname(args.log)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    # Load date lookups
    ecuador_dates = {}
    context_dates = {}
    if args.flu_filtrado_csv:
        ecuador_dates = load_ecuador_dates(args.flu_filtrado_csv, args.ecuador_date_source)
        print(f"Loaded {len(ecuador_dates)} Ecuador dates from {args.flu_filtrado_csv}", file=sys.stderr)
    if args.context_metadata:
        context_dates = load_context_dates(args.context_metadata)
        print(f"Loaded {len(context_dates)} context dates from {args.context_metadata}", file=sys.stderr)

    dates_path = os.path.join(args.outdir, "dates_from_headers.tsv")
    build_dates_tsv_from_alignment(args.aln, dates_path, ecuador_dates, context_dates)

    cmd = [
        args.treetime_exe,
        "clock",
        "--tree", args.tree,
        "--aln", args.aln,
        "--dates", dates_path,
        "--outdir", args.outdir,
    ]
    if args.seed is not None:
        cmd.extend(["--seed", str(args.seed)])

    with open(args.log, "w", encoding="utf-8", newline="") as log_handle:
        proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        raise SystemExit(proc.returncode)

    done_dir = os.path.dirname(args.done)
    if done_dir:
        os.makedirs(done_dir, exist_ok=True)
    with open(args.done, "w", encoding="utf-8") as handle:
        handle.write("done\n")


if __name__ == "__main__":
    main()