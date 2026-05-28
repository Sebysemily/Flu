#!/usr/bin/env python3
import argparse
import csv
import os
import sys

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import parse_collection_date, pick_ecuador_date  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build TreeTime dates TSV using metadata/H5N1_context.csv and config/flu_filtrado.csv"
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata/H5N1_context.csv")
    parser.add_argument("--ecuador-metadata", required=True, help="Path to config/flu_filtrado.csv")
    parser.add_argument("--ecuador-date-source", default="collection")
    parser.add_argument("--out", required=True, help="Output dates TSV")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    seen = set()
    rows = []

    # 1. Read GISAID context metadata
    with open(args.metadata, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            file_name = row.get("file_name")
            date_val = row.get("collection_date")
            if file_name and date_val:
                file_name = file_name.strip()
                date_val = date_val.strip()
                if file_name and date_val and file_name not in seen:
                    rows.append((file_name, date_val))
                    seen.add(file_name)

    # 2. Read Ecuador metadata
    with open(args.ecuador_metadata, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            epi_isl = row.get("EPI_ISL") or row.get("Epi_Isl")
            if not epi_isl:
                continue
            epi_isl = epi_isl.strip()
            if not epi_isl or epi_isl in seen:
                continue

            date_val = pick_ecuador_date(row, args.ecuador_date_source)
            if date_val:
                rows.append((epi_isl, date_val))
                seen.add(epi_isl)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("name\tdate\n")
        for name, date in sorted(rows):
            handle.write(f"{name}\t{date}\n")

    print(f"Wrote {len(rows)} date entries to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
