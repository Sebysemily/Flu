#!/usr/bin/env python3
import argparse
import csv
import os
import sys

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import parse_collection_date  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build TreeTime dates TSV from metadata/H5N1_context.csv"
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata/H5N1_context.csv")
    parser.add_argument("--out", required=True, help="Output dates TSV")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    seen = set()
    rows = []

    with open(args.metadata, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            file_name = (row.get("file_name") or "").strip()
            date_val = (row.get("collection_date") or "").strip()
            if file_name and date_val and file_name not in seen:
                rows.append((file_name, date_val))
                seen.add(file_name)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("name\tdate\n")
        for name, date in sorted(rows):
            handle.write(f"{name}\t{date}\n")

    print(f"Wrote {len(rows)} date entries to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
