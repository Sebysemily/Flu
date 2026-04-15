#!/usr/bin/env python3
import argparse
import os
import sys

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import extract_header_date  # noqa: E402


def read_fasta_headers(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or not line.startswith(">"):
                continue
            yield line[1:].strip()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build TreeTime dates TSV using only the dates embedded in FASTA headers"
    )
    parser.add_argument("--aln", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    seen = set()
    rows = []
    invalid = []

    for name in read_fasta_headers(args.aln):
        if name in seen:
            continue
        seen.add(name)

        date_value = extract_header_date(name)
        if not date_value:
            invalid.append(name)
            continue
        rows.append((name, date_value))

    if invalid:
        preview = ", ".join(invalid[:8])
        raise ValueError(
            f"Se encontraron {len(invalid)} headers sin fecha valida en el alignment subset: {preview}"
        )

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("name\tdate\n")
        for name, date in rows:
            handle.write(f"{name}\t{date}\n")

    print(f"Wrote {len(rows)} date entries to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
