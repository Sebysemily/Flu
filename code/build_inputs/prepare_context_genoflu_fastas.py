#!/usr/bin/env python3
"""Split complete context isolates into one multi-segment FASTA per isolate for GenoFLU."""

from __future__ import annotations

import argparse
import os
import re
import sys

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

from build_inputs.process_raw_to_segments import (  # noqa: E402
    filter_complete_context_isolates,
    parse_context_isolates,
    SEGMENTS,
    wrap_seq,
)


def safe_filename(isolate: str) -> str:
    name = re.sub(r"[^\w.\-]+", "_", isolate.strip())
    return name[:180] if name else "isolate"


def write_isolate_fasta(isolate: str, segs: dict, out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8") as handle:
        for seg in SEGMENTS:
            epi_isl, seq, _hdr = segs[seg]
            handle.write(f">{isolate} {seg}|{epi_isl}\n")
            handle.write(wrap_seq(seq.upper()) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write one GenoFLU input FASTA per complete context isolate."
    )
    parser.add_argument("--context-fasta", required=True)
    parser.add_argument("--flu-filtrado", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    import pandas as pd

    filtrado_df = pd.read_csv(args.flu_filtrado, dtype=str)
    local_epi_isls = set(filtrado_df["EPI_ISL"].dropna().astype(str).str.strip())

    isolates_data = parse_context_isolates(args.context_fasta)
    complete_context, *_ = filter_complete_context_isolates(isolates_data, local_epi_isls)

    os.makedirs(args.output_dir, exist_ok=True)
    for existing in os.listdir(args.output_dir):
        if existing.endswith(".fasta"):
            os.remove(os.path.join(args.output_dir, existing))

    written = 0
    for isolate, segs in sorted(complete_context.items()):
        out_path = os.path.join(args.output_dir, f"{safe_filename(isolate)}.fasta")
        write_isolate_fasta(isolate, segs, out_path)
        written += 1

    print(f"Wrote {written} isolate FASTAs to {args.output_dir}")


if __name__ == "__main__":
    main()
