#!/usr/bin/env python3
import argparse
import os


def read_first_sequence_length(path: str) -> int:
    length = 0
    in_first = False

    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_first:
                    break
                in_first = True
                continue
            if in_first:
                length += len(line)

    if length == 0:
        raise ValueError(f"Could not read first sequence length from {path}")

    return length


def ensure_directory(path: str) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Write cp12/cp3 partition file for one aligned segment")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--segment", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    seg_len = read_first_sequence_length(args.alignment)
    if seg_len < 3:
        raise ValueError(f"Segment {args.segment} is too short for codon partitioning: {seg_len}")

    cp12 = f"DNA, {args.segment}_cp12 = 1-{seg_len}\\3,2-{seg_len}\\3"
    cp3 = f"DNA, {args.segment}_cp3 = 3-{seg_len}\\3"

    ensure_directory(args.output)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(cp12 + "\n")
        handle.write(cp3 + "\n")

    print(f"Wrote segment codon partition file: {args.output}")


if __name__ == "__main__":
    main()
