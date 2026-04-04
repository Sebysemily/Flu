#!/usr/bin/env python3
import argparse
import os
import re
from typing import Dict, List, Set, Tuple


def read_fasta(path: str) -> Tuple[List[str], Dict[str, str]]:
    order: List[str] = []
    seqs: Dict[str, str] = {}
    current_id = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(chunks)
                    if current_id in seqs:
                        raise ValueError(f"Duplicate sequence id '{current_id}' in {path}")
                    order.append(current_id)
                    seqs[current_id] = seq
                current_id = line[1:].strip()
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(f"Invalid FASTA format in {path}: sequence before header")
                chunks.append(line)

    if current_id is not None:
        seq = "".join(chunks)
        if current_id in seqs:
            raise ValueError(f"Duplicate sequence id '{current_id}' in {path}")
        order.append(current_id)
        seqs[current_id] = seq

    if not seqs:
        raise ValueError(f"No sequences found in {path}")

    return order, seqs


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def ensure_directory(path: str) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)


def find_segment_path(paths: List[str], segment: str) -> str:
    pattern = re.compile(rf"(^|[_/\\\\]){re.escape(segment)}(?:\.|[_-])", re.IGNORECASE)
    matches = [p for p in paths if pattern.search(p)]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one alignment file for segment {segment}. "
            f"Found {len(matches)} matches: {matches}"
        )
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build concatenated codon-aware partitions for the main ML workflow"
    )
    parser.add_argument(
        "--segment-order",
        required=True,
        help="Comma-separated segment order, e.g. PB2,PB1,PA,HA,NP,NA,MP,NS",
    )
    parser.add_argument(
        "--codon-segments",
        required=True,
        help="Comma-separated segments that will use cp12/cp3 partitions",
    )
    parser.add_argument("--output-alignment", required=True)
    parser.add_argument("--output-partitions", required=True)
    parser.add_argument("alignments", nargs="+", help="Per-segment aligned FASTA files")
    args = parser.parse_args()

    segment_order = [s.strip() for s in args.segment_order.split(",") if s.strip()]
    codon_segments: Set[str] = {s.strip() for s in args.codon_segments.split(",") if s.strip()}

    if not segment_order:
        raise ValueError("Segment order is empty")
    if not codon_segments:
        raise ValueError("Codon segments set is empty")

    unknown = sorted(codon_segments - set(segment_order))
    if unknown:
        raise ValueError(f"Codon segments not in segment order: {unknown}")

    segment_paths = {segment: find_segment_path(args.alignments, segment) for segment in segment_order}

    per_segment: Dict[str, Dict[str, str]] = {}
    lengths_by_segment: Dict[str, int] = {}
    ids_by_segment: Dict[str, Set[str]] = {}

    for segment in segment_order:
        _, seg_seqs = read_fasta(segment_paths[segment])
        lengths = {len(seg_seqs[seq_id]) for seq_id in seg_seqs}
        if len(lengths) != 1:
            raise ValueError(f"Segment {segment} has non-uniform aligned lengths: {sorted(lengths)}")

        per_segment[segment] = seg_seqs
        lengths_by_segment[segment] = next(iter(lengths))
        ids_by_segment[segment] = set(seg_seqs.keys())

    common_ids = set.intersection(*(ids_by_segment[s] for s in segment_order))
    if not common_ids:
        raise ValueError("No shared taxa across all segment alignments")

    first_segment = segment_order[0]
    ids_order, _ = read_fasta(segment_paths[first_segment])
    ids_order = [seq_id for seq_id in ids_order if seq_id in common_ids]

    concatenated: Dict[str, str] = {seq_id: "" for seq_id in ids_order}
    partition_lines: List[str] = []
    offset = 1

    for segment in segment_order:
        seg_len = lengths_by_segment[segment]
        start = offset
        end = offset + seg_len - 1

        if segment in codon_segments:
            cp1 = start
            cp2 = start + 1
            cp3 = start + 2
            if cp3 > end:
                raise ValueError(f"Segment {segment} is too short for codon partitioning ({seg_len} nt)")

            partition_lines.append(
                f"GTR+G, {segment}_cp12 = {cp1}-{end}\\3,{cp2}-{end}\\3"
            )
            partition_lines.append(f"GTR+G, {segment}_cp3 = {cp3}-{end}\\3")
        else:
            partition_lines.append(f"GTR+G, {segment} = {start}-{end}")

        offset = end + 1

        seg_seqs = per_segment[segment]
        for seq_id in ids_order:
            concatenated[seq_id] += seg_seqs[seq_id]

    ensure_directory(args.output_alignment)
    with open(args.output_alignment, "w", encoding="utf-8") as handle:
        for seq_id in ids_order:
            handle.write(f">{seq_id}\n")
            handle.write(wrap_sequence(concatenated[seq_id]) + "\n")

    ensure_directory(args.output_partitions)
    with open(args.output_partitions, "w", encoding="utf-8") as handle:
        for line in partition_lines:
            handle.write(line + "\n")

    print(f"Wrote concatenated alignment: {args.output_alignment}")
    print(f"Wrote codon-aware partitions: {args.output_partitions}")


if __name__ == "__main__":
    main()
