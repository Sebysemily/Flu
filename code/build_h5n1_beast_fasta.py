import argparse
import csv
import re
from collections import defaultdict

REQUIRED_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
REQUIRED_SET = set(REQUIRED_ORDER)


def clean_id(text):
    text = "" if text is None else str(text).strip()
    if not text:
        return ""
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]", "", text)
    return text


def is_true_like(value):
    raw = "" if value is None else str(value).strip().lower()
    return raw in {"1", "true", "t", "yes", "y", "si", "s", "ok"}


def norm_segment(raw):
    token = "" if raw is None else str(raw).strip().upper()
    if token == "M":
        return "MP"
    return token


def wrap_seq(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def read_fasta(path):
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def parse_final_header(header):
    parts = header.split("/")
    if len(parts) < 2:
        return None, None
    sample_id = parts[0].strip()
    segment = norm_segment(parts[1])
    return sample_id, segment


def read_flu_complete_samples(path):
    complete = set()
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if not is_true_like(row.get("confirmed_has_8_regions", "")):
                continue
            sample = row.get("Código USFQ") or row.get("Codigo USFQ") or row.get("codigo_usfq")
            sample = clean_id(sample)
            if sample:
                complete.add(sample)
    return complete


def read_context_complete_samples(path):
    complete = set()
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            count = str(row.get("segment_count", "")).strip()
            if count != "8":
                continue

            accession = clean_id(row.get("accession", ""))
            isolate = clean_id(row.get("isolate", ""))

            if accession:
                complete.add(accession)
            if isolate and isolate.upper() != "UNKNOWN" and accession:
                complete.add(f"{isolate}_{accession}")
    return complete


def main():
    parser = argparse.ArgumentParser(
        description="Build BEAST concatenated FASTA with complete 8-segment H5N1 samples"
    )
    parser.add_argument("--final-fasta", required=True)
    parser.add_argument("--flu-filtrado-csv", required=True)
    parser.add_argument("--context-metadata-tsv", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    flu_complete = read_flu_complete_samples(args.flu_filtrado_csv)
    context_complete = read_context_complete_samples(args.context_metadata_tsv)

    sample_segments = defaultdict(dict)

    for header, seq in read_fasta(args.final_fasta):
        sample_id, segment = parse_final_header(header)
        if not sample_id or segment not in REQUIRED_SET:
            continue

        prev = sample_segments[sample_id].get(segment)
        if prev is None or len(seq) > len(prev):
            sample_segments[sample_id][segment] = seq

    rows = []
    selected = []

    for sample_id in sorted(sample_segments.keys()):
        segs = sample_segments[sample_id]
        present = sorted(segs.keys())
        is_complete = REQUIRED_SET.issubset(segs.keys())
        from_flu = sample_id in flu_complete
        from_context = sample_id in context_complete

        if is_complete:
            concat = "".join(segs[s] for s in REQUIRED_ORDER)
            selected.append((sample_id, concat))

        rows.append(
            {
                "sample_id": sample_id,
                "n_segments_present": len(present),
                "segments_present": ",".join(present),
                "is_complete_8": str(is_complete),
                "in_flu_filtrado_complete": str(from_flu),
                "in_context_metadata_segment8": str(from_context),
                "selected_for_beast": str(is_complete),
            }
        )

    with open(args.output_fasta, "w", encoding="utf-8") as out_fa:
        for sample_id, seq in selected:
            out_fa.write(f">{sample_id}\n")
            out_fa.write(wrap_seq(seq) + "\n")

    with open(args.output_summary, "w", encoding="utf-8", newline="") as out_csv:
        fieldnames = [
            "sample_id",
            "n_segments_present",
            "segments_present",
            "is_complete_8",
            "in_flu_filtrado_complete",
            "in_context_metadata_segment8",
            "selected_for_beast",
        ]
        writer = csv.DictWriter(out_csv, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"BEAST FASTA: {args.output_fasta}")
    print(f"Summary: {args.output_summary}")
    print(f"Samples selected (8 complete): {len(selected)}")


if __name__ == "__main__":
    main()
