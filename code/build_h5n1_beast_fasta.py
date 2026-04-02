import argparse
import csv
import re
from collections import defaultdict

from date_normalization import parse_collection_date, pick_ecuador_date

REQUIRED_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
REQUIRED_SET = set(REQUIRED_ORDER)


def clean_id(text):
    text = "" if text is None else str(text).strip()
    if not text:
        return ""
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]", "", text)
    return text


def clean_place(text):
    text = "" if text is None else str(text).strip()
    if not text:
        return "UNKNOWN"
    text = text.replace(" ", "")
    text = text.replace("á", "a").replace("é", "e").replace("í", "i").replace("ó", "o").replace("ú", "u")
    text = text.replace("Á", "A").replace("É", "E").replace("Í", "I").replace("Ó", "O").replace("Ú", "U")
    text = text.replace("ñ", "n").replace("Ñ", "N")
    text = re.sub(r"[^A-Za-z0-9_-]", "", text)
    return text or "UNKNOWN"


def normalize_date(text):
    return parse_collection_date(text)


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


def read_flu_label_map(path, ecuador_date_source):
    labels = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample = row.get("Código USFQ") or row.get("Codigo USFQ") or row.get("codigo_usfq")
            sample = clean_id(sample)
            if not sample:
                continue

            place = clean_place(row.get("Provincia") or row.get("provincia") or row.get("province"))
            date_value = normalize_date(
                pick_ecuador_date(row, ecuador_date_source)
                or row.get("year")
                or row.get("año")
                or row.get("anio")
            )
            if not date_value:
                date_value = "UNKNOWN"

            labels[sample] = f"{sample}/{place}/{date_value}"
    return labels


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
            selection_role = clean_id(row.get("selection_role", ""))
            role_suffix = f"__{selection_role}" if selection_role else ""

            if accession:
                complete.add(accession)
                if role_suffix:
                    complete.add(f"{accession}{role_suffix}")
            if isolate and isolate.upper() != "UNKNOWN" and accession:
                complete.add(f"{isolate}_{accession}")
                if role_suffix:
                    complete.add(f"{isolate}_{accession}{role_suffix}")
    return complete


def read_context_label_map(path):
    labels = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            country = clean_place(row.get("country") or row.get("source_country") or "")
            date_value = normalize_date(row.get("collection_date") or row.get("year") or row.get("strain"))
            if not date_value:
                date_value = "UNKNOWN"

            accession = clean_id(row.get("accession", ""))
            isolate = clean_id(row.get("isolate", ""))
            selection_role = clean_id(row.get("selection_role", ""))
            role_suffix = f"__{selection_role}" if selection_role else ""

            if accession:
                key = accession
                labels[key] = f"{key}/{country}/{date_value}"
                if role_suffix:
                    key = f"{accession}{role_suffix}"
                    labels[key] = f"{key}/{country}/{date_value}"

            if isolate and isolate.upper() != "UNKNOWN" and accession:
                key = f"{isolate}_{accession}"
                labels[key] = f"{key}/{country}/{date_value}"
                if role_suffix:
                    key = f"{isolate}_{accession}{role_suffix}"
                    labels[key] = f"{key}/{country}/{date_value}"
    return labels


def main():
    parser = argparse.ArgumentParser(
        description="Build BEAST concatenated FASTA with complete 8-segment H5N1 samples"
    )
    parser.add_argument("--final-fasta", required=True)
    parser.add_argument("--flu-filtrado-csv", required=True)
    parser.add_argument("--context-metadata-tsv", required=True)
    parser.add_argument("--ecuador-date-source", default="reception")
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    flu_complete = read_flu_complete_samples(args.flu_filtrado_csv)
    flu_label_map = read_flu_label_map(args.flu_filtrado_csv, args.ecuador_date_source)
    context_complete = read_context_complete_samples(args.context_metadata_tsv)
    context_label_map = read_context_label_map(args.context_metadata_tsv)

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
            # Preserve original Flu headers (already enriched with place/date from input FASTA);
            # for context samples, use reconstructed label from metadata
            if from_flu and sample_id.startswith("Flu-"):
                output_header = sample_id
            else:
                output_header = context_label_map.get(sample_id, sample_id)
            selected.append((sample_id, output_header, concat))

        rows.append(
            {
                "sample_id": sample_id,
                "n_segments_present": len(present),
                "segments_present": ",".join(present),
                "is_complete_8": str(is_complete),
                "in_flu_filtrado_complete": str(from_flu),
                "in_context_metadata_segment8": str(from_context),
                "selected_for_beast": str(is_complete),
                "output_header": sample_id if (from_flu and sample_id.startswith("Flu-")) else context_label_map.get(sample_id, sample_id),
            }
        )

    with open(args.output_fasta, "w", encoding="utf-8") as out_fa:
        for sample_id, output_header, seq in selected:
            out_fa.write(f">{output_header}\n")
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
            "output_header",
        ]
        writer = csv.DictWriter(out_csv, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"BEAST FASTA: {args.output_fasta}")
    print(f"Summary: {args.output_summary}")
    print(f"Samples selected (8 complete): {len(selected)}")


if __name__ == "__main__":
    main()
