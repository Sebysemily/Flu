import argparse
import os
from collections import Counter, defaultdict


SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def read_fasta(path):
    header = None
    chunks = []

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        yield header, "".join(chunks)


def wrap_seq(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def parse_header(header):
    parts = header.split("/")
    if len(parts) != 4:
        return None

    sample, segment, place, year = [part.strip() for part in parts]
    if not sample or not segment or not place or not year:
        return None

    return sample, segment.upper(), place, year


def main():
    parser = argparse.ArgumentParser(
        description="Separa el FASTA final en un archivo por segmento y quita el segmento del encabezado"
    )
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--summary-csv", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.summary_csv), exist_ok=True)

    grouped = defaultdict(list)
    invalid_headers = []

    for header, seq in read_fasta(args.input_fasta):
        parsed = parse_header(header)
        if parsed is None:
            invalid_headers.append(header)
            continue

        sample, segment, place, year = parsed
        grouped[segment].append((f"{sample}/{place}/{year}", seq))

    counts = Counter()
    for segment in SEGMENTS:
        out_path = os.path.join(args.output_dir, f"H5N1_{segment}.fasta")
        with open(out_path, "w") as out_handle:
            for out_header, seq in grouped.get(segment, []):
                out_handle.write(f">{out_header}\n")
                out_handle.write(wrap_seq(seq) + "\n")
                counts[segment] += 1

    with open(args.summary_csv, "w") as summary_handle:
        summary_handle.write("segment,n_sequences\n")
        for segment in SEGMENTS:
            summary_handle.write(f"{segment},{counts.get(segment, 0)}\n")
        if invalid_headers:
            summary_handle.write(f"INVALID_HEADERS,{len(invalid_headers)}\n")

    print(f"FASTA por segmento generados en: {args.output_dir}")
    print(f"Resumen generado: {args.summary_csv}")
    for segment in SEGMENTS:
        print(f"{segment}: {counts.get(segment, 0)} secuencias")
    if invalid_headers:
        print(f"Encabezados invalidos omitidos: {len(invalid_headers)}")


if __name__ == "__main__":
    main()