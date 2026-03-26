import argparse
import os
import re
import unicodedata
import yaml

import pandas as pd


DEFAULT_PER_SAMPLE_DIR = os.path.join("data", "assembled", "ecuador_intermediate_per_sample")
DEFAULT_AUDIT_CSV = os.path.join("data", "assembled", "ecuador_intermediate_audit.csv")
CONFIG_FILE = os.path.join("config", "config.yml")
DEFAULT_METADATA_CSV = os.path.join("config", "flu_filtrado.csv")
try:
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE) as fh:
            _cfg = yaml.safe_load(fh) or {}
            DEFAULT_METADATA_CSV = _cfg.get("flu_filtrado", DEFAULT_METADATA_CSV)
except Exception:
    # Fall back to builtin default if config can't be read
    pass

DEFAULT_OUTPUT_FASTA = os.path.join("data", "exports", "H5N1_EC.fasta")
DEFAULT_SUMMARY_CSV = os.path.join("data", "exports", "H5N1_EC_summary.csv")

PLACE_ALIASES = {
    "azuay": "Azuay",
    "bolivar": "Bolivar",
    "canar": "Canar",
    "carchi": "Carchi",
    "chimborazo": "Chimborazo",
    "cotopaxi": "Cotopaxi",
    "eloro": "ElOro",
    "esmeraldas": "Esmeraldas",
    "galapagos": "Galapagos",
    "guayas": "Guayas",
    "imbabura": "Imbabura",
    "loja": "Loja",
    "losrios": "LosRios",
    "manabi": "Manabi",
    "moronasantiago": "MoronaSantiago",
    "napo": "Napo",
    "orellana": "Orellana",
    "pastaza": "Pastaza",
    "pichincha": "Pichincha",
    "santaelena": "SantaElena",
    "santodomingo": "SantoDomingo",
    "santodomingodelostsachilas": "SantoDomingoDeLosTsachilas",
    "sucumbios": "Sucumbios",
    "tungurahua": "Tungurahua",
    "zamorachinchipe": "ZamoraChinchipe",
}

INVALID_PLACE_KEYS = {
    "agrocalidad",
    "maate",
}


def clean_ascii(text):
    text = "" if text is None else str(text)
    text = text.strip()
    if text.lower() in {"nan", "none", "na", "n/a"}:
        return "UNKNOWN"
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_\-]", "", text)
    return text or "UNKNOWN"


def normalize_place(text):
    cleaned = clean_ascii(text)
    if cleaned == "UNKNOWN":
        return cleaned

    compact_key = re.sub(r"[^A-Za-z0-9]", "", cleaned).lower()
    if compact_key in INVALID_PLACE_KEYS:
        return "UNKNOWN"
    if compact_key in PLACE_ALIASES:
        return PLACE_ALIASES[compact_key]

    tokens = re.split(r"[^A-Za-z0-9]+", cleaned)
    tokens = [token for token in tokens if token]
    if not tokens:
        return "UNKNOWN"

    return "".join(token[:1].upper() + token[1:].lower() for token in tokens)


def parse_year(date_value):
    if pd.isna(date_value):
        return "UNKNOWN"

    date_str = str(date_value).strip()
    if not date_str:
        return "UNKNOWN"

    # Try strict datetime parsing first; fall back to regex year extraction.
    parsed = pd.to_datetime(date_str, errors="coerce", dayfirst=False)
    if not pd.isna(parsed):
        return str(parsed.year)

    match = re.search(r"(19|20)\d{2}", date_str)
    return match.group(0) if match else "UNKNOWN"


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


def pick_column(df, candidates):
    lower_map = {c.lower().strip(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.lower().strip()
        if key in lower_map:
            return lower_map[key]
    return None


def parse_header_sample_segment(header):
    # Expected format from current pipeline: Flu-XXXX|A_SEG
    if "|" not in header:
        return None, None

    sample, right = header.split("|", 1)
    sample = sample.strip()

    segment = None
    m = re.search(r"_(PB2|PB1|PA|HA|NP|NA|MP|NS)$", right.strip().upper())
    if m:
        segment = m.group(1)

    return sample, segment


def build_metadata_map(metadata_csv):
    df = pd.read_csv(metadata_csv, dtype=str)

    sample_col = pick_column(df, ["Codigo USFQ", "Código USFQ"])
    province_col = pick_column(df, ["Provincia"])
    received_col = pick_column(df, ["Fecha recepcion", "Fecha recepción"])

    if sample_col is None:
        raise ValueError("No se encontro la columna de muestra (Codigo USFQ/Código USFQ) en metadata")
    if province_col is None:
        raise ValueError("No se encontro la columna Provincia en metadata")
    if received_col is None:
        raise ValueError("No se encontro la columna Fecha recepcion/Fecha recepción en metadata")

    metadata = {}
    for _, row in df.iterrows():
        sample = str(row.get(sample_col, "")).strip()
        if not sample:
            continue

        province = normalize_place(row.get(province_col, ""))
        year = parse_year(row.get(received_col, ""))

        metadata[sample] = {
            "province": province,
            "year": year,
        }

    return metadata


def build_assembled_set(audit_csv):
    # keep_default_na=False is critical so the segment label "NA" is not coerced to NaN.
    audit_df = pd.read_csv(audit_csv, dtype=str, keep_default_na=False)

    sample_col = pick_column(audit_df, ["Codigo USFQ", "Código USFQ"])
    segment_col = pick_column(audit_df, ["segment"])
    status_col = pick_column(audit_df, ["status"])

    if sample_col is None or segment_col is None or status_col is None:
        raise ValueError("El archivo de auditoria no tiene las columnas requeridas: Código USFQ, segment, status")

    assembled = set()
    for _, row in audit_df.iterrows():
        if str(row.get(status_col, "")).strip().lower() != "assembled":
            continue
        sample = str(row.get(sample_col, "")).strip()
        segment = str(row.get(segment_col, "")).strip().upper()
        if sample and segment:
            assembled.add((sample, segment))

    return assembled


def main():
    parser = argparse.ArgumentParser(description="Consolidar secuencias ensambladas de Ecuador con encabezado estilo influenza")
    parser.add_argument("--per-sample-dir", default=DEFAULT_PER_SAMPLE_DIR)
    parser.add_argument("--audit-csv", default=DEFAULT_AUDIT_CSV)
    parser.add_argument("--metadata-csv", default=DEFAULT_METADATA_CSV)
    parser.add_argument("--output-fasta", default=DEFAULT_OUTPUT_FASTA)
    parser.add_argument("--summary-csv", default=DEFAULT_SUMMARY_CSV)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(args.summary_csv), exist_ok=True)

    metadata = build_metadata_map(args.metadata_csv)
    assembled_set = build_assembled_set(args.audit_csv)

    fasta_files = sorted(
        f for f in os.listdir(args.per_sample_dir) if f.lower().endswith(".fasta")
    )

    rows = []
    total = 0

    with open(args.output_fasta, "w") as out:
        for fasta_name in fasta_files:
            fasta_path = os.path.join(args.per_sample_dir, fasta_name)

            for header, seq in read_fasta(fasta_path):
                sample, segment = parse_header_sample_segment(header)
                if not sample or not segment:
                    continue

                # Keep only truly assembled segments (skip segments filled with N).
                if (sample, segment) not in assembled_set:
                    continue

                md = metadata.get(sample, {"province": "UNKNOWN", "year": "UNKNOWN"})
                province = md["province"]
                year = md["year"]

                out_header = f"{sample}/{segment}/{province}/{year}"

                out.write(f">{out_header}\n")
                out.write(wrap_seq(seq) + "\n")
                total += 1

                rows.append(
                    {
                        "sample": sample,
                        "segment": segment,
                        "province": province,
                        "year": year,
                        "header": out_header,
                        "length": len(seq),
                    }
                )

    pd.DataFrame(rows).to_csv(args.summary_csv, index=False)

    print(f"FASTA generado: {args.output_fasta}")
    print(f"Resumen generado: {args.summary_csv}")
    print(f"Regiones ensambladas escritas: {total}")


if __name__ == "__main__":
    main()
