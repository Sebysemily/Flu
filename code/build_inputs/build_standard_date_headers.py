import argparse
import os
import re
import unicodedata
import yaml

import pandas as pd

from date_normalization import (
    parse_collection_date,
    pick_ecuador_date,
    validate_no_missing_dates,
)


CONFIG_FILE = os.path.join("config", "config.yml")
DEFAULT_INPUT_FASTAS = [os.path.join("data", "input", "H5N1_EC_gisaid_from_mira.fasta")]
DEFAULT_METADATA_CSV = os.path.join("config", "flu_filtrado.csv")
DEFAULT_ECUADOR_DATE_SOURCE = "reception"
try:
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE) as fh:
            _cfg = yaml.safe_load(fh) or {}
            DEFAULT_METADATA_CSV = _cfg.get("flu_filtrado", DEFAULT_METADATA_CSV)
            DEFAULT_ECUADOR_DATE_SOURCE = _cfg.get("ecuador_date_source", DEFAULT_ECUADOR_DATE_SOURCE)
except Exception:
    pass

DEFAULT_OUTPUT_FASTA = os.path.join("data", "standard_header_input_fasta", "H5N1_EC.fasta")
DEFAULT_SUMMARY_CSV = None

SEGMENTS = {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"}

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


def normalize_date(date_value):
    parsed = parse_collection_date(date_value)
    return parsed if parsed else "UNKNOWN"


def normalize_sample_id(text):
    # Match numbered IDs: Flu-0743 → Flu-0743
    match_num = re.search(r"Flu-0*(\d+)", str(text))
    if match_num:
        return f"Flu-{int(match_num.group(1)):04d}"
    # Match GISAID-key IDs: Flu-EPI_ISL_17973443 → Flu-EPI_ISL_17973443
    match_gisaid = re.search(r"(Flu-EPI_ISL_\d+)", str(text))
    if match_gisaid:
        return match_gisaid.group(1)
    return None


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


def parse_gisaid_header(header):
    parts = [part.strip() for part in header.split("|")]
    if len(parts) != 3:
        return None

    virus_name, second, third = parts
    sample = normalize_sample_id(virus_name)
    if not sample:
        return None

    second_upper = second.upper()
    third_upper = third.upper()
    if second_upper in SEGMENTS and re.fullmatch(r"EPI_ISL_\d+", third):
        segment = second_upper
        epi_isl = third
    elif re.fullmatch(r"EPI_ISL_\d+", second) and third_upper in SEGMENTS:
        segment = third_upper
        epi_isl = second
    else:
        return None

    return {
        "sample": sample,
        "segment": segment,
        "epi_isl": epi_isl,
        "virus_name": virus_name,
    }


def build_metadata_map(metadata_csv, ecuador_date_source):
    df = pd.read_csv(metadata_csv, dtype=str, keep_default_na=False)
    source_value = (ecuador_date_source or "reception").strip().lower()

    sample_col = pick_column(df, ["Codigo USFQ", "Código USFQ"])
    province_col = pick_column(df, ["Provincia"])
    collection_col = pick_column(df, ["Fecha coleccion", "Fecha colección"])
    received_col = pick_column(df, ["Fecha recepcion", "Fecha recepción"])

    if sample_col is None:
        raise ValueError("No se encontro la columna de muestra (Codigo USFQ/Código USFQ) en metadata")
    if province_col is None:
        raise ValueError("No se encontro la columna Provincia en metadata")
    if source_value == "collection" and collection_col is None:
        raise ValueError(
            "ecuador_date_source=collection pero no se encontro la columna Fecha colección/Fecha coleccion en metadata"
        )
    if source_value == "reception" and received_col is None:
        raise ValueError(
            "ecuador_date_source=reception pero no se encontro la columna Fecha recepción/Fecha recepcion en metadata"
        )
    if source_value not in {"collection", "reception"}:
        raise ValueError(
            f"ecuador_date_source invalido: {ecuador_date_source}. Usa 'collection' o 'reception'."
        )

    validation_rows = []
    for _, row in df.iterrows():
        sample = normalize_sample_id(row.get(sample_col, ""))
        if not sample:
            continue
        validation_rows.append(
            {
                "sample": sample,
                "date_value": pick_ecuador_date(row, ecuador_date_source),
            }
        )
    validate_no_missing_dates(
        validation_rows,
        label_key="sample",
        date_key="date_value",
        context="headers FASTA de Ecuador",
    )

    metadata = {}
    for _, row in df.iterrows():
        sample = normalize_sample_id(row.get(sample_col, ""))
        if not sample:
            continue

        province = normalize_place(row.get(province_col, ""))
        date_value = normalize_date(pick_ecuador_date(row, ecuador_date_source))

        metadata[sample] = {
            "province": province,
            "date": date_value,
            "year": date_value[:4] if date_value != "UNKNOWN" else "UNKNOWN",
        }

    return metadata


def main():
    parser = argparse.ArgumentParser(
        description="Construir FASTA de Ecuador con headers sample/segment/province/date desde un FASTA GISAID"
    )
    parser.add_argument("--input-fasta", nargs="+", default=DEFAULT_INPUT_FASTAS)
    parser.add_argument("--metadata-csv", default=DEFAULT_METADATA_CSV)
    parser.add_argument("--ecuador-date-source", default=DEFAULT_ECUADOR_DATE_SOURCE)
    parser.add_argument("--output-fasta", default=DEFAULT_OUTPUT_FASTA)
    parser.add_argument("--summary-csv", default=DEFAULT_SUMMARY_CSV)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_fasta), exist_ok=True)
    if args.summary_csv:
        os.makedirs(os.path.dirname(args.summary_csv), exist_ok=True)

    metadata = build_metadata_map(args.metadata_csv, args.ecuador_date_source)

    parsed_records = {}
    invalid_headers = []
    missing_metadata = set()
    conflicting_records = []
    for fasta_path in args.input_fasta:
        for header, seq in read_fasta(fasta_path):
            parsed = parse_gisaid_header(header)
            if parsed is None:
                invalid_headers.append(f"{fasta_path}: {header}")
                continue
            if parsed["sample"] not in metadata:
                missing_metadata.add(parsed["sample"])
                continue

            key = (parsed["sample"], parsed["segment"])
            previous = parsed_records.get(key)
            if previous is not None:
                previous_parsed, previous_seq, previous_header, previous_path = previous
                if previous_parsed["epi_isl"] != parsed["epi_isl"] or previous_seq.upper() != seq.upper():
                    conflicting_records.append(
                        f"{parsed['sample']} {parsed['segment']}: {previous_path} vs {fasta_path}"
                    )
                continue

            parsed_records[key] = (parsed, seq, header, fasta_path)

    if invalid_headers:
        examples = "; ".join(invalid_headers[:5])
        raise SystemExit(
            "Headers invalidos en FASTA GISAID. Esperado: "
            "virus_name|segment|EPI_ISL o virus_name|EPI_ISL|segment. Ejemplos: "
            + examples
        )
    if missing_metadata:
        raise SystemExit(
            "Muestras del FASTA GISAID ausentes de flu_filtrado.csv: "
            + ", ".join(sorted(missing_metadata))
        )
    if conflicting_records:
        raise SystemExit(
            "Secuencias GISAID duplicadas con conflicto por muestra/segmento: "
            + "; ".join(sorted(conflicting_records)[:10])
        )
    if not parsed_records:
        raise SystemExit("No se encontraron secuencias validas en: " + ", ".join(args.input_fasta))

    rows = []
    with open(args.output_fasta, "w") as out:
        for key in sorted(parsed_records):
            parsed, seq, source_header, source_file = parsed_records[key]
            md = metadata[parsed["sample"]]
            out_header = f"{parsed['sample']}/{parsed['segment']}/{md['province']}/{md['date']}"

            out.write(f">{out_header}\n")
            out.write(wrap_seq(seq) + "\n")

            rows.append(
                {
                    "sample": parsed["sample"],
                    "segment": parsed["segment"],
                    "epi_isl": parsed["epi_isl"],
                    "virus_name": parsed["virus_name"],
                    "province": md["province"],
                    "date": md["date"],
                    "year": md["year"],
                    "header": out_header,
                    "source_header": source_header,
                    "source_file": source_file,
                    "length": len(seq),
                }
            )

    if args.summary_csv:
        pd.DataFrame(rows).to_csv(args.summary_csv, index=False)

    print(f"FASTA generado: {args.output_fasta}")
    if args.summary_csv:
        print(f"Resumen generado: {args.summary_csv}")
    print(f"Regiones escritas: {len(rows)}")


if __name__ == "__main__":
    main()
