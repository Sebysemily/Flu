import argparse
import os
import re
import time
import unicodedata
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd


SEGMENT_MAP = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "MP",
    "8": "NS",
    "PB2": "PB2",
    "PB1": "PB1",
    "PA": "PA",
    "HA": "HA",
    "NP": "NP",
    "NA": "NA",
    "MP": "MP",
    "NS": "NS",
}

PLACE_ALIASES = {
    "argentina": "Argentina",
    "bolivia": "Bolivia",
    "brazil": "Brazil",
    "chile": "Chile",
    "colombia": "Colombia",
    "costarica": "CostaRica",
    "dominicanrepublic": "DominicanRepublic",
    "ecuador": "Ecuador",
    "elsalvador": "ElSalvador",
    "frenchguiana": "FrenchGuiana",
    "newzealand": "NewZealand",
    "peru": "Peru",
    "puertorico": "PuertoRico",
    "saudiarabia": "SaudiArabia",
    "southafrica": "SouthAfrica",
    "southkorea": "SouthKorea",
    "srilanka": "SriLanka",
    "trinidadandtobago": "TrinidadAndTobago",
    "unitedarabemirates": "UnitedArabEmirates",
    "unitedkingdom": "UnitedKingdom",
    "unitedstates": "UnitedStates",
    "uruguay": "Uruguay",
    "venezuela": "Venezuela",
}

INVALID_PLACE_KEYS = {
    "agrocalidad",
    "maate",
}


def clean_ascii(text):
    text = "" if text is None else str(text)
    text = text.strip()
    if text.lower() in {"", "nan", "none", "na", "n/a"}:
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


def parse_year(value):
    if value is None:
        return "UNKNOWN"
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "n/a"}:
        return "UNKNOWN"

    found = re.search(r"(19|20)\d{2}", text)
    return found.group(0) if found else "UNKNOWN"


def wrap_seq(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def sanitize_dna(seq):
    """Normalize sequence to uppercase A/C/G/T/N with gaps preserved."""
    return re.sub(r"[^ACGTN-]", "N", str(seq).upper())


def read_fasta(path):
    header = None
    chunks = []

    with open(path) as fh:
        for raw_line in fh:
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


def chunked(values, size):
    for i in range(0, len(values), size):
        yield values[i : i + size]


def parse_multi_fasta(text):
    header = None
    chunks = []

    for raw_line in text.splitlines():
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


def fetch_ncbi_fasta_batch(accessions, retries=1, pause_seconds=0.15, timeout=15):
    if not accessions:
        return {}

    params = {
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(params)

    for attempt in range(1, retries + 1):
        try:
            with urlopen(url, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
            if not text.strip() or text.lstrip().startswith("Error"):
                raise ValueError("NCBI devolvio respuesta vacia o error")

            parsed = list(parse_multi_fasta(text))
            if not parsed:
                raise ValueError("Respuesta FASTA invalida")

            # Agrupar posibles múltiples registros devueltos por accession base
            by_accession = {}
            for header, seq in parsed:
                seq = sanitize_dna(seq).replace("-", "")
                if not seq:
                    continue

                token = header.split()[0]
                token_base = token.split(".")[0]
                by_accession.setdefault(token_base, []).append((header, seq))

            time.sleep(pause_seconds)
            result = {}
            for acc in accessions:
                base = acc.split(".")[0]
                hits = by_accession.get(base, [])
                if not hits:
                    result[acc] = ([], "No se encontro secuencia en respuesta de NCBI")
                else:
                    result[acc] = (hits, None)
            return result

        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            if attempt == retries:
                return {acc: ([], str(exc)) for acc in accessions}
            time.sleep(max(1.0, pause_seconds * 2))

    return {acc: ([], "Fallo desconocido") for acc in accessions}


def normalize_segment(segment_value):
    raw = "" if segment_value is None else str(segment_value).strip().upper()
    return SEGMENT_MAP.get(raw, f"SEG_{clean_ascii(raw)}")


def parse_accession_list(value):
    text = "" if value is None else str(value).strip()
    if not text:
        return []
    parts = re.split(r"[;,\s]+", text)
    seen = set()
    out = []
    for part in parts:
        token = part.strip()
        if not token or token in seen:
            continue
        seen.add(token)
        out.append(token)
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Descarga secuencias contextuales por accession, las formatea con encabezado de influenza y combina con el FASTA de Ecuador"
    )
    parser.add_argument("--ecuador-fasta", required=True)
    parser.add_argument("--context-metadata-tsv", required=True)
    parser.add_argument("--context-fasta-out", required=True)
    parser.add_argument("--context-summary-out", required=True)
    parser.add_argument("--final-fasta-out", required=True)
    parser.add_argument("--batch-size", type=int, default=120)
    parser.add_argument("--http-timeout", type=int, default=15)
    parser.add_argument("--retries", type=int, default=1)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.context_fasta_out), exist_ok=True)
    os.makedirs(os.path.dirname(args.context_summary_out), exist_ok=True)
    os.makedirs(os.path.dirname(args.final_fasta_out), exist_ok=True)

    df = pd.read_csv(args.context_metadata_tsv, sep="\t", dtype=str)

    required = ["accession", "isolate", "segment", "collection_date", "source_country", "country"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Faltan columnas en metadata contextual: {', '.join(missing)}")

    # Download all contextual accessions in batches using all_accessions as primary source.
    accessions = []
    for _, row in df.iterrows():
        row_accessions = parse_accession_list(row.get("all_accessions", ""))
        primary_accession = str(row.get("accession", "")).strip()
        if primary_accession and primary_accession not in row_accessions:
            row_accessions.append(primary_accession)
        if not row_accessions and primary_accession:
            row_accessions = [primary_accession]
        accessions.extend(row_accessions)
    unique_accessions = sorted(set(accessions))

    downloaded = {}
    for batch in chunked(unique_accessions, max(1, args.batch_size)):
        downloaded.update(
            fetch_ncbi_fasta_batch(
                batch,
                retries=max(1, args.retries),
                timeout=max(5, args.http_timeout),
            )
        )

    records = []
    written = 0
    written_headers = set()

    with open(args.context_fasta_out, "w") as out_ctx:
        for _, row in df.iterrows():
            primary_accession = str(row.get("accession", "")).strip()
            if not primary_accession:
                continue

            row_accessions = parse_accession_list(row.get("all_accessions", ""))
            if primary_accession and primary_accession not in row_accessions:
                row_accessions.append(primary_accession)
            if not row_accessions:
                row_accessions = [primary_accession]

            isolate = clean_ascii(row.get("isolate", ""))
            # requested/annotated segment in TSV (fallback)
            requested_segment = normalize_segment(row.get("segment", ""))
            year = parse_year(row.get("collection_date", ""))

            source_country = normalize_place(row.get("source_country", ""))
            country = normalize_place(row.get("country", ""))
            place = source_country if source_country != "UNKNOWN" else country

            # Ensure sequence names are unique across context records.
            display_isolate = primary_accession if isolate == "UNKNOWN" else f"{isolate}_{primary_accession}"

            expected_count = None
            raw_segment_count = str(row.get("segment_count", "")).strip()
            if raw_segment_count.isdigit():
                expected_count = int(raw_segment_count)

            downloaded_segments = set()
            row_records_start = len(records)

            for accession in row_accessions:
                hits, error = downloaded.get(accession, ([], "No descargado"))
                status = "downloaded" if error is None and hits else "failed"

                if status == "downloaded":
                    for hit_header, hit_seq in hits:
                        # Deduce segment from returned header when possible
                        seg_match = re.search(r"segment\s*(\d)", hit_header, re.IGNORECASE)
                        gene_match = re.search(r"\b(PB2|PB1|PA|HA|NP|NA|MP|NS)\b", hit_header, re.IGNORECASE)
                        if seg_match:
                            seg = normalize_segment(seg_match.group(1))
                        elif gene_match:
                            seg = normalize_segment(gene_match.group(1).upper())
                        else:
                            seg = requested_segment

                        out_header = f"{display_isolate}/{seg}/{place}/{year}"
                        if out_header in written_headers:
                            continue

                        clean_hit_seq = sanitize_dna(hit_seq).replace("-", "")
                        if not clean_hit_seq:
                            continue

                        out_ctx.write(f">{out_header}\n")
                        out_ctx.write(wrap_seq(clean_hit_seq) + "\n")
                        written += 1
                        written_headers.add(out_header)
                        downloaded_segments.add(seg)

                        records.append(
                            {
                                "accession": accession,
                                "primary_accession": primary_accession,
                                "status": "downloaded",
                                "error": "",
                                "segment": seg,
                                "isolate": isolate,
                                "place": place,
                                "year": year,
                                "header": out_header,
                                "length": len(clean_hit_seq),
                                "ncbi_header": hit_header,
                                "expected_segment_count": expected_count,
                                "downloaded_segment_count": None,
                                "count_match": None,
                            }
                        )
                else:
                    out_header = f"{display_isolate}/{requested_segment}/{place}/{year}"
                    records.append(
                        {
                            "accession": accession,
                            "primary_accession": primary_accession,
                            "status": "failed",
                            "error": "" if error is None else error,
                            "segment": requested_segment,
                            "isolate": isolate,
                            "place": place,
                            "year": year,
                            "header": out_header,
                            "length": 0,
                            "ncbi_header": "",
                            "expected_segment_count": expected_count,
                            "downloaded_segment_count": None,
                            "count_match": None,
                        }
                    )

            downloaded_segment_count = len(downloaded_segments)
            count_match = None if expected_count is None else (downloaded_segment_count == expected_count)
            for i in range(row_records_start, len(records)):
                records[i]["downloaded_segment_count"] = downloaded_segment_count
                records[i]["count_match"] = count_match

    summary_df = pd.DataFrame(records)
    summary_df.to_csv(args.context_summary_out, index=False)

    mismatch_count = 0
    if not summary_df.empty and "count_match" in summary_df.columns:
        row_level = (
            summary_df[["primary_accession", "expected_segment_count", "downloaded_segment_count", "count_match"]]
            .drop_duplicates(subset=["primary_accession"])
            .dropna(subset=["expected_segment_count"])
        )
        mismatch_count = int((row_level["count_match"] == False).sum())

    # Final merge: Ecuador assembled FASTA + contextual FASTA.
    total_final = 0
    with open(args.final_fasta_out, "w") as out_final:
        for src in [args.ecuador_fasta, args.context_fasta_out]:
            if not os.path.exists(src):
                continue
            for header, seq in read_fasta(src):
                clean_seq = sanitize_dna(seq).replace("-", "")
                if not clean_seq:
                    continue
                out_final.write(f">{header}\n")
                out_final.write(wrap_seq(clean_seq) + "\n")
                total_final += 1

    print(f"Contexto descargado y formateado: {args.context_fasta_out}")
    print(f"Resumen de descarga: {args.context_summary_out}")
    print(f"FASTA final combinado: {args.final_fasta_out}")
    print(f"Registros de contexto escritos: {written}")
    print(f"Filas con desajuste segment_count vs descargado: {mismatch_count}")
    print(f"Registros totales en FASTA final: {total_final}")


if __name__ == "__main__":
    main()
