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

            by_accession = {}
            for header, seq in parsed:
                seq = re.sub(r"[^ACGTN]", "N", seq.upper())
                if not seq:
                    continue

                token = header.split()[0]
                token_base = token.split(".")[0]
                by_accession[token] = (header, seq)
                by_accession[token_base] = (header, seq)

            time.sleep(pause_seconds)
            result = {}
            for acc in accessions:
                base = acc.split(".")[0]
                hit = by_accession.get(acc) or by_accession.get(base)
                if hit is None:
                    result[acc] = (None, None, "No se encontro secuencia en respuesta de NCBI")
                else:
                    result[acc] = (hit[0], hit[1], None)
            return result

        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            if attempt == retries:
                return {acc: (None, None, str(exc)) for acc in accessions}
            time.sleep(max(1.0, pause_seconds * 2))

    return {acc: (None, None, "Fallo desconocido") for acc in accessions}


def normalize_segment(segment_value):
    raw = "" if segment_value is None else str(segment_value).strip().upper()
    return SEGMENT_MAP.get(raw, f"SEG_{clean_ascii(raw)}")


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

    # Download all contextual accessions in batches to reduce runtime and HTTP overhead.
    accessions = []
    for v in df["accession"].fillna("").astype(str).tolist():
        vv = v.strip()
        if vv:
            accessions.append(vv)
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

    with open(args.context_fasta_out, "w") as out_ctx:
        for _, row in df.iterrows():
            accession = str(row.get("accession", "")).strip()
            if not accession:
                continue

            isolate = clean_ascii(row.get("isolate", ""))
            segment = normalize_segment(row.get("segment", ""))
            year = parse_year(row.get("collection_date", ""))

            source_country = normalize_place(row.get("source_country", ""))
            country = normalize_place(row.get("country", ""))
            place = source_country if source_country != "UNKNOWN" else country

            remote_header, seq, error = downloaded.get(accession, (None, None, "No descargado"))
            status = "downloaded" if error is None else "failed"

            # Ensure sequence names are unique across context records.
            display_isolate = accession if isolate == "UNKNOWN" else f"{isolate}_{accession}"
            out_header = f"{display_isolate}/{segment}/{place}/{year}"

            if status == "downloaded":
                out_ctx.write(f">{out_header}\n")
                out_ctx.write(wrap_seq(seq) + "\n")
                written += 1

            records.append(
                {
                    "accession": accession,
                    "status": status,
                    "error": "" if error is None else error,
                    "segment": segment,
                    "isolate": isolate,
                    "place": place,
                    "year": year,
                    "header": out_header,
                    "length": 0 if seq is None else len(seq),
                    "ncbi_header": "" if remote_header is None else remote_header,
                }
            )

    pd.DataFrame(records).to_csv(args.context_summary_out, index=False)

    # Final merge: Ecuador assembled FASTA + contextual FASTA.
    total_final = 0
    with open(args.final_fasta_out, "w") as out_final:
        for src in [args.ecuador_fasta, args.context_fasta_out]:
            if not os.path.exists(src):
                continue
            for header, seq in read_fasta(src):
                out_final.write(f">{header}\n")
                out_final.write(wrap_seq(seq) + "\n")
                total_final += 1

    print(f"Contexto descargado y formateado: {args.context_fasta_out}")
    print(f"Resumen de descarga: {args.context_summary_out}")
    print(f"FASTA final combinado: {args.final_fasta_out}")
    print(f"Registros de contexto escritos: {written}")
    print(f"Registros totales en FASTA final: {total_final}")


if __name__ == "__main__":
    main()
