#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import sys
import unicodedata
import pandas as pd

# Allow importing date_normalization from parent directory
_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import parse_collection_date, pick_ecuador_date  # noqa: E402

SEGMENTS = {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"}

MAATE_METADATA = {
    "EPI_ISL_17973443": {"province": "Guayas", "date": "2023-01-11"},
    "EPI_ISL_17973458": {"province": "Guayas", "date": "2023-01-11"},
    "EPI_ISL_18137671": {"province": "Guayas", "date": "2023-05-18"},
}

NORTH_AMERICA_PLACES = {
    "alabama", "alaska", "arizona", "arkansas", "california", "colorado", "connecticut",
    "delaware", "florida", "georgia", "hawaii", "idaho", "illinois", "indiana", "iowa",
    "kansas", "kentucky", "louisiana", "maine", "maryland", "massachusetts", "michigan",
    "minnesota", "mississippi", "missouri", "montana", "nebraska", "nevada", "new_hampshire",
    "new_jersey", "new_mexico", "new_york", "north_carolina", "north_dakota", "ohio", "oklahoma",
    "oregon", "pennsylvania", "rhode_island", "south_carolina", "south_dakota", "tennessee",
    "texas", "utah", "vermont", "virginia", "washington", "west_virginia", "wisconsin",
    "wyoming", "usa", "united_states", "honduras", "panama"
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
    "panama": "Panama",
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
    if compact_key in PLACE_ALIASES:
        return PLACE_ALIASES[compact_key]

    tokens = re.split(r"[^A-Za-z0-9]+", cleaned)
    tokens = [token for token in tokens if token]
    if not tokens:
        return "UNKNOWN"

    return "".join(token[:1].upper() + token[1:].lower() for token in tokens)


def read_fasta(path):
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
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


def sanitize_dna(seq):
    return re.sub(r"[^ACGTN-]", "N", str(seq).upper())


def extract_metadata_from_isolate(isolate):
    parts = isolate.split("/")
    if len(parts) < 3 or parts[0] != "A":
        return "UNKNOWN", "UNKNOWN", "regional_context"

    year_raw = parts[-1]
    year_match = re.search(r"\d{4}", year_raw)
    year = year_match.group(0) if year_match else "UNKNOWN"
    date_val = parse_collection_date(year) if year != "UNKNOWN" else "UNKNOWN"

    if len(parts) >= 5:
        place_raw = parts[2]
    elif len(parts) == 4:
        place_raw = parts[1]
    else:
        place_raw = parts[1]

    place = normalize_place(place_raw)

    place_clean = place.lower().replace("_", "").replace(" ", "")
    north_america_set = {s.lower().replace("_", "").replace(" ", "") for s in NORTH_AMERICA_PLACES}
    if place_clean in north_america_set:
        context_type = "american_anchor"
    else:
        context_type = "regional_context"

    return place, date_val, context_type


def parse_ecuador_header(header):
    # Parses GISAID-style headers for Ecuador sequences, e.g.:
    # A/booby/Ecuador/Flu-0008/2023|NS|EPI_ISL_20450066
    parts = [part.strip() for part in header.split("|")]
    if len(parts) != 3:
        return None
    virus_name, second, third = parts
    second_upper = second.upper()
    if second_upper in SEGMENTS and re.fullmatch(r"EPI_ISL_\d+", third):
        segment = second_upper
        epi_isl = third
    elif re.fullmatch(r"EPI_ISL_\d+", second) and third_upper in SEGMENTS:
        segment = third_upper
        epi_isl = second
    else:
        return None
    return segment, epi_isl


def main():
    parser = argparse.ArgumentParser(
        description="Unified process: ingest Ecuador and context FASTAs, merge, filter, and write per-segment FASTAs with EPI_ISL headers."
    )
    parser.add_argument("--ecuador-fastas", nargs="+", required=True)
    parser.add_argument("--context-fasta", required=True)
    parser.add_argument("--metadata-csv", required=True)
    parser.add_argument("--ecuador-date-source", default="collection")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--metadata-out", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.metadata_out), exist_ok=True)

    # 1. Parse Ecuador Metadata (flu_filtrado.csv)
    df_filt = pd.read_csv(args.metadata_csv, dtype=str)
    
    # Identify local EPI_ISLs to exclude context duplicates
    local_epi_isls = set()
    if "EPI_ISL" in df_filt.columns:
        local_epi_isls = set(df_filt["EPI_ISL"].dropna().str.strip())
        print(f"Loaded {len(local_epi_isls)} local EPI_ISLs from metadata.")

    # 2. Read Ecuador Sequences
    # Store by segment: { segment: { epi_isl: seq } }
    ecuador_by_segment = {seg: {} for seg in SEGMENTS}
    ecuador_seen = set()

    for fasta_path in args.ecuador_fastas:
        if not os.path.exists(fasta_path):
            continue
        for header, seq in read_fasta(fasta_path):
            parsed = parse_ecuador_header(header)
            if not parsed:
                continue
            seg, epi_isl = parsed
            if epi_isl not in local_epi_isls:
                raise ValueError(f"ERROR: Ecuador sequence {epi_isl} (header: {header}) not found in metadata file {args.metadata_csv}")
            ecuador_by_segment[seg][epi_isl] = sanitize_dna(seq).replace("-", "")
            ecuador_seen.add(epi_isl)

    print(f"Loaded Ecuador local sequences: {len(ecuador_seen)} unique isolates.")

    # 3. Read and Group GISAID Context sequences by isolate
    isolates_data = {} # isolate_name -> { segment -> (epi_isl, seq, original_header) }
    for header, seq in read_fasta(args.context_fasta):
        parts = [p.strip() for p in header.split("|")]
        if len(parts) < 3:
            continue
        isolate_name = parts[0]
        segment = None
        epi_isl = None
        for p in parts[1:]:
            p_upper = p.upper()
            if p_upper in SEGMENTS:
                segment = p_upper
            elif p.startswith("EPI_ISL_"):
                epi_isl = p
        
        if not segment or not epi_isl:
            for s in SEGMENTS:
                if s in header:
                    segment = s
            match = re.search(r"EPI_ISL_\d+", header)
            if match:
                epi_isl = match.group(0)

        if not segment or not epi_isl:
            continue

        isolates_data.setdefault(isolate_name, {})[segment] = (epi_isl, seq, header)

    print(f"Loaded context isolates: {len(isolates_data)}")

    # 4. Filter Context sequences to complete isolates
    complete_context = {}
    context_dates = {}
    context_places = {}
    context_types = {}
    
    for isolate, segs in isolates_data.items():
        if len(segs) == 8:
            # Exclude duplicate Ecuador isolates
            is_local = False
            for seg, (epi_isl, seq, orig_hdr) in segs.items():
                if epi_isl in local_epi_isls:
                    is_local = True
                    break
            if is_local:
                continue

            place, _, context_type = extract_metadata_from_isolate(isolate)
            place_clean = place.lower().replace("_", "").replace(" ", "")
            us_places = {s.lower().replace("_", "").replace(" ", "") for s in NORTH_AMERICA_PLACES if s not in {"honduras", "panama"}}
            if place_clean in us_places:
                place = "USA"

            epi_isl_rep = segs["HA"][0]
            is_maate = epi_isl_rep in MAATE_METADATA
            if is_maate:
                date_value = MAATE_METADATA[epi_isl_rep]["date"]
            else:
                extracted_date = None
                for seg_name, (epi_isl, seq, orig_hdr) in segs.items():
                    hdr_parts = [p.strip() for p in orig_hdr.split("|")]
                    if len(hdr_parts) >= 3:
                        raw_date = hdr_parts[2]
                        parsed_date = parse_collection_date(raw_date)
                        if parsed_date:
                            extracted_date = parsed_date
                            break
                date_value = extracted_date if extracted_date else "UNKNOWN"

            if date_value != "UNKNOWN" and date_value:
                complete_context[isolate] = segs
                context_dates[isolate] = date_value
                context_places[isolate] = place
                context_types[isolate] = context_type

    print(f"Complete context isolates kept: {len(complete_context)}")

    # 5. Write Context Metadata CSV
    # file_name, collection_date, country, expected_role
    context_meta_rows = []
    
    # We write a row for each sequence (each segment has a unique EPI_ISL)
    for isolate, segs in sorted(complete_context.items()):
        place = context_places[isolate]
        date_value = context_dates[isolate]
        context_type = context_types[isolate]
        
        for seg, (epi_isl, seq, orig_hdr) in segs.items():
            context_meta_rows.append({
                "file_name": epi_isl,
                "collection_date": date_value,
                "country": place,
                "expected_role": context_type
            })

    with open(args.metadata_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["file_name", "collection_date", "country", "expected_role"])
        writer.writeheader()
        writer.writerows(context_meta_rows)
    print(f"Wrote context metadata: {args.metadata_out}")

    # 6. Write Per-segment fastas containing only EPI_ISL as headers
    for seg in SEGMENTS:
        out_fasta = os.path.join(args.output_dir, f"H5N1_{seg}.fasta")
        with open(out_fasta, "w") as f:
            # A. Write Ecuador sequences
            for epi_isl, seq in sorted(ecuador_by_segment[seg].items()):
                f.write(f">{epi_isl}\n")
                f.write(wrap_seq(seq) + "\n")
            
            # B. Write Context sequences
            for isolate, segs in sorted(complete_context.items()):
                epi_isl, seq, orig_hdr = segs[seg]
                clean_seq = sanitize_dna(seq).replace("-", "")
                f.write(f">{epi_isl}\n")
                f.write(wrap_seq(clean_seq) + "\n")
        
        print(f"Wrote segment file: {out_fasta}")


if __name__ == "__main__":
    main()
