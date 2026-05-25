#!/usr/bin/env python3
import argparse
import os
import re
import sys
import unicodedata
import pandas as pd

# Allow importing date_normalization from parent directory
_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import parse_collection_date  # noqa: E402

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
    """Parse place, date and context_type from a standard GISAID isolate name.

    Expected formats (slash-delimited, as used by GISAID):
      - A/host/place/id/year        (5 parts, host present)
      - A/place/id/year             (4 parts, no explicit host)

    The place field is always the token immediately before the isolate id
    (i.e. parts[-3] when len >= 4, or parts[-2] when len == 3).
    Year is always the last field.
    """
    parts = isolate.split("/")
    if len(parts) < 3 or parts[0] != "A":
        return "UNKNOWN", "UNKNOWN", "regional_context"

    # Year is always the last token
    year_raw = parts[-1]
    year_match = re.search(r"\d{4}", year_raw)
    year = year_match.group(0) if year_match else "UNKNOWN"
    date_val = parse_collection_date(year) if year != "UNKNOWN" else "UNKNOWN"

    # Place is the third-to-last token in a 5-part name (A/host/place/id/year)
    # or the second token in a 4-part name (A/place/id/year).
    if len(parts) >= 5:
        place_raw = parts[2]   # A / host / PLACE / id / year
    elif len(parts) == 4:
        place_raw = parts[1]   # A / PLACE / id / year
    else:
        place_raw = parts[1]   # fallback

    place = normalize_place(place_raw)

    # Classify: North America → american_anchor, everywhere else → regional_context
    place_clean = place.lower().replace("_", "").replace(" ", "")
    north_america_set = {s.lower().replace("_", "").replace(" ", "") for s in NORTH_AMERICA_PLACES}
    if place_clean in north_america_set:
        context_type = "american_anchor"
    else:
        context_type = "regional_context"

    return place, date_val, context_type


def main():
    parser = argparse.ArgumentParser(description="Ingest GISAID context fasta, validate 8 segments, format headers, merge with Ecuador.")
    parser.add_argument("--ecuador-fasta", required=True)
    parser.add_argument("--context-fasta-in", required=True,
                        help="Single concatenated GISAID context FASTA (e.g. data/context/gisaid_epiflu_sequence.fasta)")
    parser.add_argument("--context-fasta-out", required=True)
    parser.add_argument("--context-summary-out", default=None)
    parser.add_argument("--final-fasta-out", required=True)
    parser.add_argument("--metadata-csv", default=None)
    parser.add_argument("--max-per-country-month", type=int, default=None,
                        help="Maximum number of sequences to keep per country per month per year")
    parser.add_argument("--exclude-accessions", nargs="*", default=[],
                        help="List of EPI_ISL accessions to exclude from context sequences")
    args = parser.parse_args()

    if not os.path.exists(args.context_fasta_in):
        print(f"Error: context fasta not found: {args.context_fasta_in}")
        sys.exit(1)

    print(f"Leyendo contexto desde {args.context_fasta_in}")

    # Load local EPI_ISLs from metadata CSV to exclude them
    local_epi_isls = set()
    if args.metadata_csv and os.path.exists(args.metadata_csv):
        try:
            df_filt = pd.read_csv(args.metadata_csv, dtype=str)
            if "EPI_ISL" in df_filt.columns:
                local_epi_isls = set(df_filt["EPI_ISL"].dropna().str.strip())
                print(f"Loaded {len(local_epi_isls)} local EPI_ISLs to exclude from context.")
        except Exception as e:
            print(f"Warning: Could not read {args.metadata_csv} for duplicate exclusion: {e}")

    # Load excluded accessions
    exclude_accessions = set(args.exclude_accessions or [])
    if exclude_accessions:
        print(f"Loaded {len(exclude_accessions)} accessions to manually exclude from context.")

    # Parse and group by isolate name
    isolates_data = {} # isolate_name -> { segment -> (epi_isl, seq, original_header) }

    for header, seq in read_fasta(args.context_fasta_in):
        parts = [p.strip() for p in header.split("|")]
        if len(parts) < 3:
            continue

        # GISAID format usually: isolate_name|segment|EPI_ISL
        isolate_name = parts[0]

        # Identify segment
        segment = None
        epi_isl = None
        for p in parts[1:]:
            p_upper = p.upper()
            if p_upper in SEGMENTS:
                segment = p_upper
            elif p.startswith("EPI_ISL_"):
                epi_isl = p

        if not segment or not epi_isl:
            # Try fallback matching
            for s in SEGMENTS:
                if s in header:
                    segment = s
            match = re.search(r"EPI_ISL_\d+", header)
            if match:
                epi_isl = match.group(0)

        if not segment or not epi_isl:
            continue

        isolates_data.setdefault(isolate_name, {})[segment] = (epi_isl, seq, header)

    print(f"Isolados leídos: {len(isolates_data)}")

    # Filter to only isolates with all 8 segments
    complete_isolates = {}
    isolate_dates = {}
    isolate_places = {}
    isolate_types = {}
    discarded_count = 0
    excluded_local_count = 0
    for isolate, segs in isolates_data.items():
        if len(segs) == 8:
            is_local = False
            is_excluded = False
            for seg, (epi_isl, seq, orig_hdr) in segs.items():
                if epi_isl in local_epi_isls:
                    is_local = True
                    break
                if epi_isl in exclude_accessions:
                    is_excluded = True
                    break
            if is_local:
                excluded_local_count += 1
            elif is_excluded:
                continue
            else:
                place, _, context_type = extract_metadata_from_isolate(isolate)
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
                    if not extracted_date:
                        discarded_count += 1
                        continue
                    date_value = extracted_date
                
                if date_value == "UNKNOWN" or not date_value:
                    discarded_count += 1
                else:
                    complete_isolates[isolate] = segs
                    isolate_dates[isolate] = date_value
                    isolate_places[isolate] = place
                    isolate_types[isolate] = context_type
        else:
            discarded_count += 1

    print(f"Isolados con los 8 segmentos completos: {len(complete_isolates)}")
    print(f"Isolados descartados por incompletos (falta algún segmento): {discarded_count}")
    if excluded_local_count > 0:
        print(f"Isolados excluidos por coincidir con muestras locales (flu_filtrado): {excluded_local_count}")

    # Subsample complete context isolates if requested (max per country per month per year)
    if args.max_per_country_month is not None:
        subsampled_isolates = {}
        counts = {}
        # Sort isolates deterministically by name
        for isolate in sorted(complete_isolates.keys()):
            segs = complete_isolates[isolate]
            place = isolate_places[isolate]
            date_value = isolate_dates[isolate]
            context_type = isolate_types[isolate]
            
            epi_isl_rep = segs["HA"][0]
            is_maate = epi_isl_rep in MAATE_METADATA
            
            # Keep all local Ecuador core sequences
            if is_maate:
                subsampled_isolates[isolate] = segs
                continue
                
            # Parse year and month
            year = "UNKNOWN"
            month = "UNKNOWN"
            if date_value != "UNKNOWN" and len(date_value) >= 7:
                year = date_value[:4]
                month = date_value[5:7]
            
            # Normalize US states to USA for country-level filtering
            place_clean = place.lower().replace("_", "").replace(" ", "")
            north_america_set = {s.lower().replace("_", "").replace(" ", "") for s in NORTH_AMERICA_PLACES}
            country = "USA" if place_clean in north_america_set else place
            
            key = (country, year, month)
            counts[key] = counts.get(key, 0) + 1
            
            if counts[key] <= args.max_per_country_month:
                subsampled_isolates[isolate] = segs
                
        print(f"Subsampled context isolates from {len(complete_isolates)} to {len(subsampled_isolates)} using limit of {args.max_per_country_month} per country/month.")
        complete_isolates = subsampled_isolates

    # Process and write formatted context FASTA
    records = []
    written_count = 0
    os.makedirs(os.path.dirname(args.context_fasta_out), exist_ok=True)
    
    with open(args.context_fasta_out, "w") as out_ctx:
        for isolate, segs in sorted(complete_isolates.items()):
            place = isolate_places[isolate]
            date_value = isolate_dates[isolate]
            context_type = isolate_types[isolate]
            
            # Use the EPI_ISL of HA segment as the representative for the isolate ID
            epi_isl_rep = segs["HA"][0]
            clean_isolate = clean_ascii(isolate)
            
            is_maate = epi_isl_rep in MAATE_METADATA
            if is_maate:
                place = MAATE_METADATA[epi_isl_rep]["province"]
                date_value = MAATE_METADATA[epi_isl_rep]["date"]
            
            for seg, (epi_isl, seq, orig_hdr) in segs.items():
                if is_maate:
                    out_header = f"Flu-{epi_isl_rep}/{seg}/{place}/{date_value}"
                else:
                    out_header = f"{clean_isolate}_{epi_isl_rep}__{context_type}/{seg}/{place}/{date_value}"
                
                clean_seq = sanitize_dna(seq).replace("-", "")
                
                out_ctx.write(f">{out_header}\n")
                out_ctx.write(wrap_seq(clean_seq) + "\n")
                written_count += 1
                
                records.append({
                    "accession": epi_isl,
                    "isolate": isolate,
                    "segment": seg,
                    "place": place,
                    "date": date_value,
                    "year": date_value[:4] if date_value != "UNKNOWN" else "UNKNOWN",
                    "selection_role": "core" if is_maate else context_type,
                    "header": out_header,
                    "length": len(clean_seq),
                })

    if args.context_summary_out:
        os.makedirs(os.path.dirname(args.context_summary_out), exist_ok=True)
        pd.DataFrame(records).to_csv(args.context_summary_out, index=False)
    print(f"Contexto formateado escrito en {args.context_fasta_out} ({written_count} secuencias).")
    if args.context_summary_out:
        print(f"Resumen guardado en {args.context_summary_out}.")

    # Merge with Ecuador standard FASTA
    total_final = 0
    os.makedirs(os.path.dirname(args.final_fasta_out), exist_ok=True)
    with open(args.final_fasta_out, "w") as out_final:
        # 1. Write Ecuador core
        if os.path.exists(args.ecuador_fasta):
            for header, seq in read_fasta(args.ecuador_fasta):
                clean_seq = sanitize_dna(seq).replace("-", "")
                out_final.write(f">{header}\n")
                out_final.write(wrap_seq(clean_seq) + "\n")
                total_final += 1
        # 2. Write Context
        for header, seq in read_fasta(args.context_fasta_out):
            clean_seq = sanitize_dna(seq).replace("-", "")
            out_final.write(f">{header}\n")
            out_final.write(wrap_seq(clean_seq) + "\n")
            total_final += 1

    print(f"FASTA final combinado en {args.final_fasta_out} ({total_final} secuencias).")



if __name__ == "__main__":
    main()
