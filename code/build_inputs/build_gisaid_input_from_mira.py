import glob
import os
import re
from collections import defaultdict

import pandas as pd
import yaml


CONFIG_FILE = "config/config.yml"
OUTDIR = os.path.join("data", "standard_header_input_fasta")
STANDARD_GISAID_FASTA = os.path.join("data", "input", "H5N1_EC_gisaid_from_mira.fasta")

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
SEGMENT_SET = set(SEGMENTS)


def pick_column(df, candidates):
    lower_map = {c.lower().strip(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.lower().strip()
        if key in lower_map:
            return lower_map[key]
    return None


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


def clean_header_token(text):
    text = "" if text is None else str(text).strip()
    if not text or text.lower() in {"nan", "none", "na", "n/a"}:
        return "unknown"
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_\-]", "", text)
    return text or "unknown"


def normalize_host(species):
    species_key = clean_header_token(species).lower()
    host_aliases = {
        "aves_de_corral": "chicken",
        "pato": "duck",
        "patos": "duck",
        "condor": "condor",
        "fragata": "frigatebird",
        "fregata_magnificens": "frigatebird",
        "pardela_gris": "shearwater",
        "ardenna_grisea": "shearwater",
        "puffinus_sp": "shearwater",
        "piquero_peruano": "booby",
        "sula_nebouxii": "booby",
        "thalasseus_elegans": "tern",
        "gallinazo": "vulture",
    }
    return host_aliases.get(species_key, species_key)


def build_metadata_map(df):
    sample_col = pick_column(df, ["Código USFQ", "Codigo USFQ"])
    species_col = pick_column(df, ["Especie"])
    collection_col = pick_column(df, ["Fecha colección", "Fecha coleccion"])
    epi_col = pick_column(df, ["EPI_ISL"])

    missing_cols = [
        name
        for name, col in [
            ("Código USFQ", sample_col),
            ("Especie", species_col),
            ("Fecha colección", collection_col),
            ("EPI_ISL", epi_col),
        ]
        if col is None
    ]
    if missing_cols:
        raise SystemExit("Faltan columnas requeridas en flu_filtrado.csv: " + ", ".join(missing_cols))

    metadata = {}
    missing_epi = []
    bad_epi = []
    for _, row in df.iterrows():
        sample = normalize_sample_id(row.get(sample_col, ""))
        if not sample:
            continue

        epi_isl = str(row.get(epi_col, "")).strip()
        if not epi_isl:
            missing_epi.append(sample)
        elif not re.fullmatch(r"EPI_ISL_\d+", epi_isl):
            bad_epi.append(f"{sample}:{epi_isl}")

        date_value = str(row.get(collection_col, "")).strip()
        year_match = re.search(r"\d{4}", date_value)
        year = year_match.group(0) if year_match else "unknown"

        metadata[sample] = {
            "epi_isl": epi_isl,
            "virus_name": f"A/{normalize_host(row.get(species_col, ''))}/Ecuador/{sample}/{year}",
        }

    if missing_epi:
        raise SystemExit("Muestras sin EPI_ISL en flu_filtrado.csv: " + ", ".join(sorted(missing_epi)))
    if bad_epi:
        raise SystemExit("Valores EPI_ISL invalidos: " + ", ".join(sorted(bad_epi)))

    return metadata


def clean_seq(seq):
    return seq.upper().replace(" ", "").replace("\n", "").replace("-", "")


def parse_fasta(path):
    header = None
    seq_chunks = []

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    if header is not None:
        yield header, "".join(seq_chunks)


def parse_mira_header_sample_segment(header):
    sample = normalize_sample_id(header)
    if not sample:
        return None, None

    if "|" in header:
        last_field = header.rsplit("|", 1)[-1].strip().upper()
        if last_field in SEGMENT_SET:
            return sample, last_field

    token = header.strip().upper()
    for segment in SEGMENTS:
        if re.search(rf"(^|_){segment}(_|$)", token):
            return sample, segment

    return sample, None


def wrap_seq(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def standard_header(metadata, segment):
    return f"{metadata['virus_name']}|{segment}|{metadata['epi_isl']}"


def main():
    config = {}
    for filename in ["config/config.yml", "config/paths.yml"]:
        if os.path.exists(filename):
            with open(filename) as handle:
                config.update(yaml.safe_load(handle) or {})

    filtrado_csv = config.get("flu_filtrado", "metadata/flu_filtrado.csv")
    mira_base = config.get("mira_base_dir", "..")

    filtrado_df = pd.read_csv(filtrado_csv, dtype=str, keep_default_na=False)
    metadata_by_sample = build_metadata_map(filtrado_df)
    valid_samples = set(metadata_by_sample)
    print(f"Muestras cargadas de {filtrado_csv}: {len(valid_samples)}")

    os.makedirs(OUTDIR, exist_ok=True)

    mira_fastas = sorted(glob.glob(os.path.join(mira_base, "run*", "amended_consensus.fasta")))
    mira_fastas += sorted(glob.glob(os.path.join(mira_base, "run_agro", "amended_consensus.fasta")))
    mira_fastas = sorted(set(mira_fastas))

    if not mira_fastas:
        raise SystemExit(
            "No se encontraron amended_consensus.fasta de MIRA bajo "
            f"{mira_base}. Pon FASTA descargados de GISAID en data/input/ "
            "o ajusta mira_base_dir en config/config.yml."
        )

    print("FASTAs detectados:")
    for fasta in mira_fastas:
        print(" -", fasta)

    best_records = defaultdict(dict)
    sample_runs = defaultdict(set)
    raw_rows = []

    for fasta in mira_fastas:
        run_name = os.path.basename(os.path.dirname(fasta))
        for header, seq in parse_fasta(fasta):
            sample, segment = parse_mira_header_sample_segment(header)
            if sample is None or segment not in SEGMENT_SET:
                continue

            seq = clean_seq(seq)
            non_n = sum(1 for base in seq if base in "ACGT")
            raw_rows.append(
                {
                    "sample_norm": sample,
                    "run": run_name,
                    "segment": segment,
                    "header": header,
                    "seq_len": len(seq),
                    "non_n_bases": non_n,
                }
            )

            if non_n == 0:
                issue_rows.append({
                    "sample": sample,
                    "segment": segment,
                    "issue": "Segment sequence entirely Ns or empty in MIRA output; dropped."
                })
                continue

            sample_runs[sample].add(run_name)
            previous = best_records[sample].get(segment)
            if previous is None or non_n > previous["non_n_bases"]:
                best_records[sample][segment] = {
                    "run": run_name,
                    "header": header,
                    "sequence": seq,
                    "non_n_bases": non_n,
                }


    raw_df = pd.DataFrame(raw_rows)
    # raw_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_raw_segments.csv"), index=False)

    summary_rows = []
    audit_rows = []
    issue_rows = []
    all_found_samples = set(best_records)
    not_in_filtrado_samples = []
    discrepancy_rows = []
    sample_col = pick_column(filtrado_df, ["Código USFQ", "Codigo USFQ"])
    standard_records = []
    
    # Check for discrepancies with flu_filtrado.csv
    # We will build a map of what flu_filtrado.csv claims
    filtrado_claims = {s: {seg: False for seg in SEGMENTS} for s in valid_samples}
    for _, row in filtrado_df.iterrows():
        sample = normalize_sample_id(row.get(sample_col, ""))
        if sample in filtrado_claims:
            for seg in SEGMENTS:
                if str(row.get(seg, "")).strip().upper() == "SI":
                    filtrado_claims[sample][seg] = True

    # Check for false positives (in CSV but not in MIRA) and false negatives (not in CSV but in MIRA)
    for sample, claims in filtrado_claims.items():
        if sample not in best_records:
            if any(claims.values()):
                discrepancy_rows.append({"sample": sample, "issue": "Sample has 'SI' in CSV but is completely missing from MIRA FASTAs."})
            continue
        
        for seg in SEGMENTS:
            claimed = claims[seg]
            actual = seg in best_records[sample]
            
            if claimed and not actual:
                discrepancy_rows.append({"sample": sample, "issue": f"Segment {seg} marked as 'SI' in CSV but not found/assembled by MIRA."})
            elif not claimed and actual:
                discrepancy_rows.append({"sample": sample, "issue": f"Segment {seg} NOT marked as 'SI' in CSV but WAS assembled by MIRA."})
    for sample in sorted(best_records):
        if sample not in valid_samples:
            not_in_filtrado_samples.append(sample)
            continue

        metadata = metadata_by_sample[sample]
        segments_found = set(best_records[sample])
        runs = sorted(sample_runs[sample])
        row = {
            "Código USFQ": sample,
            "EPI_ISL": metadata["epi_isl"],
            "virus_name": metadata["virus_name"],
            "run": ",".join(runs),
        }

        for segment in SEGMENTS:
            if segment not in segments_found:
                row[segment] = ""
                continue

            row[segment] = "SI"
            record = best_records[sample][segment]
            header = standard_header(metadata, segment)
            seq = record["sequence"]

            standard_records.append((header, seq))
            audit_rows.append(
                {
                    "Código USFQ": sample,
                    "EPI_ISL": metadata["epi_isl"],
                    "virus_name": metadata["virus_name"],
                    "segment": segment,
                    "status": "assembled",
                    "source_run": record["run"],
                    "source_header": record["header"],
                    "sequence_length_written": len(seq),
                }
            )

        row["n_segments_assembled"] = sum(1 for segment in SEGMENTS if row.get(segment) == "SI")
        summary_rows.append(row)

        sample_filtrado = filtrado_df[filtrado_df[sample_col].map(normalize_sample_id) == sample].iloc[0]
        expected_segments = {
            segment
            for segment in SEGMENTS
            if str(sample_filtrado.get(segment, "")).strip().upper() == "SI"
        }
        missing_in_mira = expected_segments - segments_found
        if missing_in_mira:
            discrepancy_rows.append(
                {
                    "Código USFQ": sample,
                    "status": "PROCESSED_WITH_ISSUES",
                    "reason": "Segments marked SI in flu_filtrado but NOT found in MIRA: "
                    + ",".join(sorted(missing_in_mira)),
                    "found_in_mira": "PARTIAL",
                    "missing_segments": ",".join(sorted(missing_in_mira)),
                    "filled_with_Ns": "NO",
                }
            )

    for sample in not_in_filtrado_samples:
        issue_rows.append(
            {
                "Código USFQ": sample,
                "status": "SKIPPED",
                "reason": "Not in metadata/flu_filtrado.csv",
                "found_in_mira": "YES",
                "missing_segments": "",
                "filled_with_Ns": "",
            }
        )

    for sample in sorted(valid_samples):
        if sample not in all_found_samples:
            issue_rows.append(
                {
                    "Código USFQ": sample,
                    "status": "NOT_FOUND",
                    "reason": "Not found in MIRA amended_consensus.fasta files",
                    "found_in_mira": "NO",
                    "missing_segments": "ALL",
                    "filled_with_Ns": "NO",
                }
            )

    issue_rows.extend(discrepancy_rows)

    summary_df = pd.DataFrame(summary_rows)
    summary_columns = [
        "Código USFQ",
        "EPI_ISL",
        "virus_name",
        "run",
        *SEGMENTS,
        "n_segments_assembled",
    ]
    summary_df = summary_df.reindex(columns=summary_columns)
    # summary_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_summary.csv"), index=False)
    audit_columns = [
        "Código USFQ",
        "EPI_ISL",
        "virus_name",
        "segment",
        "status",
        "source_run",
        "source_header",
        "sequence_length_written",
    ]
    issue_columns = [
        "Código USFQ",
        "status",
        "reason",
        "found_in_mira",
        "missing_segments",
        "filled_with_Ns",
    ]
    # pd.DataFrame(audit_rows, columns=audit_columns).to_csv(
    #     os.path.join(OUTDIR, "ecuador_intermediate_audit.csv"),
    #     index=False,
    # )
    # pd.DataFrame(issue_rows, columns=issue_columns).to_csv(
    #     os.path.join(OUTDIR, "ecuador_intermediate_issues.csv"),
    #     index=False,
    # )

    with open(STANDARD_GISAID_FASTA, "w") as out:
        for header, seq in standard_records:
            out.write(f">{header}\n")
            out.write(wrap_seq(seq) + "\n")

    print("Listo.")
    print("Archivos generados:")
    print(f" - {STANDARD_GISAID_FASTA}")
    print()
    print("RESUMEN DE VALIDACION:")
    print(f"Muestras procesadas: {len(summary_df)}")
    print(f"Total secuencias (todos los segmentos combinados): {len(standard_records)}")

    # Print out warnings
    if issue_rows:
        print("\n--- ADVERTENCIAS DE CALIDAD (Secuencias vacías) ---")
        for i in issue_rows:
            if 'segment' in i:
                print(f"[{i['sample']} - {i['segment']}] {i['issue']}")

    if discrepancy_rows:
        print("\n--- INCONGRUENCIAS CON flu_filtrado.csv ---")
        for d in discrepancy_rows:
            print(f"[{d.get('sample', d.get('Código USFQ'))}] {d.get('issue', d.get('reason'))}")

    print("\nResumen por Muestra en 'ecuador_intermediate_summary.csv'")
    print(f"Muestras procesadas con discrepancias: {len(discrepancy_rows)}")
    print(f"Muestras no encontradas en MIRA: {len([r for r in issue_rows if r['status'] == 'NOT_FOUND'])}")
    print(f"Muestras descartadas por metadata: {len([r for r in issue_rows if r['status'] == 'SKIPPED'])}")


if __name__ == "__main__":
    main()
