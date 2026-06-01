#!/usr/bin/env python3
"""Build metadata/context_base.csv (GISAID only) and metadata/H5N1_context.csv."""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

from build_inputs.context_panel_metadata import (  # noqa: E402
    PANEL_COLUMNS,
    build_gisaid_context_rows,
)
from build_inputs.genoflu_parse import (  # noqa: E402
    SEGMENTS,
    lineages_by_epi_isl,
    parse_genoflu_results,
)
from build_inputs.process_raw_to_segments import (  # noqa: E402
    build_context_epi_maps,
    build_ecuador_metadata_rows,
    dedupe_metadata_rows,
    filter_complete_context_isolates,
    parse_context_isolates,
)

# flu_filtrado fields kept out of H5N1_context (clinical / lab workflow only)
FILTRADO_COLUMNS_EXCLUDE = {
    "Especie",
    "Tipo de muestra",
    "Estado salud animal",
    "observaciones",
    "Fecha colección",
    "Fecha recepción",
    "Semana Epidemiológica",
    "Provincia",
    "Ciudad",
    "Localidad",
    "Nombre_animal",
    "Género",
    "Edad",
    "MÉDICO CONTACTO",
    "Diagnóstico Kit Allplex",
    "Diagnóstico Kit Bioperfectuss/LunaSARS",
    "RT-PCR gen M",
    "Gel amplicon sequencing",
    "Secuenciamiento",
    "Código USFQ",
    "EPI_ISL",
    "Código procedencia",
    "run",
}



def get_clean_host(species: str) -> str:
    if not species or pd.isna(species):
        return "unknown"
    s = str(species).strip().lower().replace(" ", "_")
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
    return host_aliases.get(s, s)


def build_panel_rows(context_fasta: str, filtrado_df: pd.DataFrame, date_source: str) -> list[dict[str, str]]:
    local_epi_isls = set()
    if "EPI_ISL" in filtrado_df.columns:
        local_epi_isls = set(filtrado_df["EPI_ISL"].dropna().str.strip())

    ecuador_rows = build_ecuador_metadata_rows(filtrado_df, date_source)
    
    # Map local Ecuador sequences to their species/host
    epi_to_host = {}
    if "EPI_ISL" in filtrado_df.columns and "Especie" in filtrado_df.columns:
        for _, row in filtrado_df.iterrows():
            epi = str(row["EPI_ISL"]).strip() if pd.notna(row["EPI_ISL"]) else ""
            spec = str(row["Especie"]).strip() if pd.notna(row["Especie"]) else ""
            if epi and spec:
                epi_to_host[epi] = get_clean_host(spec)
                
    for row in ecuador_rows:
        row["host"] = epi_to_host.get(row["file_name"], "unknown")

    context_rows = build_gisaid_context_rows(context_fasta, local_epi_isls)

    panel_rows = dedupe_metadata_rows(ecuador_rows + context_rows)
    panel_rows.sort(key=lambda row: (row["expected_role"], row["file_name"]))
    return panel_rows


def filtrado_merge_columns(filtrado_df: pd.DataFrame) -> list[str]:
    return [
        col
        for col in filtrado_df.columns
        if col not in PANEL_COLUMNS and col not in FILTRADO_COLUMNS_EXCLUDE
    ]


def apply_context_lineages(
    row: dict[str, str],
    file_name: str,
    epi_to_segments: dict[str, set[str]],
    epi_to_isolate: dict[str, str],
    lineages_by_epi: dict[str, dict[str, str]],
    lineages_by_strain: dict[str, dict[str, str]],
) -> None:
    for segment in sorted(epi_to_segments.get(file_name, ())):
        row[segment] = "SI"

    lineages = lineages_by_epi.get(file_name)
    if not lineages:
        isolate = epi_to_isolate.get(file_name)
        if isolate:
            lineages = lineages_by_strain.get(isolate)

    if not lineages:
        return

    for seg in SEGMENTS:
        lineage = lineages.get(seg, "")
        if lineage:
            row[f"{seg}_lineage"] = lineage

    genotype = lineages.get("_genotype", "")
    if genotype:
        row["genotype"] = genotype



def write_panel_csv(rows: list[dict[str, str]], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    pd.DataFrame(rows, columns=PANEL_COLUMNS).to_csv(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build metadata/context_base.csv and metadata/H5N1_context.csv."
    )
    parser.add_argument("--context-fasta", required=True)
    parser.add_argument("--flu-filtrado", required=True)
    parser.add_argument("--genoflu-context-results", required=True)
    parser.add_argument("--ecuador-date-source", default="collection")
    parser.add_argument("--context-base-out", default="metadata/context_base.csv")
    parser.add_argument("--metadata-out", default="metadata/H5N1_context.csv")
    args = parser.parse_args()

    filtrado_df = pd.read_csv(args.flu_filtrado, dtype=str)
    local_epi_isls = set(filtrado_df["EPI_ISL"].dropna().str.strip()) if "EPI_ISL" in filtrado_df.columns else set()

    context_rows = build_gisaid_context_rows(args.context_fasta, local_epi_isls)
    write_panel_csv(context_rows, args.context_base_out)

    n_ecuador_gisaid = sum(1 for row in context_rows if row["expected_role"] == "flu_costa")
    print(
        f"Wrote {args.context_base_out}: {len(context_rows)} GISAID rows "
        f"(Ecuador coastal={n_ecuador_gisaid})"
    )

    panel_rows = build_panel_rows(args.context_fasta, filtrado_df, args.ecuador_date_source)

    filtrado_by_epi = {
        str(row["EPI_ISL"]).strip(): row
        for _, row in filtrado_df.iterrows()
        if pd.notna(row.get("EPI_ISL")) and str(row["EPI_ISL"]).strip()
    }

    isolates_data = parse_context_isolates(args.context_fasta)
    complete_context, *_ = filter_complete_context_isolates(isolates_data, local_epi_isls)
    epi_to_segments, epi_to_isolate = build_context_epi_maps(complete_context)

    lineages_by_strain = parse_genoflu_results(args.genoflu_context_results)
    lineages_by_epi = lineages_by_epi_isl(lineages_by_strain)

    merge_cols = filtrado_merge_columns(filtrado_df)
    output_columns = PANEL_COLUMNS + merge_cols + [f"{seg}_lineage" for seg in SEGMENTS]


    rows: list[dict[str, str]] = []
    for panel_row in panel_rows:
        file_name = panel_row["file_name"]
        row = {col: "" for col in output_columns}
        for col in PANEL_COLUMNS:
            row[col] = panel_row[col]

        if file_name in filtrado_by_epi:
            filt_row = filtrado_by_epi[file_name]
            for col in merge_cols:
                value = filt_row.get(col, "")
                row[col] = "" if pd.isna(value) else str(value).strip()
        else:
            apply_context_lineages(
                row, file_name, epi_to_segments, epi_to_isolate, lineages_by_epi, lineages_by_strain
            )

        rows.append(row)

    out_df = pd.DataFrame(rows, columns=output_columns)
    os.makedirs(os.path.dirname(args.metadata_out) or ".", exist_ok=True)
    out_df.to_csv(args.metadata_out, index=False)

    n_ecuador = sum(1 for row in rows if row["file_name"] in filtrado_by_epi)
    n_context_si = sum(
        1 for row in rows if row["file_name"] not in filtrado_by_epi and row["file_name"] in epi_to_segments
    )
    n_context_geno = sum(
        1 for row in rows if row["file_name"] not in filtrado_by_epi and row.get("genotype", "")
    )
    print(
        f"Wrote {args.metadata_out}: {len(rows)} rows "
        f"(Ecuador local={n_ecuador}, GISAID context={len(rows) - n_ecuador}, "
        f"context SI={n_context_si}, context with genotype={n_context_geno})"
    )


if __name__ == "__main__":
    main()
