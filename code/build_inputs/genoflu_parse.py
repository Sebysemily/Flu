#!/usr/bin/env python3
"""Shared helpers to parse GenoFLU-multi TSV output."""

from __future__ import annotations


import pandas as pd


SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def parse_genoflu_results(path: str) -> dict[str, dict[str, str]]:
    """Return strain -> {"_genotype": genotype} from a GenoFLU-multi results TSV.

    The special key ``_genotype`` holds the combined genotype assignment (e.g.
    ``B3.13``) taken directly from the ``Genotype`` column.  Per-segment
    lineage keys are no longer returned.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except pd.errors.EmptyDataError:
        return {}

    if df.empty:
        return {}

    if "Genotype" not in df.columns:
        return {}

    strain_to_lineages: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        strain = str(row["Strain"]).strip()
        genotype = str(row.get("Genotype", "")).strip() if pd.notna(row.get("Genotype")) else ""
        strain_to_lineages[strain] = {"_genotype": genotype}

    return strain_to_lineages


def lineages_by_epi_isl(strain_to_lineages: dict[str, dict[str, str]]) -> dict[str, dict[str, str]]:
    """Map EPI_ISL suffixes to segment lineages (Ecuador / per-segment strains)."""
    by_epi: dict[str, dict[str, str]] = {}
    for strain, lineages in strain_to_lineages.items():
        epi = strain.split("|")[-1].strip()
        if not epi.startswith("EPI_ISL_"):
            continue
        by_epi.setdefault(epi, {}).update(lineages)
    return by_epi
