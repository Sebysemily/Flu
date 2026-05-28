#!/usr/bin/env python3
"""Shared helpers to parse GenoFLU-multi TSV output."""

from __future__ import annotations

import re

import pandas as pd

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def parse_genoflu_results(path: str) -> dict[str, dict[str, str]]:
    """Return strain -> {segment: lineage} from a GenoFLU-multi results TSV."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except pd.errors.EmptyDataError:
        return {}

    if df.empty:
        return {}

    list_col = next((c for c in df.columns if c.startswith("Genotype List Used")), None)
    if not list_col:
        return {}

    strain_to_lineages: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        strain = str(row["Strain"]).strip()
        val = row[list_col]
        if pd.notna(val) and val:
            parts = [p.strip() for p in re.split(r"[;,]", str(val)) if p.strip()]
            for part in parts:
                if ":" in part:
                    seg, lineage = part.split(":", 1)
                elif "=" in part:
                    seg, lineage = part.split("=", 1)
                else:
                    continue
                seg = seg.strip().upper()
                lineage = lineage.strip()
                if seg in SEGMENTS:
                    strain_to_lineages.setdefault(strain, {})[seg] = lineage
        strain_to_lineages.setdefault(strain, {})

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
