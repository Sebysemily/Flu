#!/usr/bin/env python3
"""
build_main_panel_metadata.py
----------------------------
Reads a Newick tree, extracts all tip labels, and classifies each taxon
into a metadata CSV with the following columns:

  id       – the taxon name stripped of the /__{tag} suffix and the trailing
              /location/date path  (i.e. the stable accession / sample ID)
  type_c   – one of:
               american_anchor  : label contains __american_anchor
               regional_context : label contains __regional_context
               flu_epi_isl      : Flu- prefix AND contains EPI_ISL in the id
               flu_sierra       : Flu- prefix, sierra sequences (default Flu-)
  date     – ISO date (YYYY-MM-DD) extracted from the 3rd slash-delimited
              field of the full label (the part after /location/).
              Compatible with TreeTime date format.

Usage:
    python build_main_panel_metadata.py \\
        --tree  results/phylogeny/iq-tree/HA/HA_final.treefile \\
        --output metadata/main_panel.csv
"""

import argparse
import re
import os
import csv


# ── Tag patterns (order matters: more specific first) ──────────────────
TAG_RULES = [
    ("american_anchor",  re.compile(r"__american_anchor", re.IGNORECASE)),
    ("regional_context", re.compile(r"__regional_context", re.IGNORECASE)),
    # coastal Ecuador: EPI_ISL-tagged Flu- sequences OR the hardcoded Flu-0406 sample
    ("flu_epi_isl",      re.compile(r"^Flu-.*epi_isl|^Flu-0406", re.IGNORECASE)),
    ("flu_sierra",       re.compile(r"^Flu-", re.IGNORECASE)),
]


def classify(label: str) -> str:
    for type_c, pat in TAG_RULES:
        if pat.search(label):
            return type_c
    return "unknown"


def extract_id(label: str) -> str:
    """
    Return the stable sample ID: the portion of the label before the first
    '__' tag token.  For Flu- sequences that have no tag, strip the trailing
    /location/date path.

    Examples
    --------
    'Flu-0582/Cotopaxi/2022-11-25'                                 → 'Flu-0582'
    'ABlue-winged_Teal...EPI_ISL_18133416__american_anchor/USA/…'  → 'ABlue-winged_Teal...EPI_ISL_18133416'
    'Achicken...EPI_ISL_17353838__regional_context/Magdalena/…'    → 'Achicken...EPI_ISL_17353838'
    """
    # Split on '__' first (for tagged sequences)
    if "__" in label:
        return label.split("__")[0]
    # For Flu- sequences: id is the first slash-segment
    return label.split("/")[0]


def extract_date(label: str) -> str:
    """
    The date is always the last slash-delimited field of the full label,
    and follows the ISO YYYY-MM-DD format used by TreeTime.

    For incomplete dates (YYYY or YYYY-MM) the value is kept as-is so
    TreeTime can handle ambiguous precision.
    """
    parts = label.split("/")
    candidate = parts[-1].strip()
    # Validate it looks like a date (at minimum YYYY)
    if re.match(r"^\d{4}", candidate):
        return candidate
    return ""


def extract_taxa_from_newick(tree_path: str) -> list:
    """Extract all tip labels from a Newick string via regex."""
    with open(tree_path, "r") as fh:
        content = fh.read()
    taxa = re.findall(r"(?<=[(,])([^():,;]+):", content)
    return [t.strip() for t in taxa if t.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Build main-panel metadata CSV from a Newick tree."
    )
    parser.add_argument("--tree", required=True,
                        help="Path to input Newick tree file.")
    parser.add_argument("--output", required=True,
                        help="Path to output metadata CSV.")
    args = parser.parse_args()

    if not os.path.exists(args.tree):
        raise FileNotFoundError(f"Tree not found: {args.tree}")

    taxa = extract_taxa_from_newick(args.tree)
    if not taxa:
        raise ValueError(f"No taxa extracted from {args.tree}. Check Newick format.")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    rows = []
    for label in taxa:
        rows.append({
            "label":  label,          # full original tip label (used as join key in R)
            "id":     extract_id(label),
            "type_c": classify(label),
            "date":   extract_date(label),
        })

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["label", "id", "type_c", "date"])
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    from collections import Counter
    counts = Counter(r["type_c"] for r in rows)
    print(f"Metadata written to {args.output}  ({len(rows)} taxa)")
    for k, v in sorted(counts.items()):
        print(f"  {k:20s}: {v}")


if __name__ == "__main__":
    main()
