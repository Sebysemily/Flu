#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


TIP_RE = re.compile(r"(?<=[(,])([^():,;]+):")
CODE_RE = re.compile(r"^(Flu-\d{4})/")


TREE_SPECS = [
    (
        Path(f"results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE"),
        Path(f"results/phylogeny/raxml/{segment}/itol"),
    )
    for segment in ("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
] + [
    (
        Path("results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE"),
        Path("results/phylogeny/raxml/full_concat/itol"),
    ),
    (
        Path("results/phylogeny/raxml/np_mp_ns/H5N1_NP_MP_NS.raxml.supportTBE"),
        Path("results/phylogeny/raxml/np_mp_ns/itol"),
    ),
]


GROUPS = {
    "ecuador_clade1": {
        "label": "Ecuador Clade 1 / main complete-genome group",
        "color": "#1F78B4",
    },
    "ecuador_clade2_guayas": {
        "label": "Ecuador possible Clade 2 / Guayas-associated group",
        "color": "#FF1F1F",
    },
    "ecuador_clade3_condor": {
        "label": "Ecuador possible Clade 3 / Fundacion Condor",
        "color": "#FFD400",
    },
    "ecuador_2024_2025": {
        "label": "Ecuador 2024-2025 signal",
        "color": "#984EA3",
    },
    "ecuador_pichincha0694": {
        "label": "Ecuador Pichincha 0465/0694 signal",
        "color": "#00A6D6",
    },
    "ecuador_other_flu": {
        "label": "Ecuador Flu other",
        "color": "#7570B3",
    },
    "american_anchor": {
        "label": "American anchor",
        "color": "#F4A6C8",
    },
    "regional_context": {
        "label": "Regional context",
        "color": "#66A61E",
    },
    "usa_context": {
        "label": "USA neighbor/distal",
        "color": "#D8A21B",
    },
    "eurasian_anchor": {
        "label": "Eurasian anchor",
        "color": "#666666",
    },
    "other": {
        "label": "Other",
        "color": "#BDBDBD",
    },
}


GUAYAS_2023_CODES = {
    "Flu-0205",
    "Flu-0206",
    "Flu-0330",
    "Flu-0402",
    "Flu-0403",
    "Flu-0406",
    "Flu-0407",
    "Flu-0739",
}
PICHINCHA_0694_CODES = {"Flu-0465", "Flu-0694"}

HA_CLADE2_CODES = {
    "Flu-0023",
    "Flu-0402",
    "Flu-0406",
    "Flu-0407",
    "Flu-0450",
    "Flu-0465",
    "Flu-0472",
    "Flu-0578",
}

MP_CLADE2_CODES = {
    "Flu-0023",
    "Flu-0055",
    "Flu-0205",
    "Flu-0206",
    "Flu-0402",
    "Flu-0403",
    "Flu-0406",
    "Flu-0407",
    "Flu-0413",
    "Flu-0578",
}

MP_CLADE3_CODES = {
    "Flu-0679",
    "Flu-0680",
    "Flu-0688",
    "Flu-0711",
    "Flu-0713",
    "Flu-0714",
    "Flu-0715",
    "Flu-0716",
    "Flu-0717",
}


def parse_tips(tree_path: Path) -> list[str]:
    text = tree_path.read_text(encoding="utf-8").strip()
    return TIP_RE.findall(text)


def read_condor_codes(metadata_path: Path) -> set[str]:
    if not metadata_path.exists():
        return set()

    code_column = "C\u00f3digo USFQ"
    condor_markers = ("condor", "c\u00f3ndor")
    with metadata_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        return {
            row.get(code_column, "").strip()
            for row in rows
            if any(marker in row.get("Especie", "").lower() for marker in condor_markers)
        }


def classify_tip(tip: str, condor_codes: set[str], tree_prefix: str) -> str:
    if "__american_anchor" in tip:
        return "american_anchor"
    if "__usa_neighbor" in tip or "__usa_distal" in tip:
        return "usa_context"
    if "__eurasian_anchor" in tip:
        return "eurasian_anchor"
    if "__regional_context" in tip:
        return "regional_context"

    match = CODE_RE.match(tip)
    if not match:
        return "other"

    code = match.group(1)
    if tree_prefix == "H5N1_HA" and code in HA_CLADE2_CODES:
        return "ecuador_clade2_guayas"
    if tree_prefix == "H5N1_MP":
        if code in MP_CLADE3_CODES:
            return "ecuador_clade3_condor"
        if code in MP_CLADE2_CODES:
            return "ecuador_clade2_guayas"
        return "ecuador_clade1"
    if code in GUAYAS_2023_CODES or "/Guayas/" in tip:
        return "ecuador_clade2_guayas"
    if code in condor_codes:
        return "ecuador_clade3_condor"
    if code in PICHINCHA_0694_CODES:
        return "ecuador_pichincha0694"
    if re.search(r"/202[45]-", tip):
        return "ecuador_2024_2025"
    return "ecuador_clade1"


def write_colorstrip(rows: list[tuple[str, str]], out_path: Path, dataset_label: str) -> None:
    present_groups = [group for group in GROUPS if any(row_group == group for _, row_group in rows)]
    legend_shapes = "\t".join(["1"] * len(present_groups))
    legend_colors = "\t".join(GROUPS[group]["color"] for group in present_groups)
    legend_labels = "\t".join(GROUPS[group]["label"] for group in present_groups)

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{dataset_label}",
        "COLOR\t#000000",
        "STRIP_WIDTH\t35",
        "MARGIN\t5",
        "BORDER_WIDTH\t1",
        "BORDER_COLOR\t#FFFFFF",
        f"LEGEND_TITLE\t{dataset_label}",
        f"LEGEND_SHAPES\t{legend_shapes}",
        f"LEGEND_COLORS\t{legend_colors}",
        f"LEGEND_LABELS\t{legend_labels}",
        "DATA",
    ]
    lines.extend(
        f"{tip}\t{GROUPS[group]['color']}\t{GROUPS[group]['label']}" for tip, group in rows
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_tree_colors(rows: list[tuple[str, str]], out_path: Path) -> None:
    lines = [
        "TREE_COLORS",
        "SEPARATOR TAB",
        "DATA",
    ]
    for tip, group in rows:
        color = GROUPS[group]["color"]
        style = "bold" if group.startswith("ecuador_") else "normal"
        lines.append(f"{tip}\tlabel\t{color}\t{style}")
        if group.startswith("ecuador_") or group in {"american_anchor", "regional_context", "usa_context"}:
            lines.append(f"{tip}\tbranch\t{color}\tnormal")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def output_paths_for_tree(tree_path: Path, outdir: Path) -> tuple[Path, Path]:
    prefix = tree_path.name.replace(".raxml.supportTBE", "")
    return (
        outdir / f"{prefix}.itol_colorstrip_groups.txt",
        outdir / f"{prefix}.itol_tree_colors.txt",
    )


def build_annotations(tree_path: Path, outdir: Path, condor_codes: set[str]) -> None:
    if not tree_path.exists():
        raise FileNotFoundError(tree_path)

    outdir.mkdir(parents=True, exist_ok=True)
    prefix = tree_path.name.replace(".raxml.supportTBE", "")
    rows = [(tip, classify_tip(tip, condor_codes, prefix)) for tip in sorted(parse_tips(tree_path))]
    colorstrip_path, tree_colors_path = output_paths_for_tree(tree_path, outdir)
    write_colorstrip(rows, colorstrip_path, f"{prefix} groups")
    write_tree_colors(rows, tree_colors_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build iTOL colorstrip and tree-color TXT files for all ML H5N1 trees."
    )
    parser.add_argument(
        "--metadata",
        default="config/flu_filtrado.csv",
        help="Metadata CSV used to identify Fundacion Condor samples.",
    )
    args = parser.parse_args()

    condor_codes = read_condor_codes(Path(args.metadata))
    for tree_path, outdir in TREE_SPECS:
        build_annotations(tree_path, outdir, condor_codes)


if __name__ == "__main__":
    main()
