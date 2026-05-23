#!/usr/bin/env python3
import argparse
import csv
import re
from collections import Counter
from pathlib import Path


TIP_RE = re.compile(r"(?<=[(,])([^():,;]+):")
CODE_RE = re.compile(r"^(Flu-\d{4})/")


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
    with metadata_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        return {
            row.get("Código USFQ", "").strip()
            for row in rows
            if "condor" in row.get("Especie", "").lower()
            or "cóndor" in row.get("Especie", "").lower()
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
    if match:
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

    return "other"


def write_colorstrip(rows: list[tuple[str, str]], out_path: Path, dataset_label: str) -> None:
    present_groups = [g for g in GROUPS if any(group == g for _, group in rows)]
    legend_shapes = "\t".join(["1"] * len(present_groups))
    legend_colors = "\t".join(GROUPS[g]["color"] for g in present_groups)
    legend_labels = "\t".join(GROUPS[g]["label"] for g in present_groups)

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
    lines.extend(f"{tip}\t{GROUPS[group]['color']}\t{GROUPS[group]['label']}" for tip, group in rows)
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


def write_group_table(rows: list[tuple[str, str]], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxon", "group", "label", "color"])
        for tip, group in rows:
            writer.writerow([tip, group, GROUPS[group]["label"], GROUPS[group]["color"]])


def write_support_percent_tree(tree_path: Path, out_path: Path) -> None:
    text = tree_path.read_text(encoding="utf-8")
    text = re.sub(
        r"\)(0(?:\.\d+)?|1(?:\.0+)?):",
        lambda match: ")" + f"{float(match.group(1)) * 100:.1f}" + ":",
        text,
    )
    out_path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build iTOL annotation files for the H5N1 full-concat TBE tree."
    )
    parser.add_argument(
        "--tree",
        default="results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE",
    )
    parser.add_argument(
        "--outdir",
        default="results/phylogeny/raxml/full_concat/itol",
    )
    parser.add_argument(
        "--metadata",
        default="config/flu_filtrado.csv",
        help="Metadata CSV used to identify Fundacion Condor samples.",
    )
    args = parser.parse_args()

    tree_path = Path(args.tree)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tips = sorted(parse_tips(tree_path))
    condor_codes = read_condor_codes(Path(args.metadata))
    prefix = tree_path.name.replace(".raxml.supportTBE", "")
    rows = [(tip, classify_tip(tip, condor_codes, prefix)) for tip in tips]
    write_colorstrip(rows, outdir / f"{prefix}.itol_colorstrip_groups.txt", f"{prefix} groups")
    write_tree_colors(rows, outdir / f"{prefix}.itol_tree_colors.txt")
    write_group_table(rows, outdir / f"{prefix}.itol_group_assignments.tsv")
    write_support_percent_tree(tree_path, outdir / f"{prefix}.raxml.supportTBE.percent.nwk")

    counts = Counter(group for _, group in rows)
    with (outdir / f"{prefix}.itol_group_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "label", "color", "n_tips"])
        for group in GROUPS:
            if counts[group]:
                writer.writerow([group, GROUPS[group]["label"], GROUPS[group]["color"], counts[group]])


if __name__ == "__main__":
    main()
