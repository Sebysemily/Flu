#!/usr/bin/env python3
import argparse
import csv
import os
import re
from io import StringIO
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd
from Bio import Phylo


ECUADOR_CORE = [
    "Flu-0316", "Flu-0317", "Flu-0580", "Flu-0582", "Flu-0583", "Flu-0584", "Flu-0586",
    "Flu-0589", "Flu-0592", "Flu-0593", "Flu-0596", "Flu-0599", "Flu-0600", "Flu-0604", "Flu-0608",
    "Flu-0610", "Flu-0611", "Flu-0613", "Flu-0614", "Flu-0619", "Flu-0621", "Flu-0622", "Flu-0623",
    "Flu-0630", "Flu-0641", "Flu-0652", "Flu-0653", "Flu-0654",
]

REGIONAL_BLACKLIST_TOKENS = [
    "__eurasian_anchor", "__american_anchor", "__usa_",
]
ACCESSION_RE = re.compile(r"([A-Z]{1,2}\d{5,8}\.\d+)")
USA_DISTAL_QUOTA = 5


def read_tree(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read().strip()
    return Phylo.read(StringIO(text), "newick")


def read_complete_ids(path: str) -> Set[str]:
    df = pd.read_csv(path)
    # selected_for_beast uses True/False strings in this project
    mask = df["selected_for_beast"].astype(str).str.lower() == "true"
    # Use output_header column if available (tree tips match this format), fallback to sample_id
    id_col = "output_header" if "output_header" in df.columns else "sample_id"
    return set(df.loc[mask, id_col].astype(str))


def get_terminals(tree) -> Set[str]:
    return {t.name for t in tree.get_terminals() if t.name}


def flu_base_id(label: str) -> str:
    text = "" if label is None else str(label)
    if not text.startswith("Flu-"):
        return text
    return text.split("/", 1)[0]


def is_regional_context(label: str) -> bool:
    if "__regional_context" not in label:
        return False
    for token in REGIONAL_BLACKLIST_TOKENS:
        if token in label:
            return False
    return True


def is_usa_distal(label: str) -> bool:
    return "__usa_distal" in label


def nearest_candidates(
    tree,
    seeds: Iterable[str],
    candidates: Iterable[str],
    n_take: int,
    exclude: Set[str],
) -> List[Tuple[str, float, str]]:
    seed_list = [s for s in seeds if s in get_terminals(tree)]
    if not seed_list:
        return []

    scored: List[Tuple[str, float, str]] = []
    for cand in candidates:
        if cand in exclude:
            continue
        # distance to nearest seed
        dists = []
        for seed in seed_list:
            try:
                dists.append(tree.distance(cand, seed))
            except Exception:
                continue
        if not dists:
            continue
        nearest_seed = seed_list[dists.index(min(dists))]
        scored.append((cand, min(dists), nearest_seed))

    scored.sort(key=lambda x: x[1])
    return scored[:n_take]


def write_panel(path: str, rows: List[Tuple[str, str, str, float]]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxon", "role", "lineage", "distance_to_seed"])
        for taxon, role, lineage, dist in rows:
            writer.writerow([taxon, role, lineage, "" if dist is None else f"{dist:.6f}"])


def read_country_map(path: str) -> Dict[str, str]:
    df = pd.read_csv(path, sep=None, engine="python")
    country_col = "country" if "country" in df.columns else None
    accession_col = "accession" if "accession" in df.columns else None
    if country_col is None or accession_col is None:
        return {}
    out: Dict[str, str] = {}
    for _, row in df.iterrows():
        acc = str(row[accession_col]).strip()
        country = str(row[country_col]).strip()
        if acc and acc != "nan" and country and country != "nan":
            out[acc] = country
    return out


def taxon_country(label: str, country_map: Dict[str, str]) -> str:
    match = ACCESSION_RE.search(label)
    if match:
        return country_map.get(match.group(1), "UNKNOWN")
    return country_map.get(label, "UNKNOWN")


def cluster_members(
    tree,
    seeds: Iterable[str],
    candidates: Iterable[str],
    max_cluster_dist: float,
) -> List[Tuple[str, float, str]]:
    seed_list = [s for s in seeds if s in get_terminals(tree)]
    members: List[Tuple[str, float, str]] = []
    for cand in candidates:
        dists: List[float] = []
        for seed in seed_list:
            try:
                dists.append(tree.distance(cand, seed))
            except Exception:
                continue
        if not dists:
            continue
        min_dist = min(dists)
        if min_dist <= max_cluster_dist:
            nearest_seed = seed_list[dists.index(min_dist)]
            members.append((cand, min_dist, nearest_seed))
    members.sort(key=lambda x: x[1])
    return members


def main() -> None:
    parser = argparse.ArgumentParser(description="Build main BEAST panel from hybrid cluster/country selection")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--beast-summary", required=True)
    parser.add_argument("--context-metadata", required=True)
    parser.add_argument("--panel-main-out", required=True)
    parser.add_argument("--audit-out", required=True)
    parser.add_argument("--max-cluster-dist", type=float, default=0.08)
    parser.add_argument("--n-per-country", type=int, default=4)
    parser.add_argument("--n-total", type=int, default=60)
    args = parser.parse_args()

    tree = read_tree(args.tree)
    complete_ids = read_complete_ids(args.beast_summary)
    country_map = read_country_map(args.context_metadata)
    tree_tips = get_terminals(tree)

    flu_tip_map: Dict[str, str] = {}
    for tip in tree_tips:
        if tip.startswith("Flu-"):
            flu_tip_map[flu_base_id(tip)] = tip

    core_available = []
    core_missing = []
    for core in ECUADOR_CORE:
        tip_label = flu_tip_map.get(core, core)
        if core in complete_ids and tip_label in tree_tips:
            core_available.append(tip_label)
        else:
            core_missing.append(core)

    seed_main = list(core_available)

    regional_candidates = [
        t for t in tree_tips
        if t in complete_ids and is_regional_context(t)
    ]
    usa_distal_candidates = [
        t for t in tree_tips
        if t in complete_ids and is_usa_distal(t)
    ]

    members = cluster_members(
        tree=tree,
        seeds=seed_main,
        candidates=regional_candidates,
        max_cluster_dist=args.max_cluster_dist,
    )

    selected_context: List[Tuple[str, float, str, str]] = []
    selected_ids: Set[str] = set()

    # 1) Country diversity first.
    by_country: Dict[str, List[Tuple[str, float, str]]] = {}
    for taxon, dist, seed in members:
        country = taxon_country(taxon, country_map)
        by_country.setdefault(country, []).append((taxon, dist, seed))
    for country in by_country:
        by_country[country].sort(key=lambda x: x[1])

    for country, rows in sorted(by_country.items(), key=lambda x: x[0]):
        take = 0
        for taxon, dist, seed in rows:
            if taxon in selected_ids:
                continue
            if take >= args.n_per_country:
                break
            selected_context.append((taxon, dist, seed, "country_diversity"))
            selected_ids.add(taxon)
            take += 1
            if len(selected_context) >= args.n_total:
                break
        if len(selected_context) >= args.n_total:
            break

    # 2) Complete by nearest phylogenetic distance.
    if len(selected_context) < args.n_total:
        for taxon, dist, seed in members:
            if taxon in selected_ids:
                continue
            selected_context.append((taxon, dist, seed, "distance_only"))
            selected_ids.add(taxon)
            if len(selected_context) >= args.n_total:
                break

    selected_usa_distal = nearest_candidates(
        tree=tree,
        seeds=seed_main,
        candidates=usa_distal_candidates,
        n_take=USA_DISTAL_QUOTA,
        exclude=set(),
    )

    panel_main_rows: List[Tuple[str, str, str, float]] = []
    for t in core_available:
        panel_main_rows.append((t, "ecuador_core", "main_cluster", None))
    for taxon, dist, seed, reason in selected_context:
        panel_main_rows.append((taxon, "regional_context", reason, dist))
    for taxon, dist, seed in selected_usa_distal:
        panel_main_rows.append((taxon, "usa_distal", "distance_to_ecuador_core", dist))

    write_panel(args.panel_main_out, panel_main_rows)

    audit_dir = os.path.dirname(args.audit_out)
    if audit_dir:
        os.makedirs(audit_dir, exist_ok=True)
    with open(args.audit_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["core_available", len(core_available)])
        writer.writerow(["core_missing", len(core_missing)])
        writer.writerow(["cluster_members", len(members)])
        writer.writerow(["panel_main_total", len(panel_main_rows)])
        writer.writerow(["max_cluster_dist", args.max_cluster_dist])
        writer.writerow(["n_per_country", args.n_per_country])
        writer.writerow(["n_total", args.n_total])
        writer.writerow(["usa_distal_quota", USA_DISTAL_QUOTA])
        writer.writerow(["selected_country_diversity", sum(1 for x in selected_context if x[3] == "country_diversity")])
        writer.writerow(["selected_distance_only", sum(1 for x in selected_context if x[3] == "distance_only")])
        writer.writerow(["selected_usa_distal", len(selected_usa_distal)])
        # country coverage report
        country_counts: Dict[str, int] = {}
        for taxon, _, _, _ in selected_context:
            c = taxon_country(taxon, country_map)
            country_counts[c] = country_counts.get(c, 0) + 1
        for c, n in sorted(country_counts.items()):
            writer.writerow([f"country_count::{c}", n])
        if core_missing:
            writer.writerow(["core_missing_ids", ";".join(core_missing)])


if __name__ == "__main__":
    main()
