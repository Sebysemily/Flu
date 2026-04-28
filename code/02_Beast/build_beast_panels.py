#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from io import StringIO
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd
from Bio import Phylo

# Allow importing date_normalization from the parent code/ directory
_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import parse_collection_date  # noqa: E402


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

BEAST_PANEL_PROFILE = {
    "label": "relaxed",
    "n_per_country": 12,
    "n_total": 110,
    "min_mrca_support": 50.0,
    "max_per_country_month": 4,
    "usa_distal_quota": 0,
    "additional_american_anchor_quota": 5,
    "forced_american_anchor_accession": "OQ968009",
}


def read_tree(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read().strip()
    return Phylo.read(StringIO(text), "newick")


def canonical_tip(label: str) -> str:
    """Lowercase only the country segment (second slash-delimited field) of a tip label
    so that /USA/ and /Usa/ compare as equal."""
    if not label:
        return label
    parts = label.split("/")
    if len(parts) >= 2:
        parts[1] = parts[1].lower()
    return "/".join(parts)


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


def is_american_anchor(label: str) -> bool:
    return "__american_anchor" in label


def normalize_support(v) -> float:
    """Normalize TBE support to 0-100 scale regardless of whether source is 0-1 or 0-100."""
    if v is None:
        return 0.0
    f = float(v)
    return f * 100.0 if f <= 1.0 else f


def clade_seed_fraction(clade, seeds_set: Set[str]) -> float:
    """Fraction of clade terminals that are Ecuador core seeds."""
    terminals = {t.name for t in clade.get_terminals() if t.name}
    if not terminals:
        return 0.0
    return len(terminals & seeds_set) / len(terminals)


def mrca_candidates(
    tree,
    seeds: Iterable[str],
    candidates: Iterable[str],
    min_support: float,
    exclude: Set[str],
) -> List[Tuple[str, float, float, str]]:
    """Score each candidate by its MRCA with the seed set.

    Returns list of (taxon, seed_fraction, patristic_dist_to_nearest_seed, nearest_seed)
    sorted primary by seed_fraction DESC, tiebreak by patristic distance ASC.

    Candidates are excluded when:
    - Their shared MRCA with the seed set is the tree root (phylogenetically outside seed clade).
    - Their shared MRCA has support < min_support.
    """
    tree_tips = get_terminals(tree)
    seed_list = [s for s in seeds if s in tree_tips]
    if not seed_list:
        return []
    seeds_set = set(seed_list)

    scored: List[Tuple[str, float, float, str]] = []
    for cand in candidates:
        if cand in exclude or cand not in tree_tips:
            continue
        try:
            mrca = tree.common_ancestor([cand] + seed_list)
        except Exception:
            continue
        # Exclude candidates whose only shared ancestor with the seeds is the root
        if mrca is tree.root:
            continue
        support_raw = getattr(mrca, "confidence", None)
        if support_raw is not None and normalize_support(support_raw) < min_support:
            continue
        sf = clade_seed_fraction(mrca, seeds_set)
        min_dist = float("inf")
        nearest_seed = seed_list[0]
        for seed in seed_list:
            try:
                d = tree.distance(cand, seed)
                if d < min_dist:
                    min_dist = d
                    nearest_seed = seed
            except Exception:
                continue
        if min_dist == float("inf"):
            continue
        scored.append((cand, sf, min_dist, nearest_seed))

    scored.sort(key=lambda x: (-x[1], x[2]))
    return scored


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


def read_date_map(path: str) -> Dict[str, str]:
    """Returns accession -> normalized ISO YYYY-MM-DD from the context metadata TSV."""
    df = pd.read_csv(path, sep=None, engine="python")
    if "collection_date" not in df.columns or "accession" not in df.columns:
        return {}
    out: Dict[str, str] = {}
    for _, row in df.iterrows():
        acc = str(row["accession"]).strip()
        raw_date = str(row.get("collection_date", "")).strip()
        parsed = parse_collection_date(raw_date)
        if acc and acc != "nan":
            out[acc] = parsed or ""
    return out


def taxon_ym_bucket(label: str, date_map: Dict[str, str]) -> str:
    """Returns YYYY-MM bucket for a taxon label using its GenBank accession, or 'UNKNOWN'."""
    match = ACCESSION_RE.search(label)
    if match:
        iso = date_map.get(match.group(1), "")
        if iso and len(iso) >= 7:
            return iso[:7]
    return "UNKNOWN"


def nearest_candidates(
    tree,
    seeds: Iterable[str],
    candidates: Iterable[str],
    n_take: int,
    exclude: Set[str],
) -> List[Tuple[str, float, str]]:
    """Distance-based selector used for out-of-clade groups (usa_distal, american_anchor)
    where MRCA with Ecuador seeds is expected to be the root."""
    tree_tips = get_terminals(tree)
    seed_list = [s for s in seeds if s in tree_tips]
    if not seed_list:
        return []
    scored: List[Tuple[str, float, str]] = []
    for cand in candidates:
        if cand in exclude or cand not in tree_tips:
            continue
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


def select_regional_context(
    candidates_scored: List[Tuple[str, float, float, str]],
    country_map: Dict[str, str],
    date_map: Dict[str, str],
    n_per_country: int,
    n_total: int,
    max_per_country_month: int,
) -> List[Tuple[str, float, str, str]]:
    """Select up to n_total regional context candidates with country/month spread.

    Pass 1 – diversity: fill (country × YYYY-MM) buckets, respecting:
      - n_per_country:         max samples per country.
      - max_per_country_month: max samples per (country, YYYY-MM) bucket.

    Pass 2 – fill by MRCA score while still respecting n_per_country,
    relaxing only max_per_country_month.

    Returns list of (taxon, dist, nearest_seed, reason).
    """
    selected: List[Tuple[str, float, str, str]] = []
    selected_ids: Set[str] = set()
    country_count: Dict[str, int] = {}
    country_month_count: Dict[Tuple[str, str], int] = {}

    # Pass 1: diversity (country + country/month caps)
    for taxon, sf, dist, seed in candidates_scored:
        if len(selected) >= n_total:
            break
        if taxon in selected_ids:
            continue
        country = taxon_country(taxon, country_map)
        ym = taxon_ym_bucket(taxon, date_map)
        if country_count.get(country, 0) >= n_per_country:
            continue
        if country_month_count.get((country, ym), 0) >= max_per_country_month:
            continue
        selected.append((taxon, dist, seed, "country_month_diversity"))
        selected_ids.add(taxon)
        country_count[country] = country_count.get(country, 0) + 1
        country_month_count[(country, ym)] = country_month_count.get((country, ym), 0) + 1

    # Pass 2: fill by MRCA score, still respecting n_per_country
    for taxon, sf, dist, seed in candidates_scored:
        if len(selected) >= n_total:
            break
        if taxon in selected_ids:
            continue
        country = taxon_country(taxon, country_map)
        if country_count.get(country, 0) >= n_per_country:
            continue
        selected.append((taxon, dist, seed, "mrca_score_fill"))
        selected_ids.add(taxon)
        country_count[country] = country_count.get(country, 0) + 1

    return selected
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build main BEAST panel using MRCA-clade selection with country/month spread"
    )
    parser.add_argument("--tree", required=True)
    parser.add_argument("--context-metadata", required=True)
    parser.add_argument("--panel-main-out", required=True)
    parser.add_argument("--audit-out", required=True)
    parser.add_argument(
        "--country-month-audit-out",
        default=None,
        help="Optional path for country/month coverage TSV (separate from main audit)",
    )
    args = parser.parse_args()

    tree = read_tree(args.tree)
    country_map = read_country_map(args.context_metadata)
    date_map = read_date_map(args.context_metadata)
    tree_tips = get_terminals(tree)
    panel_profile = BEAST_PANEL_PROFILE

    # Use all tree tips as candidates (concat tree already guarantees 8-segment completeness)
    complete_ids = {canonical_tip(t) for t in tree_tips}
    complete_filter_source = "tree_tips_all"

    flu_tip_map: Dict[str, str] = {}
    for tip in tree_tips:
        if tip.startswith("Flu-"):
            flu_tip_map[flu_base_id(tip)] = tip

    core_available = []
    core_missing = []
    for core in ECUADOR_CORE:
        tip_label = flu_tip_map.get(core, core)
        if tip_label in tree_tips:
            core_available.append(tip_label)
        else:
            core_missing.append(core)

    seed_main = list(core_available)

    # ── Regional context: MRCA-scored + country/month spread ──────────────────
    regional_candidates = [
        t
        for t in tree_tips
        if canonical_tip(t) in complete_ids and is_regional_context(t)
    ]
    regional_scored = mrca_candidates(
        tree=tree,
        seeds=seed_main,
        candidates=regional_candidates,
        min_support=panel_profile["min_mrca_support"],
        exclude=set(),
    )
    selected_context = select_regional_context(
        candidates_scored=regional_scored,
        country_map=country_map,
        date_map=date_map,
        n_per_country=panel_profile["n_per_country"],
        n_total=panel_profile["n_total"],
        max_per_country_month=panel_profile["max_per_country_month"],
    )

    # ── USA distal: distance-based (intentionally out-of-clade, MRCA would be root) ──
    usa_distal_candidates = [
        t for t in tree_tips if canonical_tip(t) in complete_ids and is_usa_distal(t)
    ]
    selected_usa_distal = nearest_candidates(
        tree=tree,
        seeds=seed_main,
        candidates=usa_distal_candidates,
        n_take=panel_profile["usa_distal_quota"],
        exclude=set(),
    )

    # ── American anchor: distance-based (same reason as usa_distal) ──────────
    american_anchor_candidates = [
        t for t in tree_tips if canonical_tip(t) in complete_ids and is_american_anchor(t)
    ]
    selected_american_anchor: List[Tuple[str, float, str]] = []

    # Force the requested accession into american anchor regardless of distance rank.
    forced_anchor_accession = panel_profile["forced_american_anchor_accession"]
    forced_anchor = next(
        (t for t in american_anchor_candidates if forced_anchor_accession in t),
        None,
    )
    if forced_anchor:
        forced_pick = nearest_candidates(
            tree=tree,
            seeds=seed_main,
            candidates=[forced_anchor],
            n_take=1,
            exclude=set(),
        )
        if forced_pick:
            selected_american_anchor.extend(forced_pick)

    selected_american_anchor.extend(
        nearest_candidates(
            tree=tree,
            seeds=seed_main,
            candidates=american_anchor_candidates,
            n_take=panel_profile["additional_american_anchor_quota"],
            exclude={taxon for taxon, _, _ in selected_american_anchor},
        )
    )

    # ── Assemble panel rows ───────────────────────────────────────────────────
    panel_main_rows: List[Tuple[str, str, str, float]] = []
    for t in core_available:
        panel_main_rows.append((t, "ecuador_core", "main_cluster", None))
    for taxon, dist, seed, reason in selected_context:
        panel_main_rows.append((taxon, "regional_context", reason, dist))
    for taxon, dist, seed in selected_usa_distal:
        panel_main_rows.append((taxon, "usa_distal", "distance_to_ecuador_core", dist))
    for taxon, dist, seed in selected_american_anchor:
        panel_main_rows.append((taxon, "american_anchor", "distance_to_ecuador_core", dist))

    write_panel(args.panel_main_out, panel_main_rows)

    # ── Main audit ────────────────────────────────────────────────────────────
    audit_dir = os.path.dirname(args.audit_out)
    if audit_dir:
        os.makedirs(audit_dir, exist_ok=True)
    with open(args.audit_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["complete_filter_source", complete_filter_source])
        writer.writerow(["tree_tips_total", len(tree_tips)])
        writer.writerow(["panel_profile", panel_profile["label"]])
        writer.writerow(["ecuador_core_configured", len(ECUADOR_CORE)])
        writer.writerow(["core_available", len(core_available)])
        writer.writerow(["core_missing", len(core_missing)])
        writer.writerow(["regional_candidates_total", len(regional_candidates)])
        writer.writerow(["regional_candidates_passed_mrca", len(regional_scored)])
        writer.writerow([
            "regional_candidates_filtered_mrca",
            len(regional_candidates) - len(regional_scored),
        ])
        writer.writerow(["panel_main_total", len(panel_main_rows)])
        writer.writerow(["min_mrca_support", panel_profile["min_mrca_support"]])
        writer.writerow(["max_per_country_month", panel_profile["max_per_country_month"]])
        writer.writerow(["n_per_country", panel_profile["n_per_country"]])
        writer.writerow(["n_total", panel_profile["n_total"]])
        writer.writerow(["usa_distal_quota", panel_profile["usa_distal_quota"]])
        writer.writerow([
            "additional_american_anchor_quota",
            panel_profile["additional_american_anchor_quota"],
        ])
        writer.writerow([
            "forced_american_anchor_accession",
            forced_anchor_accession,
        ])
        writer.writerow([
            "regional_blacklist_tokens",
            ";".join(REGIONAL_BLACKLIST_TOKENS),
        ])
        writer.writerow([
            "selected_country_month_diversity",
            sum(1 for x in selected_context if x[3] == "country_month_diversity"),
        ])
        writer.writerow([
            "selected_mrca_score_fill",
            sum(1 for x in selected_context if x[3] == "mrca_score_fill"),
        ])
        writer.writerow(["selected_usa_distal", len(selected_usa_distal)])
        writer.writerow(["selected_american_anchor", len(selected_american_anchor)])
        writer.writerow([
            "forced_american_anchor_present",
            any(
                forced_anchor_accession
                and forced_anchor_accession in taxon
                for taxon, _, _ in selected_american_anchor
            ),
        ])
        country_counts: Dict[str, int] = {}
        for taxon, _, _, _ in selected_context:
            c = taxon_country(taxon, country_map)
            country_counts[c] = country_counts.get(c, 0) + 1
        for c, n in sorted(country_counts.items()):
            writer.writerow([f"country_count::{c}", n])
        if core_missing:
            writer.writerow(["core_missing_ids", ";".join(core_missing)])

    # ── Country/month coverage audit (optional separate file) ─────────────────
    cm_out = args.country_month_audit_out
    if cm_out:
        cm_dir = os.path.dirname(cm_out)
        if cm_dir:
            os.makedirs(cm_dir, exist_ok=True)
        coverage: Dict[Tuple[str, str], int] = {}
        for taxon, _, _, _ in selected_context:
            country = taxon_country(taxon, country_map)
            ym = taxon_ym_bucket(taxon, date_map)
            key = (country, ym)
            coverage[key] = coverage.get(key, 0) + 1
        with open(cm_out, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["country", "month", "n_selected"])
            for (country, month), n in sorted(coverage.items()):
                writer.writerow([country, month, n])


if __name__ == "__main__":
    main()
