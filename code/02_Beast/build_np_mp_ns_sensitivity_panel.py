#!/usr/bin/env python3
import argparse
import csv
import os
from io import StringIO
from typing import Dict, List, Set

from Bio import Phylo


def read_tree_tips(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as handle:
        tree = Phylo.read(StringIO(handle.read().strip()), "newick")
    return [terminal.name for terminal in tree.get_terminals() if terminal.name]


def is_usa_distal(taxon: str) -> bool:
    return "__usa_distal" in taxon


def is_usa_neighbor(taxon: str) -> bool:
    return "__usa_neighbor" in taxon


def is_american_anchor(taxon: str) -> bool:
    return "__american_anchor" in taxon


def is_non_american_anchor(taxon: str) -> bool:
    return "_anchor" in taxon and not is_american_anchor(taxon)


def role_for(taxon: str, selected_american_anchors: Set[str]) -> str:
    if taxon in selected_american_anchors:
        return "american_anchor"
    if "__regional_context" in taxon:
        return "regional_context"
    if taxon.startswith("Flu-"):
        return "ecuador"
    return "context"


def lineage_for(
    taxon: str,
    forced_american_anchor: str,
    selected_american_anchors: Set[str],
) -> str:
    if taxon in selected_american_anchors:
        if forced_american_anchor and forced_american_anchor in taxon:
            return "forced_american_anchor"
        return "limited_american_anchor"
    return "included_by_group_filter"


def write_panel(
    path: str,
    tips: List[str],
    selected_american_anchors: Set[str],
    forced_american_anchor: str,
) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxon", "role", "lineage", "distance_to_seed"])
        for taxon in tips:
            writer.writerow(
                [
                    taxon,
                    role_for(taxon, selected_american_anchors),
                    lineage_for(taxon, forced_american_anchor, selected_american_anchors),
                    "",
                ]
            )


def write_audit(path: str, metrics: Dict[str, object]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, value])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build an NP+MP+NS RTT sensitivity panel by group filters, "
            "without MRCA/clade pruning."
        )
    )
    parser.add_argument("--tree", required=True)
    parser.add_argument("--panel-out", required=True)
    parser.add_argument("--audit-out", required=True)
    parser.add_argument(
        "--forced-american-anchor-accession",
        default="OQ968009",
    )
    parser.add_argument(
        "--american-anchor-total",
        type=int,
        default=6,
        help="Total american anchors to keep, including the forced accession.",
    )
    args = parser.parse_args()

    if args.american_anchor_total < 1:
        raise ValueError("--american-anchor-total must be at least 1")

    tips = read_tree_tips(args.tree)
    american_anchors = [tip for tip in tips if is_american_anchor(tip)]
    if not american_anchors:
        print("Warning: No American anchors found in the tree. Skipping American anchor selection.")
        selected_american_anchors = set()
        forced_anchor = None
    else:
        forced_anchor = next(
            (
                tip
                for tip in american_anchors
                if args.forced_american_anchor_accession in tip
            ),
            None,
        )
        if forced_anchor is None:
            raise ValueError(
                "Forced american anchor accession "
                f"{args.forced_american_anchor_accession!r} was not found in the tree"
            )

        additional_anchors = [
            tip
            for tip in sorted(american_anchors)
            if tip != forced_anchor
        ][: args.american_anchor_total - 1]
        selected_american_anchors = set([forced_anchor, *additional_anchors])

    included: List[str] = []
    excluded_counts = {
        "excluded_usa_distal": 0,
        "excluded_usa_neighbor": 0,
        "excluded_non_american_anchor": 0,
        "excluded_extra_american_anchor": 0,
    }

    for tip in tips:
        if is_usa_distal(tip):
            excluded_counts["excluded_usa_distal"] += 1
            continue
        if is_usa_neighbor(tip):
            excluded_counts["excluded_usa_neighbor"] += 1
            continue
        if is_american_anchor(tip):
            if tip in selected_american_anchors:
                included.append(tip)
            else:
                excluded_counts["excluded_extra_american_anchor"] += 1
            continue
        if is_non_american_anchor(tip):
            excluded_counts["excluded_non_american_anchor"] += 1
            continue
        included.append(tip)

    write_panel(
        args.panel_out,
        included,
        selected_american_anchors,
        args.forced_american_anchor_accession,
    )
    write_audit(
        args.audit_out,
        {
            "tree_tips_total": len(tips),
            "panel_taxa_total": len(included),
            "american_anchor_total_available": len(american_anchors),
            "american_anchor_total_kept": len(selected_american_anchors),
            "american_anchor_total_limit": args.american_anchor_total,
            "forced_american_anchor_accession": args.forced_american_anchor_accession,
            "forced_american_anchor_taxon": forced_anchor,
            "additional_american_anchor_selection": "lexicographic_taxon_order",
            "selected_american_anchors": ";".join(sorted(selected_american_anchors)),
            **excluded_counts,
        },
    )


if __name__ == "__main__":
    main()
