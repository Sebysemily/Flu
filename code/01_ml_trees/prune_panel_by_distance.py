#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from typing import Dict, Iterable, List, Set, Tuple
from copy import deepcopy

import pandas as pd
from Bio import Phylo, SeqIO

# Allow importing date_normalization from the parent code/ directory
_CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)
from date_normalization import parse_collection_date  # noqa: E402

def read_tree(path: str):
    return Phylo.read(path, "newick")

def get_terminals(tree) -> Set[str]:
    return {t.name for t in tree.get_terminals() if t.name}

def load_panel_metadata(path: str) -> Dict[str, Dict[str, str]]:
    """Load deduped panel metadata: file_name -> {expected_role, collection_date}."""
    meta = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("file_name") or "").strip()
            if not name or name in meta:
                continue
            meta[name] = {
                "expected_role": (row.get("expected_role") or "").strip(),
                "collection_date": (row.get("collection_date") or "").strip(),
                "country": (row.get("country") or "").strip(),
            }
    return meta

def parse_date_for_tip(label: str, meta: Dict[str, Dict[str, str]]) -> str:
    if label in meta and meta[label].get("collection_date"):
        return meta[label]["collection_date"]
    parts = label.split("/")
    if len(parts) >= 3:
        return parts[-1].strip()
    return "UNKNOWN"

def write_panel(path: str, rows: List[Tuple[str, str, str, float]]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxon", "role", "lineage", "distance_to_mrca"])
        for taxon, role, lineage, dist in rows:
            writer.writerow([taxon, role, lineage, "" if dist is None else f"{dist:.6f}"])

def main() -> None:
    parser = argparse.ArgumentParser(description="Subset tree to BEAST panel by distance to MRCA.")
    parser.add_argument("--tree", required=True, help="Input Newick tree to measure distances")
    parser.add_argument("--out-context", required=True, help="Output text file with selected context taxa")
    parser.add_argument("--audit-out", required=True, help="Output subset audit TSV")
    parser.add_argument("--distance-audit-out", required=True, help="Output distance audit TSV")
    parser.add_argument("--discarded-csv", default=None, help="Output CSV for discarded sequences")
    parser.add_argument("--metadata", required=True, help="Path to metadata/H5N1_context.csv")
    parser.add_argument("--outgroup-root-sample", required=True, help="Oldest root sample to stop temporal protection")
    
    args = parser.parse_args()

    panel_meta = load_panel_metadata(args.metadata)
    tree = read_tree(args.tree)
    tree_tips = get_terminals(tree)

    def role_for_tip(tip: str) -> str:
        return panel_meta.get(tip, {}).get("expected_role", "")

    # 1. Identify core tips based on regex flu_.*
    core_tips = sorted([t for t in tree_tips if re.match(r"^flu_.*", role_for_tip(t), re.IGNORECASE)])
    american_anchor_tips = sorted([t for t in tree_tips if role_for_tip(t) == "american_anchor"])
    regional_context_tips = sorted([t for t in tree_tips if role_for_tip(t) == "regional_context"])

    print(f"Loaded tree with {len(tree_tips)} tips.")
    print(f"Found {len(core_tips)} core flu tips (metadata roles starting with flu_).")
    print(f"Found {len(american_anchor_tips)} American anchor tips.")
    print(f"Found {len(regional_context_tips)} Regional context tips.")

    # 2. Find MRCA of core tips and calculate distances
    core_clades = [tree.find_any(name=t) for t in core_tips if tree.find_any(name=t)]
    if not core_clades:
        raise ValueError("No core clades found in tree.")
    
    mrca = tree.common_ancestor(core_clades)
    print(f"Calculated MRCA for {len(core_clades)} core tips.")
    
    tip_distances = {}
    for tip in tree_tips:
        clade = tree.find_any(name=tip)
        if clade:
            try:
                # Path length between MRCA and the tip
                dist = tree.distance(mrca, clade)
                tip_distances[tip] = dist
            except Exception:
                pass

    # 3. Temporal Protection from MRCA clade (5 per month)
    def get_month(tip):
        date_str = parse_date_for_tip(tip, panel_meta)
        match = re.search(r"(\d{4})-(\d{2})", date_str)
        return f"{match.group(1)}-{match.group(2)}" if match else "UNKNOWN"

    anchors_by_month = {}
    protected_non_usa_anchors = set()
    
    for tip in american_anchor_tips:
        # Only consider american anchors that are within the MRCA clade!
        if tip in tip_distances:
            # Protect non-USA american anchors
            country = panel_meta.get(tip, {}).get("country", "")
            if country != "USA":
                protected_non_usa_anchors.add(tip)
            
            # Bucket by month for USA anchors
            m = get_month(tip)
            if m != "UNKNOWN":
                if m not in anchors_by_month:
                    anchors_by_month[m] = []
                anchors_by_month[m].append(tip)

    protected_temporal_anchors = set()
    for month_key, month_tips in anchors_by_month.items():
        valid_tips = [t for t in month_tips if t in tip_distances]
        valid_tips.sort(key=lambda t: tip_distances[t])
        limit = min(5, len(valid_tips))
        for i in range(limit):
            protected_temporal_anchors.add(valid_tips[i])

    print(f"Protected {len(protected_non_usa_anchors)} non-USA american anchors.")
    print(f"Protected {len(protected_temporal_anchors)} temporal anchors (up to 5 per month from MRCA clade).")

    # 4. Select remaining candidates to reach exactly 300 total tips
    protected_taxa_set = set(core_tips).union(protected_temporal_anchors).union(protected_non_usa_anchors)
    
    candidates = []
    for tip in tree_tips:
        if tip not in protected_taxa_set and tip in tip_distances:
            role = role_for_tip(tip)
            if role in ("american_anchor", "regional_context"):
                candidates.append((tip, role, tip_distances[tip]))
                
    # Sort candidates by distance (ascending)
    candidates.sort(key=lambda x: x[2])
    
    n_to_select = max(0, 300 - len(protected_taxa_set))
    selected_candidates = candidates[:n_to_select]
    print(f"Total protected taxa (core + protected anchors): {len(protected_taxa_set)}")
    print(f"Selected closest candidates to reach 300: {len(selected_candidates)} (out of {len(candidates)} candidates)")

    # 5. Build final panel rows
    panel_rows = []
    
    # Core Ecuador
    for tip in core_tips:
        panel_rows.append((tip, role_for_tip(tip), "core_flu", tip_distances.get(tip, None)))
            
    # Temporal Protected Anchors & Non-USA anchors
    for tip in protected_temporal_anchors.union(protected_non_usa_anchors):
        if tip not in core_tips:
            panel_rows.append((tip, "american_anchor", "temporal_or_non_usa", tip_distances.get(tip, None)))

    # Selected Candidates
    for tip, role, dist in selected_candidates:
        panel_rows.append((tip, role, "closest_candidates", dist))

    print(f"Final panel: {len(panel_rows)} taxa.")
    panel_taxa_set = {row[0] for row in panel_rows}

    # 6. Write context taxa list
    print(f"Writing context taxa list to {args.out_context}...")
    out_dir = os.path.dirname(args.out_context)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out_context, "w") as f:
        for tip in panel_taxa_set:
            if not tip.startswith("flu_") and tip not in core_tips:
                f.write(f"{tip}\n")
    print(f"Wrote context taxa to {args.out_context}.")

    # 7. Write distance audit
    out_dist_dir = os.path.dirname(args.distance_audit_out)
    if out_dist_dir:
        os.makedirs(out_dist_dir, exist_ok=True)
    with open(args.distance_audit_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["tree_tips_total", len(tree_tips)])
        writer.writerow(["core_flu_tips", len(core_tips)])
        writer.writerow(["american_anchors_total", len(american_anchor_tips)])
        writer.writerow(["american_anchors_temporal", len(protected_temporal_anchors)])
        writer.writerow(["regional_context_total", len(regional_context_tips)])
        writer.writerow(["panel_total", len(panel_rows)])

    # 9. Write subset audit
    out_audit_dir = os.path.dirname(args.audit_out)
    if out_audit_dir:
        os.makedirs(out_audit_dir, exist_ok=True)

    with open(args.audit_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["requested_taxa", len(panel_taxa_set)])
        writer.writerow(["core_taxa", len(core_tips)])
        writer.writerow(["context_taxa", len(panel_taxa_set) - len(core_tips)])

    # 9. Write discarded CSV if requested
    if args.discarded_csv:
        disc_dir = os.path.dirname(args.discarded_csv)
        if disc_dir:
            os.makedirs(disc_dir, exist_ok=True)
            
        discarded_rows = []
        for tip in sorted(tree_tips):
            if tip in panel_taxa_set:
                continue
            
            dist_val = tip_distances.get(tip, None)
            reason = f"exceeded distance threshold or no protection (distance: {dist_val:.6f})" if dist_val is not None else "no path to MRCA"
            discarded_rows.append([tip, reason, "" if dist_val is None else f"{dist_val:.6f}"])
            
        with open(args.discarded_csv, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=",")
            writer.writerow(["taxon", "reason", "distance_to_mrca"])
            for row in discarded_rows:
                writer.writerow(row)
        print(f"Wrote discarded sequences to {args.discarded_csv} ({len(discarded_rows)} sequences).")

if __name__ == "__main__":
    main()
