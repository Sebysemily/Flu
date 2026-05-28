#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
import heapq
from io import StringIO
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
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read().strip()
    return Phylo.read(StringIO(text), "newick")

def get_terminals(tree) -> Set[str]:
    return {t.name for t in tree.get_terminals() if t.name}

def calculate_distances_to_core(tree, core_names: Set[str]) -> Dict[str, float]:
    adj = {}
    name_to_clade = {}
    
    def traverse(clade):
        if clade.name:
            name_to_clade[clade.name] = clade
        if clade not in adj:
            adj[clade] = []
        for child in clade.clades:
            dist = child.branch_length if child.branch_length is not None else 0.0
            if child not in adj:
                adj[child] = []
            adj[clade].append((child, dist))
            adj[child].append((clade, dist))
            traverse(child)
            
    traverse(tree.root)
    
    pq = []
    visited = {}
    
    for name in core_names:
        if name in name_to_clade:
            clade = name_to_clade[name]
            visited[clade] = 0.0
            heapq.heappush(pq, (0.0, id(clade), clade))
            
    while pq:
        d, _, u = heapq.heappop(pq)
        if d > visited.get(u, float("inf")):
            continue
        for v, weight in adj.get(u, []):
            new_d = d + weight
            if new_d < visited.get(v, float("inf")):
                visited[v] = new_d
                heapq.heappush(pq, (new_d, id(v), v))
                
    tip_distances = {}
    for name, clade in name_to_clade.items():
        if clade in visited:
            tip_distances[name] = visited[clade]
            
    return tip_distances

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


CORE_ROLES = {"flu_costa", "flu_sierra", "flu_epi_isl"}

def write_panel(path: str, rows: List[Tuple[str, str, str, float]]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxon", "role", "lineage", "distance_to_seed"])
        for taxon, role, lineage, dist in rows:
            writer.writerow([taxon, role, lineage, "" if dist is None else f"{dist:.6f}"])

def main() -> None:
    parser = argparse.ArgumentParser(description="Consolidated script to select, subset and prune tree to BEAST panel by distance.")
    parser.add_argument("--alignment", required=True, help="Input alignment to subset")
    parser.add_argument("--tree", required=True, help="Input Newick tree to prune")
    parser.add_argument(
        "--panel-main-out",
        default=None,
        help="Optional TSV audit of panel taxa (role/lineage). Omit to rely on --out-tree tips only.",
    )
    parser.add_argument("--out-alignment", required=True, help="Output subset FASTA alignment")
    parser.add_argument("--out-tree", required=True, help="Output pruned Newick tree")
    parser.add_argument("--audit-out", required=True, help="Output subset audit TSV")
    parser.add_argument("--distance-audit-out", required=True, help="Output distance audit TSV")
    parser.add_argument("--discarded-csv", default=None, help="Output CSV for discarded sequences")
    parser.add_argument("--country-month-audit-out", default=None, help="Output country/month audit TSV")
    parser.add_argument("--n-closest", type=int, default=200, help="Target number of closest candidates to select")
    parser.add_argument("--max-distance", type=float, default=0.08, help="Maximum patristic distance threshold")
    parser.add_argument("--protect-anchors-per-month", type=int, default=2, help="Number of oldest American anchors to protect per month up to 2022-05")
    parser.add_argument("--protect-regional-per-month", type=int, default=3, help="Number of regional context sequences to protect per month closest to core")
    parser.add_argument("--metadata", required=True, help="Path to metadata/H5N1_context.csv")
    args = parser.parse_args()

    panel_meta = load_panel_metadata(args.metadata)

    tree = read_tree(args.tree)
    tree_tips = get_terminals(tree)

    def role_for_tip(tip: str) -> str:
        return panel_meta.get(tip, {}).get("expected_role", "")

    # 1. Identify groups from unified metadata
    core_tips = sorted([t for t in tree_tips if role_for_tip(t) in CORE_ROLES])
    american_anchor_tips = sorted([t for t in tree_tips if role_for_tip(t) == "american_anchor"])
    regional_context_tips = sorted([t for t in tree_tips if role_for_tip(t) == "regional_context"])

    core_costa = [t for t in core_tips if role_for_tip(t) == "flu_costa"]
    core_sierra = [t for t in core_tips if role_for_tip(t) == "flu_sierra"]
    core_epi = [t for t in core_tips if role_for_tip(t) == "flu_epi_isl"]

    print(f"Loaded tree with {len(tree_tips)} tips.")
    print(f"Found {len(core_tips)} Ecuador core tips (metadata roles).")
    print(f"  flu_costa: {len(core_costa)}, flu_sierra: {len(core_sierra)}, flu_epi_isl: {len(core_epi)}")
    print(f"Found {len(american_anchor_tips)} American anchor tips.")
    print(f"Found {len(regional_context_tips)} Regional context tips.")

    tips_missing_role = sorted(
        t for t in tree_tips
        if t not in core_tips and t not in american_anchor_tips and t not in regional_context_tips
    )
    if tips_missing_role:
        print(f"Warning: {len(tips_missing_role)} tips have no expected_role in metadata (excluded from panel logic).")

    # 2. Calculate patristic distance to closest core sequence for all anchor and regional context tips using Dijkstra's algorithm
    tip_distances = calculate_distances_to_core(tree, set(core_tips))

    # Group anchors by month and protect the closest to core per month up to May 2022
    anchors_by_month = {}
    for tip in american_anchor_tips:
        date_str = parse_date_for_tip(tip, panel_meta)
        match = re.search(r"(\d{4})-(\d{2})", date_str)
        month_key = f"{match.group(1)}-{match.group(2)}" if match else "UNKNOWN"

        if month_key not in anchors_by_month:
            anchors_by_month[month_key] = []
        anchors_by_month[month_key].append(tip)

    protected_anchors = set()
    for month_key, month_tips in anchors_by_month.items():
        if month_key != "UNKNOWN" and month_key <= "2022-05":
            valid_tips = [t for t in month_tips if t in tip_distances]
            valid_tips.sort(key=lambda t: tip_distances[t])
            limit = min(args.protect_anchors_per_month, len(valid_tips))
            for i in range(limit):
                protected_anchors.add(valid_tips[i])

    if len(protected_anchors) > 0:
        print(f"Protected {len(protected_anchors)} closest American anchor tips per month up to 2022-05: {sorted(protected_anchors)}")

    # Group regional contexts by month and protect the closest to core per month
    regional_by_month = {}
    for tip in regional_context_tips:
        date_str = parse_date_for_tip(tip, panel_meta)
        match = re.search(r"(\d{4})-(\d{2})", date_str)
        month_key = f"{match.group(1)}-{match.group(2)}" if match else "UNKNOWN"

        if month_key not in regional_by_month:
            regional_by_month[month_key] = []
        regional_by_month[month_key].append(tip)

    protected_regional = set()
    for month_key, month_tips in regional_by_month.items():
        if month_key != "UNKNOWN":
            valid_tips = [t for t in month_tips if t in tip_distances]
            valid_tips.sort(key=lambda t: tip_distances[t])
            limit = min(args.protect_regional_per_month, len(valid_tips))
            for i in range(limit):
                protected_regional.add(valid_tips[i])

    if len(protected_regional) > 0:
        print(f"Protected {len(protected_regional)} closest Regional context tips per month: {sorted(protected_regional)}")

    # 3. Identify remaining candidates (non-protected anchors and regional contexts)
    candidates = []
    for tip in american_anchor_tips:
        if tip not in protected_anchors and tip in tip_distances:
            candidates.append((tip, "american_anchor", tip_distances[tip]))
    for tip in regional_context_tips:
        if tip not in protected_regional and tip in tip_distances:
            candidates.append((tip, "regional_context", tip_distances[tip]))

    # 4. Filter by max distance and take top N closest
    passed_max_dist = [c for c in candidates if c[2] <= args.max_distance]
    passed_max_dist.sort(key=lambda x: x[2]) # sort by distance ascending
    
    selected_candidates = passed_max_dist[:args.n_closest]

    print(f"Total scored context tips: {len(tip_distances)}")
    print(f"Candidates within max_distance ({args.max_distance}): {len(passed_max_dist)}")
    print(f"Selected top {len(selected_candidates)} closest candidates.")

    # 5. Build final panel rows
    panel_rows = []
    
    # Core Ecuador — preserve flu_costa / flu_sierra from H5N1_context
    for tip in core_tips:
        panel_rows.append((tip, role_for_tip(tip), "main_cluster", None))
        
    # Protected Anchors
    for tip in sorted(protected_anchors):
        panel_rows.append((tip, "american_anchor", "closest_protected", tip_distances.get(tip, None)))
        
    # Protected Regional Contexts
    for tip in sorted(protected_regional):
        panel_rows.append((tip, "regional_context", "closest_protected", tip_distances.get(tip, None)))
        
    # Selected Context candidates
    for tip, role, dist in selected_candidates:
        panel_rows.append((tip, role, "distance_to_ecuador_core", dist))

    if args.panel_main_out:
        write_panel(args.panel_main_out, panel_rows)
        print(f"Wrote optional panel audit TSV to {args.panel_main_out} ({len(panel_rows)} taxa).")
    print(f"Final panel: {len(panel_rows)} taxa (tip list = pruned tree output).")

    # Set of final panel taxon names
    panel_taxa_set = {row[0] for row in panel_rows}

    # 7. Subset alignment
    print(f"Subsetting alignment {args.alignment}...")
    records = list(SeqIO.parse(args.alignment, "fasta"))
    kept_records = [r for r in records if r.id in panel_taxa_set]
    
    out_aln_dir = os.path.dirname(args.out_alignment)
    if out_aln_dir:
        os.makedirs(out_aln_dir, exist_ok=True)
    SeqIO.write(kept_records, args.out_alignment, "fasta")
    print(f"Wrote subset alignment to {args.out_alignment} ({len(kept_records)} records).")

    # 8. Prune tree
    print("Pruning tree...")
    pruned_tree = deepcopy(tree)
    for tip in list(pruned_tree.get_terminals()):
        if tip.name not in panel_taxa_set:
            try:
                pruned_tree.prune(tip)
            except ValueError:
                # Ignore edge-case if node already pruned through parent collapse
                pass
    
    out_tree_dir = os.path.dirname(args.out_tree)
    if out_tree_dir:
        os.makedirs(out_tree_dir, exist_ok=True)
    Phylo.write(pruned_tree, args.out_tree, "newick")
    print(f"Wrote pruned tree to {args.out_tree} ({len(pruned_tree.get_terminals())} terminals).")

    # 9. Write distance audit
    out_dist_dir = os.path.dirname(args.distance_audit_out)
    if out_dist_dir:
        os.makedirs(out_dist_dir, exist_ok=True)
    with open(args.distance_audit_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["tree_tips_total", len(tree_tips)])
        writer.writerow(["ecuador_core_tips", len(core_tips)])
        writer.writerow(["flu_costa_tips", len(core_costa)])
        writer.writerow(["flu_sierra_tips", len(core_sierra)])
        writer.writerow(["flu_epi_isl_tips", len(core_epi)])
        writer.writerow(["tips_missing_metadata_role", len(tips_missing_role)])
        writer.writerow(["american_anchors_total", len(american_anchor_tips)])
        writer.writerow(["american_anchors_protected", len(protected_anchors)])
        writer.writerow(["regional_context_total", len(regional_context_tips)])
        writer.writerow(["regional_context_protected", len(protected_regional)])
        writer.writerow(["candidates_total", len(candidates)])
        writer.writerow(["candidates_passed_distance_cutoff", len(passed_max_dist)])
        writer.writerow(["candidates_selected", len(selected_candidates)])
        writer.writerow(["panel_main_total", len(panel_rows)])
        writer.writerow(["max_distance_cutoff", args.max_distance])
        writer.writerow(["n_closest_target", args.n_closest])

    # 10. Write subset audit
    out_audit_dir = os.path.dirname(args.audit_out)
    if out_audit_dir:
        os.makedirs(out_audit_dir, exist_ok=True)
        
    aln_ids = {r.id for r in records}
    missing_in_tree = sorted(panel_taxa_set - tree_tips)
    missing_in_alignment = sorted(panel_taxa_set - aln_ids)

    with open(args.audit_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["requested_taxa", len(panel_taxa_set)])
        writer.writerow(["kept_alignment", len(kept_records)])
        writer.writerow(["kept_tree", len(pruned_tree.get_terminals())])
        writer.writerow(["missing_in_tree", ";".join(missing_in_tree)])
        writer.writerow(["missing_in_alignment", ";".join(missing_in_alignment)])

    # 11. Write country/month audit if requested
    if args.country_month_audit_out:
        cm_dir = os.path.dirname(args.country_month_audit_out)
        if cm_dir:
            os.makedirs(cm_dir, exist_ok=True)
        with open(args.country_month_audit_out, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["country", "month", "n_selected"])
            coverage = {}
            for tip, role, dist in selected_candidates:
                date_str = parse_date_for_tip(tip, panel_meta)
                ym = date_str[:7] if len(date_str) >= 7 else "UNKNOWN"
                country = panel_meta.get(tip, {}).get("country", "UNKNOWN")
                if country == "UNKNOWN":
                    parts = tip.split("/")
                    country = parts[1] if len(parts) >= 2 else "UNKNOWN"
                coverage[(country, ym)] = coverage.get((country, ym), 0) + 1
            for (country, month), n in sorted(coverage.items()):
                writer.writerow([country, month, n])

    # 12. Write discarded CSV if requested
    if args.discarded_csv:
        disc_dir = os.path.dirname(args.discarded_csv)
        if disc_dir:
            os.makedirs(disc_dir, exist_ok=True)
            
        cand_dist = tip_distances
        passed_set = {cand for cand, _, _ in passed_max_dist}
        
        discarded_rows = []
        for tip in sorted(tree_tips):
            if tip in panel_taxa_set:
                continue
            
            dist_val = cand_dist.get(tip, None)
            if tip not in cand_dist:
                reason = "no path to core sequence in tree"
            elif tip not in passed_set:
                reason = f"exceeded max_distance (distance: {dist_val:.6f} > {args.max_distance})"
            else:
                reason = f"not in top {args.n_closest} closest (distance: {dist_val:.6f})"
                
            discarded_rows.append([tip, reason, "" if dist_val is None else f"{dist_val:.6f}"])
            
        with open(args.discarded_csv, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=",")
            writer.writerow(["taxon", "reason", "distance_to_seed"])
            for row in discarded_rows:
                writer.writerow(row)
        print(f"Wrote discarded sequences to {args.discarded_csv} ({len(discarded_rows)} sequences).")

if __name__ == "__main__":
    main()
