#!/usr/bin/env python3
import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from io import StringIO
from typing import Dict, Iterable, List, Optional, Set, Tuple

from Bio import Phylo


DEFAULT_TREES: List[Tuple[str, str]] = [
    ("PB2", "results/phylogeny/raxml/PB2/H5N1_PB2.raxml.supportTBE"),
    ("PB1", "results/phylogeny/raxml/PB1/H5N1_PB1.raxml.supportTBE"),
    ("PA", "results/phylogeny/raxml/PA/H5N1_PA.raxml.supportTBE"),
    ("HA", "results/phylogeny/raxml/HA/H5N1_HA.raxml.supportTBE"),
    ("NP", "results/phylogeny/raxml/NP/H5N1_NP.raxml.supportTBE"),
    ("NA", "results/phylogeny/raxml/NA/H5N1_NA.raxml.supportTBE"),
    ("MP", "results/phylogeny/raxml/MP/H5N1_MP.raxml.supportTBE"),
    ("NS", "results/phylogeny/raxml/NS/H5N1_NS.raxml.supportTBE"),
    ("concat", "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE"),
]

GROUPS = ["flu2024", "guayas2023", "pichincha0694"]
CODE_RE = re.compile(r"(Flu-\d{4})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize per-sample tree incongruence focused on group switching "
            "(TBE-based, standalone script)."
        )
    )
    parser.add_argument(
        "--tree",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Tree input in NAME=PATH format. Repeatable. Defaults include 8 segments + concat.",
    )
    parser.add_argument(
        "--out",
        default="incongruence/group_incongruence_summary.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--group-candidates-out",
        default="incongruence/group_candidates_by_baseline.tsv",
        help="Output TSV keyed by baseline group.",
    )
    parser.add_argument(
        "--incongruences-only-out",
        default="incongruence/group_incongruence_incongruences_only.tsv",
        help="Output TSV with only incongruent rows.",
    )
    parser.add_argument(
        "--min-support",
        type=float,
        default=70.0,
        help="Minimum TBE support to call a group change (default: 70).",
    )
    parser.add_argument(
        "--distance-margin",
        type=float,
        default=0.002,
        help="Minimum best-vs-second distance gap to accept low-support group switching.",
    )
    parser.add_argument(
        "--anchor-flu2024",
        default="",
        help="Comma-separated anchor IDs/codes for flu2024.",
    )
    parser.add_argument(
        "--anchor-guayas2023",
        default="Flu-0402,Flu-0403,Flu-0406,Flu-0407,Flu-0205,Flu-0206",
        help="Comma-separated anchor IDs/codes for guayas2023.",
    )
    parser.add_argument(
        "--anchor-pichincha0694",
        default="Flu-0694,Flu-0465",
        help="Comma-separated anchor IDs/codes for pichincha0694.",
    )
    parser.add_argument(
        "--secondary-min-votes",
        type=int,
        default=2,
        help=(
            "Minimum tree votes a non-baseline group must accumulate to trigger "
            "alert_in_other_group=yes (default: 2). Multi-tree signal requirement."
        ),
    )
    parser.add_argument(
        "--candidate-min-votes",
        type=int,
        default=1,
        help=(
            "Minimum tree votes a non-baseline group needs to include the sample "
            "as a secondary candidate in the group candidates TSV (default: 1). "
            "Samples with only one present tree can still appear as candidates."
        ),
    )
    return parser.parse_args()


def parse_tree_args(tree_args: List[str]) -> List[Tuple[str, str]]:
    if not tree_args:
        return list(DEFAULT_TREES)

    parsed = []
    for item in tree_args:
        if "=" not in item:
            raise ValueError(f"Invalid --tree value '{item}'. Expected NAME=PATH.")
        name, path = item.split("=", 1)
        name = name.strip()
        path = path.strip()
        if not name or not path:
            raise ValueError(f"Invalid --tree value '{item}'. Expected NAME=PATH.")
        parsed.append((name, path))
    return parsed


def read_tree(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        return Phylo.read(StringIO(handle.read().strip()), "newick")


def parent_map(tree) -> Dict:
    pmap = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            pmap[child] = clade
    return pmap


def extract_code(label: str) -> Optional[str]:
    if not label:
        return None
    match = CODE_RE.search(label)
    return match.group(1) if match else None


def parse_year(label: str) -> Optional[int]:
    parts = label.split("/")
    if not parts:
        return None
    try:
        return int(parts[-1])
    except ValueError:
        return None


def resolve_label(label: str, code_to_full_map: Dict[str, str]) -> str:
    if "/" in label:
        return label
    code = extract_code(label)
    if code and code in code_to_full_map:
        return code_to_full_map[code]
    return label


def is_full_id(label: str) -> bool:
    return label.startswith("Flu-") and "/" in label


def parse_anchor_values(raw: str) -> Set[str]:
    values: Set[str] = set()
    if not raw:
        return values
    for item in raw.split(","):
        value = item.strip()
        if value:
            values.add(value)
    return values


def median(values: List[float]) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def find_tip(tree, sample_full: str, sample_code: str):
    for tip in tree.get_terminals():
        tip_name = tip.name or ""
        if tip_name == sample_full:
            return tip

    for tip in tree.get_terminals():
        tip_name = tip.name or ""
        if tip_name == sample_code or tip_name.startswith(sample_code + "/"):
            return tip
    return None


def sister_context(tree, tip) -> Tuple[List[str], Optional[float]]:
    pmap = parent_map(tree)
    parent = pmap.get(tip)
    if parent is None:
        return [], None

    sisters: List[str] = []
    for child in parent.clades:
        if child is tip:
            continue
        for term in child.get_terminals():
            if term.name:
                sisters.append(term.name)

    support = parent.confidence if hasattr(parent, "confidence") else None
    return sisters, support


def clade_context(tree, tip, min_flu_size: int = 12) -> Tuple[List[str], Optional[float]]:
    pmap = parent_map(tree)
    node = pmap.get(tip)
    if node is None:
        return [], None

    while True:
        members = [t.name for t in node.get_terminals() if t.name and t is not tip]
        flu_members = [label for label in members if extract_code(label)]
        parent = pmap.get(node)
        if len(flu_members) >= min_flu_size or parent is None:
            support = node.confidence if hasattr(node, "confidence") else None
            return members, support
        node = parent


def normalize_support(raw_support: Optional[float]) -> Optional[float]:
    if raw_support is None:
        return None
    value = float(raw_support)
    if value <= 1.0:
        return value * 100.0
    return value


def build_sample_identity(
    trees: Dict[str, object],
    segment_tree_names: List[str],
) -> Tuple[List[Tuple[str, str]], Set[str], Dict[str, str]]:
    code_to_full_candidates: Dict[str, Counter] = defaultdict(Counter)

    for tree_name in segment_tree_names:
        tree = trees[tree_name]
        for tip in tree.get_terminals():
            if not tip.name:
                continue
            code = extract_code(tip.name)
            if not code:
                continue
            if is_full_id(tip.name):
                code_to_full_candidates[code][tip.name] += 1

    sample_rows: List[Tuple[str, str]] = []
    ambiguous_codes: Set[str] = set()

    all_codes: Set[str] = set()
    for tree_name in segment_tree_names:
        tree = trees[tree_name]
        for tip in tree.get_terminals():
            if tip.name:
                code = extract_code(tip.name)
                if code:
                    all_codes.add(code)

    for code in sorted(all_codes):
        candidates = code_to_full_candidates.get(code, Counter())
        if not candidates:
            sample_rows.append((code, code))
            continue

        ordered = sorted(candidates.items(), key=lambda kv: (-kv[1], kv[0]))
        chosen_full = ordered[0][0]
        if len(ordered) > 1:
            ambiguous_codes.add(code)
        sample_rows.append((chosen_full, code))

    code_to_full_map = {code: full for full, code in sample_rows}
    return sample_rows, ambiguous_codes, code_to_full_map


def tree_tip_index(tree) -> Tuple[Set[str], Dict[str, List[str]]]:
    labels: Set[str] = set()
    by_code: Dict[str, List[str]] = defaultdict(list)
    for tip in tree.get_terminals():
        if not tip.name:
            continue
        labels.add(tip.name)
        code = extract_code(tip.name)
        if code:
            by_code[code].append(tip.name)
    return labels, by_code


def code_in_group(code: str, full_label: str, group_anchor_values: Dict[str, Set[str]]) -> Dict[str, int]:
    score = {g: 0 for g in GROUPS}
    for group in GROUPS:
        if code in group_anchor_values[group] or full_label in group_anchor_values[group]:
            score[group] += 1
    return score


def classify_group_score(
    sample_full: str,
    sample_code: str,
    context_tips: List[str],
    group_anchor_values: Dict[str, Set[str]],
    code_to_full_map: Dict[str, str],
) -> str:
    score: Dict[str, int] = {g: 0 for g in GROUPS}

    resolved_sample = resolve_label(sample_full, code_to_full_map)
    sample_year = parse_year(resolved_sample)
    if sample_year is not None and sample_year >= 2024:
        score["flu2024"] += 6
    if "/Guayas/2023" in resolved_sample:
        score["guayas2023"] += 6
    if "/Pichincha/" in resolved_sample or sample_code == "Flu-0694":
        score["pichincha0694"] += 6

    self_score = code_in_group(sample_code, sample_full, group_anchor_values)
    for group, value in self_score.items():
        score[group] += value

    for raw_label in context_tips:
        label = resolve_label(raw_label, code_to_full_map)
        code = extract_code(label)
        if not code:
            continue

        anchor_score = code_in_group(code, label, group_anchor_values)
        for group, value in anchor_score.items():
            score[group] += value * 3

        if "/Guayas/2023" in label:
            score["guayas2023"] += 2
        if code == "Flu-0694" or "/Pichincha/" in label:
            score["pichincha0694"] += 2
        year = parse_year(label)
        if year is not None and year >= 2024:
            score["flu2024"] += 1

    max_group = max(score, key=score.get)
    max_value = score[max_group]
    if max_value <= 0:
        return "outside_main_groups"

    tied = [group for group in GROUPS if score[group] == max_value]
    if len(tied) > 1:
        return "outside_main_groups"

    return max_group


def anchors_for_group(
    group_name: str,
    group_anchor_values: Dict[str, Set[str]],
    sample_rows: List[Tuple[str, str]],
) -> Set[str]:
    anchors: Set[str] = set(group_anchor_values[group_name])
    if group_name != "flu2024":
        return anchors

    blocked = set(group_anchor_values.get("guayas2023", set())) | set(group_anchor_values.get("pichincha0694", set()))
    for sample_full, sample_code in sample_rows:
        if sample_full in blocked or sample_code in blocked:
            continue
        y = parse_year(sample_full)
        if y is not None and y >= 2024:
            anchors.add(sample_full)
    return anchors


def group_by_anchor_distance(
    tree,
    sample_tip_name: str,
    anchor_sets: Dict[str, Set[str]],
) -> Tuple[str, Dict[str, float], float]:
    labels, by_code = tree_tip_index(tree)
    group_distance_map: Dict[str, float] = {}

    for group_name, anchors in anchor_sets.items():
        dists: List[float] = []
        for anchor in anchors:
            code = extract_code(anchor)
            candidate_tips: List[str] = []
            if anchor in labels:
                candidate_tips.append(anchor)
            if code and code in labels:
                candidate_tips.append(code)
            if code and code in by_code:
                candidate_tips.extend(by_code[code])

            for cand in set(candidate_tips):
                if cand == sample_tip_name:
                    continue
                try:
                    dists.append(float(tree.distance(sample_tip_name, cand)))
                except Exception:
                    continue

        m = median(dists)
        if m is not None:
            group_distance_map[group_name] = m

    if not group_distance_map:
        return "outside_main_groups", {}, 0.0

    ranked = sorted(group_distance_map.items(), key=lambda kv: kv[1])
    best_group, best_dist = ranked[0]
    margin = 0.0
    if len(ranked) >= 2:
        margin = ranked[1][1] - best_dist
    return best_group, group_distance_map, margin


def serialize_map(items: Iterable[Tuple[str, str]]) -> str:
    return ";".join([f"{k}:{v}" for k, v in items])


def parse_kv_str(s: str) -> Dict[str, int]:
    """Parse a 'key:value;key:value' string into a dict with integer values."""
    if not s:
        return {}
    result: Dict[str, int] = {}
    for part in s.split(";"):
        if ":" in part:
            k, v = part.split(":", 1)
            try:
                result[k.strip()] = int(v.strip())
            except ValueError:
                pass
    return result


def main() -> None:
    args = parse_args()
    tree_pairs = parse_tree_args(args.tree)

    trees: Dict[str, object] = {}
    for name, path in tree_pairs:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Tree not found for {name}: {path}")
        trees[name] = read_tree(path)

    segment_tree_names = [name for name in trees if name != "concat"]
    if not segment_tree_names:
        raise ValueError("At least one non-concat tree is required.")

    sample_rows, ambiguous_codes, code_to_full_map = build_sample_identity(trees, segment_tree_names)

    group_anchor_values = {
        "flu2024": parse_anchor_values(args.anchor_flu2024),
        "guayas2023": parse_anchor_values(args.anchor_guayas2023),
        "pichincha0694": parse_anchor_values(args.anchor_pichincha0694),
    }
    anchor_sets = {
        "flu2024": anchors_for_group("flu2024", group_anchor_values, sample_rows),
        "guayas2023": anchors_for_group("guayas2023", group_anchor_values, sample_rows),
        "pichincha0694": anchors_for_group("pichincha0694", group_anchor_values, sample_rows),
    }

    sample_records: List[Dict[str, str]] = []

    for sample_full, sample_code in sample_rows:
        present_trees: List[str] = []
        tree_group_map: Dict[str, str] = {}
        tree_support_map: Dict[str, str] = {}
        tree_margin_map: Dict[str, str] = {}
        tree_neighbors: Dict[str, str] = {}
        low_support: List[str] = []
        outside_groups: List[str] = []

        baseline_votes: Counter = Counter()
        concat_group = ""
        concat_support = ""

        for tree_name, tree in trees.items():
            tip = find_tip(tree, sample_full=sample_full, sample_code=sample_code)
            if tip is None:
                continue

            present_trees.append(tree_name)

            sisters, local_support = sister_context(tree, tip)
            context_tips, context_support = clade_context(tree, tip)
            local_support_value = normalize_support(local_support)
            context_support_value = normalize_support(context_support)
            support_value = context_support_value if context_support_value is not None else local_support_value

            score_group = classify_group_score(
                sample_full=sample_full,
                sample_code=sample_code,
                context_tips=context_tips,
                group_anchor_values=group_anchor_values,
                code_to_full_map=code_to_full_map,
            )
            dist_group, dist_map, margin = group_by_anchor_distance(tree, tip.name, anchor_sets)

            if dist_group in GROUPS and (margin >= args.distance_margin or score_group == "outside_main_groups"):
                group = dist_group
            else:
                group = score_group

            tree_group_map[tree_name] = group
            tree_support_map[tree_name] = "" if support_value is None else f"{support_value:.2f}"
            tree_margin_map[tree_name] = f"{margin:.6f}"
            tree_neighbors[tree_name] = ";".join(sisters[:5])

            if tree_name != "concat" and group in GROUPS:
                baseline_votes[group] += 1

            if group == "outside_main_groups":
                outside_groups.append(tree_name)

            if support_value is None or support_value < args.min_support:
                low_support.append(tree_name)

            if tree_name == "concat":
                concat_group = group
                concat_support = tree_support_map[tree_name]

        baseline_group = baseline_votes.most_common(1)[0][0] if baseline_votes else "outside_main_groups"

        # Best non-baseline group by vote count
        non_baseline_votes = {g: v for g, v in baseline_votes.items() if g != baseline_group}
        if non_baseline_votes:
            secondary_group = max(non_baseline_votes, key=lambda g: non_baseline_votes[g])
            secondary_group_votes = non_baseline_votes[secondary_group]
        else:
            secondary_group = ""
            secondary_group_votes = 0

        inconsistent_trees: List[str] = []
        inconsistent_groups: List[Tuple[str, str]] = []
        inconsistent_trees_any_support: List[str] = []
        inconsistent_groups_any_support: List[Tuple[str, str]] = []

        for tree_name in trees:
            if tree_name == "concat":
                continue

            group = tree_group_map.get(tree_name, "missing")
            if group in ("missing", "outside_main_groups"):
                continue

            support_text = tree_support_map.get(tree_name, "")
            support_value = float(support_text) if support_text else None
            margin_text = tree_margin_map.get(tree_name, "")
            margin_value = float(margin_text) if margin_text else 0.0

            if group != baseline_group:
                inconsistent_trees_any_support.append(tree_name)
                inconsistent_groups_any_support.append((tree_name, group))

                if (support_value is not None and support_value >= args.min_support) or margin_value >= args.distance_margin:
                    inconsistent_trees.append(tree_name)
                    inconsistent_groups.append((tree_name, group))

        # alert_in_other_group is vote-based: requires the secondary group to have
        # accumulated >= secondary_min_votes tree placements (multi-tree signal).
        # Single-tree anomalies are still visible via inconsistent_trees columns.
        alert_in_other_group = "yes" if secondary_group_votes >= args.secondary_min_votes else "no"

        sample_records.append(
            {
                "sample_full_id": sample_full,
                "short_code": sample_code,
                "present_trees": ";".join(present_trees),
                "baseline_group": baseline_group,
                "baseline_votes": serialize_map(sorted(baseline_votes.items())),
                "secondary_group": secondary_group,
                "secondary_group_votes": str(secondary_group_votes) if secondary_group_votes else "",
                "alert_in_other_group": alert_in_other_group,
                "inconsistent_trees": ";".join(inconsistent_trees),
                "inconsistent_groups": serialize_map(inconsistent_groups),
                "inconsistent_trees_any_support": ";".join(inconsistent_trees_any_support),
                "inconsistent_groups_any_support": serialize_map(inconsistent_groups_any_support),
                "outside_main_groups_trees": ";".join(sorted(set(outside_groups))),
                "low_support_no_call_trees": ";".join(sorted(set(low_support))),
                "tree_group_map": serialize_map([(t, tree_group_map[t]) for t in present_trees if t in tree_group_map]),
                "tree_support_map": serialize_map([(t, tree_support_map[t]) for t in present_trees if t in tree_support_map]),
                "tree_distance_margin_map": serialize_map([(t, tree_margin_map[t]) for t in present_trees if t in tree_margin_map]),
                "representative_neighbors": serialize_map([(t, tree_neighbors[t]) for t in present_trees if t in tree_neighbors]),
                "concat_group": concat_group,
                "concat_support": concat_support,
                "ambiguity_flag": "yes" if sample_code in ambiguous_codes else "no",
            }
        )

    sample_records.sort(
        key=lambda row: (
            -(0 if not row["inconsistent_trees"] else len(row["inconsistent_trees"].split(";"))),
            row["sample_full_id"],
        )
    )

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fieldnames = [
        "sample_full_id",
        "short_code",
        "present_trees",
        "baseline_group",
        "baseline_votes",
        "secondary_group",
        "secondary_group_votes",
        "alert_in_other_group",
        "inconsistent_trees",
        "inconsistent_groups",
        "inconsistent_trees_any_support",
        "inconsistent_groups_any_support",
        "outside_main_groups_trees",
        "low_support_no_call_trees",
        "tree_group_map",
        "tree_support_map",
        "tree_distance_margin_map",
        "representative_neighbors",
        "concat_group",
        "concat_support",
        "ambiguity_flag",
    ]

    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_records)

    group_out_dir = os.path.dirname(args.group_candidates_out)
    if group_out_dir:
        os.makedirs(group_out_dir, exist_ok=True)

    group_rows: List[Dict[str, str]] = []
    for row in sample_records:
        votes_dict = parse_kv_str(row["baseline_votes"])
        base_group = row["baseline_group"]

        # Primary candidate — this is the sample's dominant group
        group_rows.append(
            {
                "baseline_group": base_group,
                "candidate_type": "primary",
                "votes_in_group": str(votes_dict.get(base_group, 0)),
                "sample_full_id": row["sample_full_id"],
                "alert_in_other_group": row["alert_in_other_group"],
                "secondary_group": row["secondary_group"],
                "secondary_group_votes": row["secondary_group_votes"],
                "inconsistent_trees": row["inconsistent_trees"],
                "inconsistent_trees_any_support": row["inconsistent_trees_any_support"],
                "outside_main_groups_trees": row["outside_main_groups_trees"],
                "present_trees": row["present_trees"],
            }
        )

        # Secondary candidates — groups where the sample has >= secondary_min_votes
        # tree placements despite not being the baseline. Even a single tree vote
        # qualifies (samples with only one present segment included).
        # alert_in_other_group uses the stricter secondary_min_votes threshold.
        for group, votes in votes_dict.items():
            if group != base_group and votes >= args.candidate_min_votes:
                group_rows.append(
                    {
                        "baseline_group": group,
                        "candidate_type": "secondary",
                        "votes_in_group": str(votes),
                        "sample_full_id": row["sample_full_id"],
                        "alert_in_other_group": row["alert_in_other_group"],
                        "secondary_group": row["secondary_group"],
                        "secondary_group_votes": row["secondary_group_votes"],
                        "inconsistent_trees": row["inconsistent_trees"],
                        "inconsistent_trees_any_support": row["inconsistent_trees_any_support"],
                        "outside_main_groups_trees": row["outside_main_groups_trees"],
                        "present_trees": row["present_trees"],
                    }
                )

    group_rows.sort(key=lambda r: (r["baseline_group"], r["candidate_type"], -int(r["votes_in_group"] or 0), r["sample_full_id"]))
    with open(args.group_candidates_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "baseline_group",
                "candidate_type",
                "votes_in_group",
                "sample_full_id",
                "alert_in_other_group",
                "secondary_group",
                "secondary_group_votes",
                "inconsistent_trees",
                "inconsistent_trees_any_support",
                "outside_main_groups_trees",
                "present_trees",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(group_rows)

    incong_out_dir = os.path.dirname(args.incongruences_only_out)
    if incong_out_dir:
        os.makedirs(incong_out_dir, exist_ok=True)

    incong_rows = [
        row
        for row in sample_records
        if row["alert_in_other_group"] == "yes" or row["outside_main_groups_trees"]
    ]

    with open(args.incongruences_only_out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(incong_rows)


if __name__ == "__main__":
    main()
