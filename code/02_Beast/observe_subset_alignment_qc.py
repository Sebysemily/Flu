#!/usr/bin/env python3
import argparse
import csv
import math
import os
import statistics
from collections import Counter


VALID_BASES = set("ACGT")
MISSING_BASES = set("N")


def read_fasta(path):
    records = {}
    header = None
    chunks = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks).upper()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        records[header] = "".join(chunks).upper()

    if not records:
        raise ValueError(f"No se encontraron secuencias en {path}")
    return records


def load_taxon_annotations(path):
    annotations = {}
    if not path:
        return annotations

    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxon = str(row.get("taxon", "")).strip()
            if not taxon:
                continue
            annotations[taxon] = {
                "role": str(row.get("role", "")).strip(),
                "lineage": str(row.get("lineage", "")).strip(),
            }
    return annotations


def longest_run(seq, char):
    best = 0
    current = 0
    for base in seq:
        if base == char:
            current += 1
            if current > best:
                best = current
        else:
            current = 0
    return best


def run_count(seq, char):
    runs = 0
    in_run = False
    for base in seq:
        if base == char:
            if not in_run:
                runs += 1
                in_run = True
        else:
            in_run = False
    return runs


def basic_metrics(seq):
    ungapped = seq.replace("-", "")
    letters = [base for base in ungapped if base.isalpha()]
    acgt = sum(base in VALID_BASES for base in letters)
    missing = sum(base in MISSING_BASES for base in letters)
    ambiguous = sum(base.isalpha() and base not in VALID_BASES and base not in MISSING_BASES for base in letters)

    leading_gaps = len(seq) - len(seq.lstrip("-"))
    trailing_gaps = len(seq) - len(seq.rstrip("-"))
    total_gaps = seq.count("-")
    internal_gaps = total_gaps - leading_gaps - trailing_gaps
    internal_seq = seq[leading_gaps: len(seq) - trailing_gaps if trailing_gaps else len(seq)]

    informative = max(1, acgt + missing + ambiguous)
    return {
        "alignment_length": len(seq),
        "ungapped_length": len(ungapped),
        "acgt_count": acgt,
        "n_count": missing,
        "ambiguous_count": ambiguous,
        "n_fraction": missing / informative,
        "ambiguous_fraction": ambiguous / informative,
        "gap_count": total_gaps,
        "gap_fraction": total_gaps / max(1, len(seq)),
        "leading_gap_count": leading_gaps,
        "trailing_gap_count": trailing_gaps,
        "internal_gap_count": internal_gaps,
        "internal_gap_fraction": internal_gaps / max(1, len(seq)),
        "gap_run_count": run_count(seq, "-"),
        "internal_gap_run_count": run_count(internal_seq, "-"),
        "longest_gap_run": longest_run(seq, "-"),
        "ungapped_sequence": ungapped,
    }


def build_alignment_profile(records):
    taxa = list(records.keys())
    seqs = list(records.values())
    alignment_length = len(seqs[0])
    if any(len(seq) != alignment_length for seq in seqs):
        raise ValueError("El archivo no es un alignment rectangular: no todas las secuencias tienen la misma longitud")

    consensus = []
    non_gap_occupancy = []
    same_base_counts = {taxon: 0 for taxon in taxa}
    rare_occupancy_counts = {taxon: 0 for taxon in taxa}
    informative_counts = {taxon: 0 for taxon in taxa}
    singleton_base_counts = {taxon: 0 for taxon in taxa}

    for pos in range(alignment_length):
        column = [records[taxon][pos] for taxon in taxa]
        acgt_counts = Counter(base for base in column if base in VALID_BASES)
        consensus_base = acgt_counts.most_common(1)[0][0] if acgt_counts else None
        consensus.append(consensus_base)

        occupancy = sum(base != "-" for base in column)
        non_gap_occupancy.append(occupancy)

        for taxon in taxa:
            base = records[taxon][pos]
            if base not in VALID_BASES:
                continue

            informative_counts[taxon] += 1
            if consensus_base == base:
                same_base_counts[taxon] += 1
            if occupancy <= 2:
                rare_occupancy_counts[taxon] += 1
            if acgt_counts.get(base, 0) == 1:
                singleton_base_counts[taxon] += 1

    profile = {}
    for taxon in taxa:
        informative = informative_counts[taxon]
        matching = same_base_counts[taxon]
        profile[taxon] = {
            "consensus_distance_fraction": (
                (informative - matching) / informative if informative else 0.0
            ),
            "rare_occupancy_fraction": (
                rare_occupancy_counts[taxon] / informative if informative else 0.0
            ),
            "singleton_base_fraction": (
                singleton_base_counts[taxon] / informative if informative else 0.0
            ),
        }
    return profile


def robust_scale(values):
    if not values:
        return 0.0, 0.0
    median = statistics.median(values)
    deviations = [abs(v - median) for v in values]
    mad = statistics.median(deviations)
    return median, mad


def robust_positive_outliers(values, threshold):
    median, mad = robust_scale(values)
    if mad > 0:
        scale = 1.4826 * mad
        return [
            (value - median) / scale if value > median else 0.0
            for value in values
        ]

    q1 = statistics.quantiles(values, n=4, method="inclusive")[0]
    q3 = statistics.quantiles(values, n=4, method="inclusive")[2]
    iqr = q3 - q1
    if iqr > 0:
        upper = q3 + 1.5 * iqr
        return [1.0 if value > upper else 0.0 for value in values]
    return [0.0 for _ in values]


def metric_median(rows, key):
    return statistics.median(row[key] for row in rows)


def write_tsv(path, rows, fieldnames):
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_float(value):
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6f}"
    return value


def main():
    parser = argparse.ArgumentParser(
        description="Observational QC for raw subset FASTA vs final aligned subset used in BEAST"
    )
    parser.add_argument("--before-alignment", required=True)
    parser.add_argument("--after-alignment", required=True)
    parser.add_argument("--taxa-tsv", default=None)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--out-outliers", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--z-threshold", type=float, default=3.5)
    args = parser.parse_args()

    before_records = read_fasta(args.before_alignment)
    after_records = read_fasta(args.after_alignment)
    before_taxa = set(before_records)
    after_taxa = set(after_records)
    if before_taxa != after_taxa:
        missing_before = sorted(after_taxa - before_taxa)
        missing_after = sorted(before_taxa - after_taxa)
        raise ValueError(
            "Los alignments no tienen el mismo set de taxa. "
            f"Solo en after: {missing_before[:5]}; solo en before: {missing_after[:5]}"
        )

    annotations = load_taxon_annotations(args.taxa_tsv)
    after_profile = build_alignment_profile(after_records)

    rows = []
    taxa = sorted(before_records)
    for taxon in taxa:
        before_basic = basic_metrics(before_records[taxon])
        after_basic = basic_metrics(after_records[taxon])
        role = annotations.get(taxon, {}).get("role", "")
        lineage = annotations.get(taxon, {}).get("lineage", "")

        row = {
            "taxon": taxon,
            "role": role,
            "lineage": lineage,
            "ungapped_sequence_changed": int(
                before_basic["ungapped_sequence"] != after_basic["ungapped_sequence"]
            ),
            "before_alignment_length": before_basic["alignment_length"],
            "after_alignment_length": after_basic["alignment_length"],
            "before_ungapped_length": before_basic["ungapped_length"],
            "after_ungapped_length": after_basic["ungapped_length"],
            "before_n_fraction": before_basic["n_fraction"],
            "after_n_fraction": after_basic["n_fraction"],
            "before_ambiguous_fraction": before_basic["ambiguous_fraction"],
            "after_ambiguous_fraction": after_basic["ambiguous_fraction"],
            "before_gap_fraction": before_basic["gap_fraction"],
            "after_gap_fraction": after_basic["gap_fraction"],
            "before_internal_gap_fraction": before_basic["internal_gap_fraction"],
            "after_internal_gap_fraction": after_basic["internal_gap_fraction"],
            "before_leading_gap_count": before_basic["leading_gap_count"],
            "after_leading_gap_count": after_basic["leading_gap_count"],
            "before_trailing_gap_count": before_basic["trailing_gap_count"],
            "after_trailing_gap_count": after_basic["trailing_gap_count"],
            "before_gap_run_count": before_basic["gap_run_count"],
            "after_gap_run_count": after_basic["gap_run_count"],
            "before_longest_gap_run": before_basic["longest_gap_run"],
            "after_longest_gap_run": after_basic["longest_gap_run"],
            "after_consensus_distance_fraction": after_profile[taxon]["consensus_distance_fraction"],
            "after_rare_occupancy_fraction": after_profile[taxon]["rare_occupancy_fraction"],
            "after_singleton_base_fraction": after_profile[taxon]["singleton_base_fraction"],
        }

        row["delta_gap_fraction"] = row["after_gap_fraction"] - row["before_gap_fraction"]
        row["delta_internal_gap_fraction"] = (
            row["after_internal_gap_fraction"] - row["before_internal_gap_fraction"]
        )
        row["delta_n_fraction"] = row["after_n_fraction"] - row["before_n_fraction"]
        row["delta_ambiguous_fraction"] = row["after_ambiguous_fraction"] - row["before_ambiguous_fraction"]
        rows.append(row)

    median_internal_gap_shift = metric_median(rows, "delta_internal_gap_fraction")
    for row in rows:
        row["internal_gap_shift_deviation"] = abs(
            row["delta_internal_gap_fraction"] - median_internal_gap_shift
        )

    metric_specs = [
        ("after_n_fraction", "high_n_fraction", 0.01),
        ("after_ambiguous_fraction", "high_ambiguous_fraction", 0.005),
        ("after_consensus_distance_fraction", "high_consensus_distance_after", 0.01),
        ("after_rare_occupancy_fraction", "high_private_insertion_burden_after", 0.005),
        ("after_singleton_base_fraction", "high_singleton_substitution_burden_after", 0.003),
        ("internal_gap_shift_deviation", "atypical_internal_gap_shift", 0.01),
    ]

    for metric, flag_name, min_value in metric_specs:
        values = [row[metric] for row in rows]
        scores = robust_positive_outliers(values, args.z_threshold)
        for row, score in zip(rows, scores):
            if (score >= args.z_threshold or score == 1.0) and row[metric] >= min_value:
                row[flag_name] = 1
            else:
                row[flag_name] = 0
            row[f"{flag_name}_score"] = score

    median_ungapped_length = metric_median(rows, "after_ungapped_length")
    for row in rows:
        low_length_value = max(0.0, median_ungapped_length - row["after_ungapped_length"])
        row["_length_deficit"] = low_length_value
    length_scores = robust_positive_outliers([row["_length_deficit"] for row in rows], args.z_threshold)
    for row, score in zip(rows, length_scores):
        min_length_deficit = max(250.0, median_ungapped_length * 0.05)
        row["short_ungapped_length"] = (
            1 if (score >= args.z_threshold or score == 1.0) and row["_length_deficit"] >= min_length_deficit else 0
        )
        row["short_ungapped_length_score"] = score

    for row in rows:
        flags = []
        for key in [
            "high_n_fraction",
            "high_ambiguous_fraction",
            "short_ungapped_length",
            "high_consensus_distance_after",
            "high_private_insertion_burden_after",
            "high_singleton_substitution_burden_after",
            "atypical_internal_gap_shift",
        ]:
            if row.get(key):
                flags.append(key)
        if row["ungapped_sequence_changed"]:
            flags.append("ungapped_sequence_changed")
        row["flag_count"] = len(flags)
        row["flags"] = ",".join(flags)
        row["outlier_score"] = row["flag_count"] + (
            3 if row["ungapped_sequence_changed"] else 0
        )
        del row["_length_deficit"]

    rows.sort(
        key=lambda row: (
            -row["outlier_score"],
            -row["flag_count"],
            -row["after_consensus_distance_fraction"],
            row["taxon"],
        )
    )

    metrics_fieldnames = [
        "taxon",
        "role",
        "lineage",
        "outlier_score",
        "flag_count",
        "flags",
        "ungapped_sequence_changed",
        "before_alignment_length",
        "after_alignment_length",
        "before_ungapped_length",
        "after_ungapped_length",
        "before_n_fraction",
        "after_n_fraction",
        "before_ambiguous_fraction",
        "after_ambiguous_fraction",
        "before_gap_fraction",
        "after_gap_fraction",
        "delta_gap_fraction",
        "before_internal_gap_fraction",
        "after_internal_gap_fraction",
        "delta_internal_gap_fraction",
        "before_leading_gap_count",
        "after_leading_gap_count",
        "before_trailing_gap_count",
        "after_trailing_gap_count",
        "before_gap_run_count",
        "after_gap_run_count",
        "before_longest_gap_run",
        "after_longest_gap_run",
        "after_consensus_distance_fraction",
        "after_rare_occupancy_fraction",
        "after_singleton_base_fraction",
        "delta_n_fraction",
        "delta_ambiguous_fraction",
        "high_n_fraction",
        "high_n_fraction_score",
        "high_ambiguous_fraction",
        "high_ambiguous_fraction_score",
        "short_ungapped_length",
        "short_ungapped_length_score",
        "high_consensus_distance_after",
        "high_consensus_distance_after_score",
        "high_private_insertion_burden_after",
        "high_private_insertion_burden_after_score",
        "high_singleton_substitution_burden_after",
        "high_singleton_substitution_burden_after_score",
        "internal_gap_shift_deviation",
        "atypical_internal_gap_shift",
        "atypical_internal_gap_shift_score",
    ]

    formatted_rows = []
    for row in rows:
        formatted = {}
        for field in metrics_fieldnames:
            formatted[field] = format_float(row[field])
        formatted_rows.append(formatted)
    write_tsv(args.out_metrics, formatted_rows, metrics_fieldnames)

    outliers = [row for row in rows if row["outlier_score"] > 0]
    outlier_fieldnames = [
        "taxon",
        "role",
        "lineage",
        "outlier_score",
        "flag_count",
        "flags",
        "after_ungapped_length",
        "after_n_fraction",
        "after_ambiguous_fraction",
        "after_consensus_distance_fraction",
        "after_rare_occupancy_fraction",
        "after_singleton_base_fraction",
        "delta_internal_gap_fraction",
    ]
    formatted_outliers = []
    for row in outliers:
        formatted = {}
        for field in outlier_fieldnames:
            formatted[field] = format_float(row[field])
        formatted_outliers.append(formatted)
    write_tsv(args.out_outliers, formatted_outliers, outlier_fieldnames)

    summary_rows = []
    flag_names = [
        "high_n_fraction",
        "high_ambiguous_fraction",
        "short_ungapped_length",
        "high_consensus_distance_after",
        "high_private_insertion_burden_after",
        "high_singleton_substitution_burden_after",
        "atypical_internal_gap_shift",
    ]
    for flag_name in flag_names:
        summary_rows.append(
            {
                "metric": flag_name,
                "n_taxa_flagged": sum(row[flag_name] for row in rows),
            }
        )
    summary_rows.extend(
        [
            {"metric": "ungapped_sequence_changed", "n_taxa_flagged": sum(row["ungapped_sequence_changed"] for row in rows)},
            {"metric": "any_outlier_flag", "n_taxa_flagged": len(outliers)},
            {"metric": "n_taxa_total", "n_taxa_flagged": len(rows)},
            {"metric": "after_alignment_length", "n_taxa_flagged": rows[0]["after_alignment_length"]},
        ]
    )
    write_tsv(args.out_summary, summary_rows, ["metric", "n_taxa_flagged"])

    top_rows = outliers[:15]
    with open(args.out_report, "w", encoding="utf-8") as handle:
        handle.write("# Observational QC: raw subset vs final aligned subset\n\n")
        handle.write("## Inputs\n\n")
        handle.write(f"- Raw subset FASTA: `{args.before_alignment}`\n")
        handle.write(f"- Final aligned subset used in BEAST: `{args.after_alignment}`\n")
        if args.taxa_tsv:
            handle.write(f"- Taxa annotations: `{args.taxa_tsv}`\n")
        handle.write("\n## Overview\n\n")
        handle.write(f"- Taxa compared: {len(rows)}\n")
        handle.write(f"- Final alignment length: {rows[0]['after_alignment_length']}\n")
        handle.write(f"- Taxa with ungapped sequence changed: {sum(row['ungapped_sequence_changed'] for row in rows)}\n")
        handle.write(f"- Candidate outliers with >=1 flag: {len(outliers)}\n")
        handle.write("\n## Flag Counts\n\n")
        for item in summary_rows:
            handle.write(f"- {item['metric']}: {item['n_taxa_flagged']}\n")

        handle.write("\n## Top Candidate Outliers\n\n")
        if not top_rows:
            handle.write("No se detectaron outliers con las reglas actuales.\n")
        else:
            handle.write("| taxon | role | score | flags | after_consensus_distance | after_private_insertion | internal_gap_delta |\n")
            handle.write("| --- | --- | ---: | --- | ---: | ---: | ---: |\n")
            for row in top_rows:
                handle.write(
                    f"| {row['taxon']} | {row['role'] or '-'} | {row['outlier_score']} | "
                    f"{row['flags'] or '-'} | "
                    f"{row['after_consensus_distance_fraction']:.4f} | "
                    f"{row['after_rare_occupancy_fraction']:.4f} | "
                    f"{row['delta_internal_gap_fraction']:.4f} |\n"
                )


if __name__ == "__main__":
    main()
