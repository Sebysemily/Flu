#!/usr/bin/env python3
import argparse
import csv
import math
import os
import statistics
from collections import Counter, defaultdict


SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
VALID_BASES = set("ACGT")
MISSING_BASES = set("N")


def read_fasta_records(path):
    records = []
    header = None
    chunks = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks).upper()))

    if not records:
        raise ValueError(f"No se encontraron secuencias en {path}")
    return records


def read_alignment(path):
    records = {}
    for header, seq in read_fasta_records(path):
        records[header] = seq
    return records


def parse_segment_header(header):
    parts = [part.strip() for part in header.split("/")]
    if len(parts) != 4:
        return None
    sample, segment, place, date_value = parts
    if not sample or not segment or not place or not date_value:
        return None
    return {
        "sample": sample,
        "segment": segment.upper(),
        "place": place,
        "date": date_value,
        "taxon": f"{sample}/{place}/{date_value}",
    }


def load_panel_annotations(path):
    annotations = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxon = str(row.get("taxon", "")).strip()
            if not taxon:
                continue
            annotations[taxon] = {
                "role": str(row.get("role", "")).strip(),
                "lineage": str(row.get("lineage", "")).strip(),
                "distance_to_seed": str(row.get("distance_to_seed", "")).strip(),
            }
    if not annotations:
        raise ValueError(f"No se encontraron taxa en {path}")
    return annotations


def basic_seq_metrics(seq):
    ungapped = seq.replace("-", "")
    letters = [base for base in ungapped if base.isalpha()]
    acgt = sum(base in VALID_BASES for base in letters)
    missing = sum(base in MISSING_BASES for base in letters)
    ambiguous = sum(base.isalpha() and base not in VALID_BASES and base not in MISSING_BASES for base in letters)
    informative = max(1, acgt + missing + ambiguous)
    return {
        "ungapped_length": len(ungapped),
        "acgt_count": acgt,
        "n_count": missing,
        "ambiguous_count": ambiguous,
        "n_fraction": missing / informative,
        "ambiguous_fraction": ambiguous / informative,
        "ungapped_sequence": ungapped,
    }


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


def basic_alignment_metrics(seq):
    ungapped = seq.replace("-", "")
    letters = [base for base in ungapped if base.isalpha()]
    acgt = sum(base in VALID_BASES for base in letters)
    missing = sum(base in MISSING_BASES for base in letters)
    ambiguous = sum(base.isalpha() and base not in VALID_BASES and base not in MISSING_BASES for base in letters)
    leading_gaps = len(seq) - len(seq.lstrip("-"))
    trailing_gaps = len(seq) - len(seq.rstrip("-"))
    total_gaps = seq.count("-")
    internal_gaps = total_gaps - leading_gaps - trailing_gaps
    informative = max(1, acgt + missing + ambiguous)
    internal_seq = seq[leading_gaps: len(seq) - trailing_gaps if trailing_gaps else len(seq)]
    return {
        "alignment_length": len(seq),
        "ungapped_sequence": ungapped,
        "ungapped_length": len(ungapped),
        "n_fraction": missing / informative,
        "ambiguous_fraction": ambiguous / informative,
        "gap_fraction": total_gaps / max(1, len(seq)),
        "internal_gap_fraction": internal_gaps / max(1, len(seq)),
        "gap_run_count": run_count(seq, "-"),
        "longest_gap_run": longest_run(seq, "-"),
        "leading_gap_count": leading_gaps,
        "trailing_gap_count": trailing_gaps,
        "internal_gap_count": internal_gaps,
    }


def build_alignment_profile(records):
    taxa = list(records.keys())
    seqs = list(records.values())
    alignment_length = len(seqs[0])
    if any(len(seq) != alignment_length for seq in seqs):
        raise ValueError("El archivo no es un alignment rectangular")

    same_base_counts = {taxon: 0 for taxon in taxa}
    rare_occupancy_counts = {taxon: 0 for taxon in taxa}
    informative_counts = {taxon: 0 for taxon in taxa}
    singleton_base_counts = {taxon: 0 for taxon in taxa}

    for pos in range(alignment_length):
        column = [records[taxon][pos] for taxon in taxa]
        acgt_counts = Counter(base for base in column if base in VALID_BASES)
        consensus_base = acgt_counts.most_common(1)[0][0] if acgt_counts else None
        occupancy = sum(base != "-" for base in column)

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
    median = statistics.median(values)
    deviations = [abs(v - median) for v in values]
    mad = statistics.median(deviations)
    return median, mad


def robust_positive_scores(values):
    median, mad = robust_scale(values)
    if mad > 0:
        scale = 1.4826 * mad
        return [(value - median) / scale if value > median else 0.0 for value in values]

    if len(values) >= 4:
        q1 = statistics.quantiles(values, n=4, method="inclusive")[0]
        q3 = statistics.quantiles(values, n=4, method="inclusive")[2]
        iqr = q3 - q1
        if iqr > 0:
            upper = q3 + 1.5 * iqr
            return [1.0 if value > upper else 0.0 for value in values]
    return [0.0 for _ in values]


def write_tsv(path, rows, fieldnames):
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_value(value):
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6f}"
    return value


def load_ecuador_summary(path):
    data = {}
    if not path:
        return data
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample = str(row.get("sample", "")).strip()
            if not sample:
                continue
            data.setdefault(sample, {"segments": [], "headers": []})
            data[sample]["segments"].append(str(row.get("segment", "")).strip().upper())
            data[sample]["headers"].append(str(row.get("header", "")).strip())
    return data


def load_context_summary(path):
    data = {}
    if not path:
        return data
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            header = str(row.get("header", "")).strip()
            parsed = parse_segment_header(header)
            if not parsed:
                continue
            entry = data.setdefault(parsed["taxon"], {
                "primary_accession": str(row.get("primary_accession", "")).strip(),
                "selection_role": str(row.get("selection_role", "")).strip(),
                "expected_segment_count": 0,
                "downloaded_segment_count": 0,
                "count_match": "",
            })
            expected = str(row.get("expected_segment_count", "")).strip()
            downloaded = str(row.get("downloaded_segment_count", "")).strip()
            if expected.isdigit():
                entry["expected_segment_count"] = max(entry["expected_segment_count"], int(expected))
            if downloaded.isdigit():
                entry["downloaded_segment_count"] = max(entry["downloaded_segment_count"], int(downloaded))
            match = str(row.get("count_match", "")).strip()
            if match:
                entry["count_match"] = match
    return data


def load_ecuador_audit(path):
    data = defaultdict(lambda: {"assembled": 0, "filled_with_N": 0, "other_status": 0})
    if not path:
        return data
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample = str(row.get("Código USFQ", "")).strip()
            status = str(row.get("status", "")).strip()
            if not sample:
                continue
            if status == "assembled":
                data[sample]["assembled"] += 1
            elif status == "filled_with_N":
                data[sample]["filled_with_N"] += 1
            else:
                data[sample]["other_status"] += 1
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Observational QC for the BEAST subset starting from source FASTA records"
    )
    parser.add_argument("--final-fasta", required=True)
    parser.add_argument("--panel-taxa", required=True)
    parser.add_argument("--ecuador-summary", default=None)
    parser.add_argument("--context-summary", default=None)
    parser.add_argument("--ecuador-audit", default=None)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--out-outliers", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--z-threshold", type=float, default=3.5)
    args = parser.parse_args()

    panel = load_panel_annotations(args.panel_taxa)
    subset_taxa = set(panel)

    final_records = read_fasta_records(args.final_fasta)
    final_by_taxon = defaultdict(dict)
    segment_lengths = defaultdict(list)
    for header, seq in final_records:
        parsed = parse_segment_header(header)
        if not parsed:
            continue
        final_by_taxon[parsed["taxon"]][parsed["segment"]] = seq
        segment_lengths[parsed["segment"]].append(basic_seq_metrics(seq)["ungapped_length"])

    missing_taxa = sorted(taxon for taxon in subset_taxa if taxon not in final_by_taxon)
    if missing_taxa:
        raise ValueError(f"Taxa del panel no encontrados en final FASTA: {missing_taxa[:5]}")

    segment_medians = {
        segment: statistics.median(lengths)
        for segment, lengths in segment_lengths.items()
        if lengths
    }

    ecuador_summary = load_ecuador_summary(args.ecuador_summary)
    context_summary = load_context_summary(args.context_summary)
    ecuador_audit = load_ecuador_audit(args.ecuador_audit)

    rows = []
    for taxon in sorted(subset_taxa):
        role = panel[taxon]["role"]
        lineage = panel[taxon]["lineage"]
        origin = "ecuador" if taxon.startswith("Flu-") else "context"
        source_sample = taxon.split("/", 1)[0]

        per_segment = final_by_taxon[taxon]
        present_segments = sorted(per_segment)
        missing_segments = [segment for segment in SEGMENTS if segment not in per_segment]

        total_acgt = 0
        total_n = 0
        total_ambig = 0
        length_ratios = []
        short_segments = []
        high_n_segments = []
        high_ambig_segments = []

        for segment in SEGMENTS:
            seq = per_segment.get(segment)
            if seq is None:
                continue
            metrics = basic_seq_metrics(seq)
            total_acgt += metrics["acgt_count"]
            total_n += metrics["n_count"]
            total_ambig += metrics["ambiguous_count"]

            segment_median = segment_medians.get(segment)
            if segment_median:
                ratio = metrics["ungapped_length"] / segment_median
                length_ratios.append(ratio)
                if ratio < 0.97:
                    short_segments.append(segment)
            if metrics["n_fraction"] >= 0.01:
                high_n_segments.append(segment)
            if metrics["ambiguous_fraction"] >= 0.005:
                high_ambig_segments.append(segment)

        total_letters = max(1, total_acgt + total_n + total_ambig)
        row = {
            "taxon": taxon,
            "origin": origin,
            "role": role,
            "lineage": lineage,
            "segments_present": len(present_segments),
            "missing_segments": ",".join(missing_segments),
            "raw_total_n_fraction": total_n / total_letters,
            "raw_total_ambiguous_fraction": total_ambig / total_letters,
            "raw_min_segment_length_ratio": min(length_ratios) if length_ratios else 0.0,
            "raw_median_segment_length_ratio": statistics.median(length_ratios) if length_ratios else 0.0,
            "raw_short_segment_count": len(short_segments),
            "raw_short_segments": ",".join(short_segments),
            "raw_high_n_segment_count": len(high_n_segments),
            "raw_high_n_segments": ",".join(high_n_segments),
            "raw_high_ambiguous_segment_count": len(high_ambig_segments),
            "raw_high_ambiguous_segments": ",".join(high_ambig_segments),
        }

        if origin == "ecuador":
            audit = ecuador_audit.get(source_sample, {"assembled": 0, "filled_with_N": 0, "other_status": 0})
            row["ecuador_mira_assembled_segments"] = audit["assembled"]
            row["ecuador_mira_filled_with_N_segments"] = audit["filled_with_N"]
            row["ecuador_mira_other_status_segments"] = audit["other_status"]
            row["context_expected_segment_count"] = ""
            row["context_downloaded_segment_count"] = ""
            row["context_count_match"] = ""
        else:
            ctx = context_summary.get(taxon, {})
            row["ecuador_mira_assembled_segments"] = ""
            row["ecuador_mira_filled_with_N_segments"] = ""
            row["ecuador_mira_other_status_segments"] = ""
            row["context_expected_segment_count"] = ctx.get("expected_segment_count", "")
            row["context_downloaded_segment_count"] = ctx.get("downloaded_segment_count", "")
            row["context_count_match"] = ctx.get("count_match", "")

        rows.append(row)

    metrics_to_flag = [
        ("raw_total_n_fraction", "high_raw_n_fraction", 0.01),
        ("raw_total_ambiguous_fraction", "high_raw_ambiguous_fraction", 0.005),
        ("raw_short_segment_count", "multiple_short_segments", 1.0),
    ]

    for metric, flag_name, min_value in metrics_to_flag:
        values = [row[metric] for row in rows]
        scores = robust_positive_scores(values)
        for row, score in zip(rows, scores):
            row[flag_name] = 1 if (score >= args.z_threshold or score == 1.0) and row[metric] >= min_value else 0
            row[f"{flag_name}_score"] = score

    min_ratio_values = [max(0.0, 1.0 - row["raw_min_segment_length_ratio"]) for row in rows]
    min_ratio_scores = robust_positive_scores(min_ratio_values)
    for row, score, deficit in zip(rows, min_ratio_scores, min_ratio_values):
        row["low_min_segment_length_ratio"] = 1 if (score >= args.z_threshold or score == 1.0) and deficit >= 0.03 else 0
        row["low_min_segment_length_ratio_score"] = score

        row["subset_missing_segments_flag"] = 1 if row["segments_present"] < 8 else 0
        row["context_segment_count_mismatch"] = 1 if str(row["context_count_match"]).lower() == "false" else 0
        row["ecuador_mira_incomplete_history"] = 1 if (
            row["origin"] == "ecuador" and str(row["ecuador_mira_filled_with_N_segments"]).isdigit()
            and int(row["ecuador_mira_filled_with_N_segments"]) > 0
        ) else 0

        flags = []
        for key in [
            "subset_missing_segments_flag",
            "high_raw_n_fraction",
            "high_raw_ambiguous_fraction",
            "multiple_short_segments",
            "low_min_segment_length_ratio",
            "context_segment_count_mismatch",
            "ecuador_mira_incomplete_history",
        ]:
            if row.get(key):
                flags.append(key)

        row["flags"] = ",".join(flags)
        row["flag_count"] = len(flags)
        row["outlier_score"] = row["flag_count"]

    rows.sort(
        key=lambda row: (
            -row["outlier_score"],
            -row["flag_count"],
            -row["raw_total_n_fraction"],
            -row["raw_short_segment_count"],
            row["taxon"],
        )
    )

    fieldnames = [
        "taxon",
        "origin",
        "role",
        "lineage",
        "outlier_score",
        "flag_count",
        "flags",
        "segments_present",
        "missing_segments",
        "raw_total_n_fraction",
        "raw_total_ambiguous_fraction",
        "raw_min_segment_length_ratio",
        "raw_median_segment_length_ratio",
        "raw_short_segment_count",
        "raw_short_segments",
        "raw_high_n_segment_count",
        "raw_high_n_segments",
        "raw_high_ambiguous_segment_count",
        "raw_high_ambiguous_segments",
        "context_expected_segment_count",
        "context_downloaded_segment_count",
        "context_count_match",
        "ecuador_mira_assembled_segments",
        "ecuador_mira_filled_with_N_segments",
        "ecuador_mira_other_status_segments",
        "subset_missing_segments_flag",
        "high_raw_n_fraction",
        "high_raw_n_fraction_score",
        "high_raw_ambiguous_fraction",
        "high_raw_ambiguous_fraction_score",
        "multiple_short_segments",
        "multiple_short_segments_score",
        "low_min_segment_length_ratio",
        "low_min_segment_length_ratio_score",
        "context_segment_count_mismatch",
        "ecuador_mira_incomplete_history",
    ]

    formatted_rows = []
    for row in rows:
        formatted = {field: format_value(row[field]) for field in fieldnames}
        formatted_rows.append(formatted)
    write_tsv(args.out_metrics, formatted_rows, fieldnames)

    outliers = [row for row in rows if row["outlier_score"] > 0]
    outlier_fields = [
        "taxon",
        "origin",
        "role",
        "lineage",
        "outlier_score",
        "flag_count",
        "flags",
        "raw_total_n_fraction",
        "raw_total_ambiguous_fraction",
        "raw_min_segment_length_ratio",
        "raw_short_segment_count",
    ]
    formatted_outliers = []
    for row in outliers:
        formatted_outliers.append({field: format_value(row[field]) for field in outlier_fields})
    write_tsv(args.out_outliers, formatted_outliers, outlier_fields)

    summary_rows = []
    summary_flags = [
        "subset_missing_segments_flag",
        "high_raw_n_fraction",
        "high_raw_ambiguous_fraction",
        "multiple_short_segments",
        "low_min_segment_length_ratio",
        "context_segment_count_mismatch",
        "ecuador_mira_incomplete_history",
    ]
    for flag in summary_flags:
        summary_rows.append({"metric": flag, "n_taxa_flagged": sum(row[flag] for row in rows)})
    summary_rows.extend(
        [
            {"metric": "any_outlier_flag", "n_taxa_flagged": len(outliers)},
            {"metric": "n_taxa_total", "n_taxa_flagged": len(rows)},
        ]
    )
    write_tsv(args.out_summary, summary_rows, ["metric", "n_taxa_flagged"])

    report_dir = os.path.dirname(args.out_report)
    if report_dir:
        os.makedirs(report_dir, exist_ok=True)
    with open(args.out_report, "w", encoding="utf-8") as handle:
        handle.write("# Source-first QC for the BEAST subset\n\n")
        handle.write("## Inputs\n\n")
        handle.write(f"- Final FASTA: `{args.final_fasta}`\n")
        handle.write(f"- Panel taxa: `{args.panel_taxa}`\n")
        if args.ecuador_summary:
            handle.write(f"- Ecuador summary: `{args.ecuador_summary}`\n")
        if args.context_summary:
            handle.write(f"- Context summary: `{args.context_summary}`\n")
        if args.ecuador_audit:
            handle.write(f"- Ecuador audit: `{args.ecuador_audit}`\n")

        handle.write("\n## Overview\n\n")
        handle.write(f"- Taxa in BEAST panel: {len(rows)}\n")
        handle.write(f"- Candidate outliers with >=1 flag: {len(outliers)}\n")
        handle.write(f"- Ecuador taxa in panel: {sum(row['origin'] == 'ecuador' for row in rows)}\n")
        handle.write(f"- Context taxa in panel: {sum(row['origin'] == 'context' for row in rows)}\n")

        handle.write("\n## Flag Counts\n\n")
        for item in summary_rows:
            handle.write(f"- {item['metric']}: {item['n_taxa_flagged']}\n")

        handle.write("\n## Top Candidate Outliers\n\n")
        if not outliers:
            handle.write("No se detectaron outliers con las reglas actuales.\n")
        else:
            handle.write("| taxon | origin | role | score | flags | raw_N | short_segments | min_len_ratio |\n")
            handle.write("| --- | --- | --- | ---: | --- | ---: | ---: | ---: |\n")
            for row in outliers[:20]:
                handle.write(
                    f"| {row['taxon']} | {row['origin']} | {row['role'] or '-'} | {row['outlier_score']} | "
                    f"{row['flags'] or '-'} | {row['raw_total_n_fraction']:.4f} | "
                    f"{row['raw_short_segment_count']} | "
                    f"{row['raw_min_segment_length_ratio']:.4f} |\n"
                )


if __name__ == "__main__":
    main()
