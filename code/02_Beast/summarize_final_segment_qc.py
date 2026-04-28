#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict


def read_tsv_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(value, default=0.0):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def segment_name_from_path(path):
    base = os.path.basename(path)
    stem = base.split(".")[0]
    return stem.replace("H5N1_", "")


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-segment QC for the final BEAST panel alignments"
    )
    parser.add_argument("--segment-metrics", nargs="+", required=True)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--out-outliers", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-report", required=True)
    args = parser.parse_args()

    per_taxon = defaultdict(
        lambda: {
            "role": "",
            "lineage": "",
            "segments_evaluated": 0,
            "segments_with_flags": 0,
            "segments_with_ungapped_change": 0,
            "total_outlier_score": 0.0,
            "max_segment_outlier_score": 0.0,
            "segments_flagged": [],
            "segments_ungapped_changed": [],
            "all_flags": set(),
            "max_after_n_fraction": 0.0,
            "max_after_ambiguous_fraction": 0.0,
            "max_after_consensus_distance_fraction": 0.0,
            "max_after_singleton_base_fraction": 0.0,
            "max_abs_delta_internal_gap_fraction": 0.0,
        }
    )

    total_rows = 0
    for metrics_path in args.segment_metrics:
        segment = segment_name_from_path(metrics_path)
        for row in read_tsv_rows(metrics_path):
            taxon = str(row.get("taxon", "")).strip()
            if not taxon:
                continue
            total_rows += 1
            entry = per_taxon[taxon]
            entry["role"] = entry["role"] or str(row.get("role", "")).strip()
            entry["lineage"] = entry["lineage"] or str(row.get("lineage", "")).strip()
            entry["segments_evaluated"] += 1

            outlier_score = as_float(row.get("outlier_score", 0.0))
            flag_count = int(as_float(row.get("flag_count", 0.0)))
            ungapped_changed = int(as_float(row.get("ungapped_sequence_changed", 0.0)))
            flags = [flag for flag in str(row.get("flags", "")).split(",") if flag]

            if flag_count > 0:
                entry["segments_with_flags"] += 1
                entry["segments_flagged"].append(segment)
            if ungapped_changed:
                entry["segments_with_ungapped_change"] += 1
                entry["segments_ungapped_changed"].append(segment)

            entry["total_outlier_score"] += outlier_score
            entry["max_segment_outlier_score"] = max(entry["max_segment_outlier_score"], outlier_score)
            entry["all_flags"].update(flags)
            entry["max_after_n_fraction"] = max(entry["max_after_n_fraction"], as_float(row.get("after_n_fraction", 0.0)))
            entry["max_after_ambiguous_fraction"] = max(
                entry["max_after_ambiguous_fraction"],
                as_float(row.get("after_ambiguous_fraction", 0.0)),
            )
            entry["max_after_consensus_distance_fraction"] = max(
                entry["max_after_consensus_distance_fraction"],
                as_float(row.get("after_consensus_distance_fraction", 0.0)),
            )
            entry["max_after_singleton_base_fraction"] = max(
                entry["max_after_singleton_base_fraction"],
                as_float(row.get("after_singleton_base_fraction", 0.0)),
            )
            entry["max_abs_delta_internal_gap_fraction"] = max(
                entry["max_abs_delta_internal_gap_fraction"],
                abs(as_float(row.get("delta_internal_gap_fraction", 0.0))),
            )

    metrics_rows = []
    for taxon, entry in per_taxon.items():
        metrics_rows.append(
            {
                "taxon": taxon,
                "role": entry["role"],
                "lineage": entry["lineage"],
                "segments_evaluated": entry["segments_evaluated"],
                "segments_with_flags": entry["segments_with_flags"],
                "segments_with_ungapped_change": entry["segments_with_ungapped_change"],
                "total_outlier_score": f"{entry['total_outlier_score']:.6f}",
                "max_segment_outlier_score": f"{entry['max_segment_outlier_score']:.6f}",
                "segments_flagged": ",".join(sorted(entry["segments_flagged"])),
                "segments_ungapped_changed": ",".join(sorted(entry["segments_ungapped_changed"])),
                "all_flags": ",".join(sorted(entry["all_flags"])),
                "max_after_n_fraction": f"{entry['max_after_n_fraction']:.6f}",
                "max_after_ambiguous_fraction": f"{entry['max_after_ambiguous_fraction']:.6f}",
                "max_after_consensus_distance_fraction": f"{entry['max_after_consensus_distance_fraction']:.6f}",
                "max_after_singleton_base_fraction": f"{entry['max_after_singleton_base_fraction']:.6f}",
                "max_abs_delta_internal_gap_fraction": f"{entry['max_abs_delta_internal_gap_fraction']:.6f}",
            }
        )

    metrics_rows.sort(
        key=lambda row: (
            -int(row["segments_with_flags"]),
            -float(row["total_outlier_score"]),
            -int(row["segments_with_ungapped_change"]),
            row["taxon"],
        )
    )

    metrics_fields = [
        "taxon",
        "role",
        "lineage",
        "segments_evaluated",
        "segments_with_flags",
        "segments_with_ungapped_change",
        "total_outlier_score",
        "max_segment_outlier_score",
        "segments_flagged",
        "segments_ungapped_changed",
        "all_flags",
        "max_after_n_fraction",
        "max_after_ambiguous_fraction",
        "max_after_consensus_distance_fraction",
        "max_after_singleton_base_fraction",
        "max_abs_delta_internal_gap_fraction",
    ]
    write_tsv(args.out_metrics, metrics_rows, metrics_fields)

    outlier_rows = [
        row for row in metrics_rows
        if int(row["segments_with_flags"]) > 0 or int(row["segments_with_ungapped_change"]) > 0
    ]
    write_tsv(
        args.out_outliers,
        outlier_rows,
        metrics_fields,
    )

    summary_rows = [
        {"metric": "n_taxa_total", "value": len(metrics_rows)},
        {"metric": "n_segment_rows_total", "value": total_rows},
        {"metric": "n_taxa_with_any_segment_flag", "value": sum(int(row["segments_with_flags"]) > 0 for row in metrics_rows)},
        {"metric": "n_taxa_with_multiple_segments_flagged", "value": sum(int(row["segments_with_flags"]) >= 2 for row in metrics_rows)},
        {"metric": "n_taxa_with_any_ungapped_change", "value": sum(int(row["segments_with_ungapped_change"]) > 0 for row in metrics_rows)},
        {"metric": "n_taxa_with_multiple_ungapped_changes", "value": sum(int(row["segments_with_ungapped_change"]) >= 2 for row in metrics_rows)},
    ]
    write_tsv(args.out_summary, summary_rows, ["metric", "value"])

    report_dir = os.path.dirname(args.out_report)
    if report_dir:
        os.makedirs(report_dir, exist_ok=True)
    with open(args.out_report, "w", encoding="utf-8") as handle:
        handle.write("# Final Segment QC Summary for the BEAST Panel\n\n")
        handle.write("## Inputs\n\n")
        for path in args.segment_metrics:
            handle.write(f"- `{path}`\n")
        handle.write("\n## Overview\n\n")
        for item in summary_rows:
            handle.write(f"- {item['metric']}: {item['value']}\n")
        handle.write("\n## Top Taxa With Segment-Level QC Signals\n\n")
        if not outlier_rows:
            handle.write("No se detectaron señales de QC en los alineamientos finales por segmento.\n")
        else:
            handle.write("| taxon | flagged_segments | ungapped_changed | total_score | flags |\n")
            handle.write("| --- | ---: | ---: | ---: | --- |\n")
            for row in outlier_rows[:20]:
                handle.write(
                    f"| {row['taxon']} | {row['segments_with_flags']} | {row['segments_with_ungapped_change']} | "
                    f"{float(row['total_outlier_score']):.2f} | {row['all_flags'] or '-'} |\n"
                )


if __name__ == "__main__":
    main()
