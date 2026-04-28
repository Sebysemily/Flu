#!/usr/bin/env python3
import argparse
import csv
import os


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


def prune_dates_tsv(dates_in, dates_out, outlier_taxa):
    with open(dates_in, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or ["name", "date"]
        rows = [
            row
            for row in reader
            if str(row.get("name", "")).strip() not in outlier_taxa
        ]

    write_tsv(dates_out, rows, fieldnames)


def read_rtt_outliers(path):
    outliers = {}
    if not os.path.exists(path):
        return outliers

    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return outliers
        for row in reader:
            if not row:
                continue
            taxon = str(row[0]).strip()
            if not taxon:
                continue
            outliers[taxon] = {
                "given_date": row[1].strip() if len(row) > 1 else "",
                "apparent_date": row[2].strip() if len(row) > 2 else "",
                "residual": row[3].strip() if len(row) > 3 else "",
            }
    return outliers


def main():
    parser = argparse.ArgumentParser(
        description="Filter BEAST panel taxa using TreeTime root-to-tip outliers"
    )
    parser.add_argument("--panel-taxa", required=True)
    parser.add_argument("--rtt-outliers", required=True)
    parser.add_argument("--filtered-panel-out", required=True)
    parser.add_argument("--exclusions-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--dates-in", default=None)
    parser.add_argument("--dates-out", default=None)
    args = parser.parse_args()

    panel_rows = read_tsv_rows(args.panel_taxa)
    rtt_outliers = read_rtt_outliers(args.rtt_outliers)

    filtered_rows = []
    excluded_rows = []

    for row in panel_rows:
        taxon = str(row.get("taxon", "")).strip()
        outlier = rtt_outliers.get(taxon)
        if outlier:
            excluded_rows.append(
                {
                    "taxon": taxon,
                    "excluded": "true",
                    "reason_code": "treetime_rtt_outlier",
                    "reason_detail": (
                        "TreeTime ignored this tip during root-to-tip rate estimation"
                    ),
                    "given_date": outlier["given_date"],
                    "apparent_date": outlier["apparent_date"],
                    "residual": outlier["residual"],
                }
            )
        else:
            filtered_rows.append(row)

    write_tsv(
        args.filtered_panel_out,
        filtered_rows,
        ["taxon", "role", "lineage", "distance_to_seed"],
    )
    write_tsv(
        args.exclusions_out,
        excluded_rows,
        [
            "taxon",
            "excluded",
            "reason_code",
            "reason_detail",
            "given_date",
            "apparent_date",
            "residual",
        ],
    )
    write_tsv(
        args.summary_out,
        [
            {"metric": "panel_taxa_input", "value": len(panel_rows)},
            {"metric": "rtt_outliers_detected", "value": len(rtt_outliers)},
            {"metric": "panel_taxa_excluded", "value": len(excluded_rows)},
            {"metric": "panel_taxa_kept", "value": len(filtered_rows)},
        ],
        ["metric", "value"],
    )

    if args.dates_in and args.dates_out:
        prune_dates_tsv(args.dates_in, args.dates_out, set(rtt_outliers.keys()))


if __name__ == "__main__":
    main()
