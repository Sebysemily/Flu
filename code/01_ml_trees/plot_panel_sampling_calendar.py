#!/usr/bin/env python3
"""Summarize HA panel sampling density by month and role."""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", required=True)
    p.add_argument("--context-taxa", required=True)
    p.add_argument("--postqc-fasta", required=True)
    p.add_argument("--costa-audit", default=None)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-png", required=True)
    return p.parse_args()


def load_ids(path: str) -> set[str]:
    ids: set[str] = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def load_fasta_ids(path: str) -> set[str]:
    ids: set[str] = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip())
    return ids


def year_month(date_str: str) -> str:
    date_str = (date_str or "").strip()
    if not date_str:
        return "unknown"
    parts = date_str.split("-")
    if len(parts) >= 2:
        return f"{parts[0]}-{parts[1]}"
    return parts[0]


def load_audit_window(path: str | None) -> tuple[str, str]:
    if not path or not Path(path).exists():
        return "", ""
    with open(path, encoding="utf-8") as handle:
        rows = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    return rows.get("window_start", ""), rows.get("window_end", "")


def month_in_window(ym: str, start: str, end: str) -> bool:
    if not start or not end or ym == "unknown":
        return False
    probe = f"{ym}-15"
    return start <= probe <= end


def main() -> None:
    args = parse_args()
    context_ids = load_ids(args.context_taxa)
    postqc_ids = load_fasta_ids(args.postqc_fasta)
    window_start, window_end = load_audit_window(args.costa_audit)

    rows = list(csv.DictReader(open(args.metadata, encoding="utf-8")))

    buckets = ["flu_costa", "flu_sierra", "regional_context", "american_anchor", "other"]
    stages = ["context_panel", "postqc_panel"]
    counts: dict[tuple[str, str, str], int] = defaultdict(int)

    for row in rows:
        name = (row.get("file_name") or "").strip()
        if not name:
            continue
        role = (row.get("expected_role") or "").strip()
        if role.startswith("flu_costa"):
            bucket = "flu_costa"
        elif role.startswith("flu_sierra") or role.startswith("flu_amazonia"):
            bucket = "flu_sierra"
        elif role == "regional_context":
            bucket = "regional_context"
        elif role == "american_anchor":
            bucket = "american_anchor"
        else:
            bucket = "other"

        ym = year_month(row.get("collection_date", ""))
        if name in context_ids:
            counts[(ym, bucket, "context_panel")] += 1
        if name in postqc_ids:
            counts[(ym, bucket, "postqc_panel")] += 1

    months = sorted({m for m, _, _ in counts.keys() if m != "unknown"})
    out_dir = Path(args.out_tsv).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.out_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["year_month", "role_group", "stage", "count"])
        for ym in months:
            for bucket in buckets:
                for stage in stages:
                    writer.writerow([ym, bucket, stage, counts.get((ym, bucket, stage), 0)])

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    x = list(range(len(months)))
    width = 0.18
    flu_roles = ["flu_costa", "flu_sierra"]
    for ax_idx, (stage, title) in enumerate(
        [("context_panel", "Panel tras augur (context_taxa)"), ("postqc_panel", "Panel post-QC (postQC fasta)")]
    ):
        ax = axes[ax_idx]
        for i, role in enumerate(flu_roles):
            vals = [counts.get((m, role, stage), 0) for m in months]
            ax.bar([xi + (i - 0.5) * width for xi in x], vals, width=width, label=role)
        rc_vals = [counts.get((m, "regional_context", stage), 0) for m in months]
        ax.bar([xi + width for xi in x], rc_vals, width=width, label="regional_context", alpha=0.85)
        if window_start and window_end:
            for xi, ym in enumerate(months):
                if month_in_window(ym, window_start, window_end):
                    ax.axvspan(xi - 0.5, xi + 0.5, color="gold", alpha=0.12)
            ax.plot([], [], color="gold", alpha=0.4, label=f"flu_costa window {window_start}..{window_end}")
        ax.set_title(title)
        ax.set_ylabel("Taxa")
        ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(months, rotation=45, ha="right")
    axes[-1].set_xlabel("year_month")
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {args.out_tsv} and {args.out_png}")


if __name__ == "__main__":
    main()
