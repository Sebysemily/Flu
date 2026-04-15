#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


MARGINAL_PATTERNS = [
    re.compile(
        r"marginal(?:\s+log)?\s+likelihood(?:\s+estimate)?\s*[:=]\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
        re.IGNORECASE,
    ),
    re.compile(
        r"log\s+marginal\s+likelihood(?:\s+estimate)?\s*[:=]\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
        re.IGNORECASE,
    ),
    re.compile(
        r"log\s+marginal\s+L\s+estimate\s*[:=]\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
        re.IGNORECASE,
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize BETS marginal likelihood outputs.")
    parser.add_argument("--out-summary", required=True, help="Per-scenario BETS summary TSV.")
    parser.add_argument(
        "--out-bayes-factors",
        required=True,
        help="Pairwise heterochronous vs isochronous summary TSV.",
    )
    parser.add_argument(
        "--scenario-log",
        action="append",
        default=[],
        help="Scenario to stdout mapping, e.g. strict_constant_heterochronous=path/to/stdout.log",
    )
    parser.add_argument(
        "--scenario-err",
        action="append",
        default=[],
        help="Scenario to stderr mapping, e.g. strict_constant_heterochronous=path/to/stderr.log",
    )
    return parser.parse_args()


def parse_mapping(items: list[str]) -> dict[str, Path]:
    mapping: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Expected scenario mapping in the form name=path, got: {item}")
        name, raw_path = item.split("=", 1)
        mapping[name] = Path(raw_path)
    return mapping


def extract_marginal_log_likelihood(text: str) -> float | None:
    matches: list[float] = []
    for pattern in MARGINAL_PATTERNS:
        matches.extend(float(match.group(1)) for match in pattern.finditer(text))
    if not matches:
        return None
    return matches[-1]


def evidence_label(delta: float) -> str:
    magnitude = abs(delta)
    if math.isclose(magnitude, 0.0, abs_tol=1e-12):
        return "none"
    if magnitude >= 10:
        return "very_strong"
    if magnitude >= 5:
        return "strong"
    if magnitude >= 3:
        return "moderate"
    return "weak"


def favored_model(delta: float) -> str:
    if math.isclose(delta, 0.0, abs_tol=1e-12):
        return "tie"
    return "heterochronous" if delta > 0 else "isochronous"


def scenario_parts(name: str) -> tuple[str, str, str]:
    parts = name.split("_")
    if len(parts) < 3:
        raise ValueError(f"Unexpected BETS scenario name: {name}")
    return parts[0], parts[1], "_".join(parts[2:])


def main() -> int:
    args = parse_args()
    scenario_logs = parse_mapping(args.scenario_log)
    scenario_errs = parse_mapping(args.scenario_err)

    summary_rows: list[dict[str, str]] = []
    scenario_values: dict[str, float | None] = {}

    for scenario in sorted(set(scenario_logs) | set(scenario_errs)):
        stdout_path = scenario_logs.get(scenario)
        stderr_path = scenario_errs.get(scenario)
        stdout_text = stdout_path.read_text(encoding="utf-8", errors="replace") if stdout_path and stdout_path.exists() else ""
        stderr_text = stderr_path.read_text(encoding="utf-8", errors="replace") if stderr_path and stderr_path.exists() else ""
        value = extract_marginal_log_likelihood(stdout_text)
        if value is None:
            value = extract_marginal_log_likelihood(stderr_text)

        clock_model, tree_prior, temporal_mode = scenario_parts(scenario)
        scenario_values[scenario] = value
        summary_rows.append(
            {
                "scenario": scenario,
                "clock_model": clock_model,
                "tree_prior": tree_prior,
                "temporal_mode": temporal_mode,
                "stdout_log": str(stdout_path) if stdout_path else "",
                "stderr_log": str(stderr_path) if stderr_path else "",
                "marginal_log_likelihood": "" if value is None else f"{value:.6f}",
                "status": "parsed" if value is not None else "missing_marginal_log_likelihood",
            }
        )

    out_summary = Path(args.out_summary)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    with out_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scenario",
                "clock_model",
                "tree_prior",
                "temporal_mode",
                "stdout_log",
                "stderr_log",
                "marginal_log_likelihood",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    comparisons = [
        ("strict_constant", "strict_constant_heterochronous", "strict_constant_isochronous"),
        ("ucln_constant", "ucln_constant_heterochronous", "ucln_constant_isochronous"),
    ]

    out_bayes = Path(args.out_bayes_factors)
    out_bayes.parent.mkdir(parents=True, exist_ok=True)
    with out_bayes.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "comparison",
                "clock_model",
                "tree_prior",
                "heterochronous_scenario",
                "isochronous_scenario",
                "log_bayes_factor_heterochronous_minus_isochronous",
                "favored_temporal_model",
                "evidence_strength",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for comparison, hetero_name, iso_name in comparisons:
            hetero = scenario_values.get(hetero_name)
            iso = scenario_values.get(iso_name)
            clock_model, tree_prior, _ = scenario_parts(hetero_name)
            if hetero is None or iso is None:
                writer.writerow(
                    {
                        "comparison": comparison,
                        "clock_model": clock_model,
                        "tree_prior": tree_prior,
                        "heterochronous_scenario": hetero_name,
                        "isochronous_scenario": iso_name,
                        "log_bayes_factor_heterochronous_minus_isochronous": "",
                        "favored_temporal_model": "",
                        "evidence_strength": "",
                        "status": "missing_marginal_log_likelihood",
                    }
                )
                continue

            delta = hetero - iso
            writer.writerow(
                {
                    "comparison": comparison,
                    "clock_model": clock_model,
                    "tree_prior": tree_prior,
                    "heterochronous_scenario": hetero_name,
                    "isochronous_scenario": iso_name,
                    "log_bayes_factor_heterochronous_minus_isochronous": f"{delta:.6f}",
                    "favored_temporal_model": favored_model(delta),
                    "evidence_strength": evidence_label(delta),
                    "status": "ok",
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
