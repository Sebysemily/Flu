#!/usr/bin/env python3
"""Sync BEAST XML MCMC settings from config/config.yml (beast.models)."""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

import yaml

DEFAULT_CONFIG = "config/config.yml"


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_xml_params(content: str) -> dict[str, int | None]:
    mcmc = re.search(r"<mcmc\b[^>]*\bchainLength=\"(\d+)\"", content)
    logs = re.findall(r"<log\b[^>]*\blogEvery=\"(\d+)\"", content)
    chk = re.search(
        r"<logCheckpoint\b[^>]*\bcheckpointFinal=\"(\d+)\"", content
    )
    chk_every = re.search(
        r"<logCheckpoint\b[^>]*\bcheckpointEvery=\"(\d+)\"", content
    )
    return {
        "chainLength": int(mcmc.group(1)) if mcmc else None,
        "logEvery": int(logs[0]) if logs else None,
        "logEvery_all": [int(x) for x in logs],
        "checkpointFinal": int(chk.group(1)) if chk else None,
        "checkpointEvery": int(chk_every.group(1)) if chk_every else None,
    }


def apply_params(
    content: str,
    chain_length: int,
    log_every: int,
    checkpoint_every: int | None = None,
) -> str:
    content = re.sub(
        r"(<mcmc\b[^>]*\bchainLength=\")(\d+)(\")",
        rf"\g<1>{chain_length}\g<3>",
        content,
        count=1,
    )
    if content.count(f'chainLength="{chain_length}"') != 1:
        raise ValueError("Expected exactly one mcmc chainLength attribute after update.")

    content = re.sub(
        r"(<log\b[^>]*\blogEvery=\")(\d+)(\")",
        rf"\g<1>{log_every}\g<3>",
        content,
    )
    content = re.sub(
        r"(<logTree\b[^>]*\blogEvery=\")(\d+)(\")",
        rf"\g<1>{log_every}\g<3>",
        content,
    )
    if checkpoint_every is None:
        checkpoint_every = max(1_000_000, chain_length // 50)
    if re.search(r"<logCheckpoint\b", content):
        content = re.sub(
            r"(<logCheckpoint\b[^>]*\bcheckpointEvery=\")(\d+)(\")",
            rf"\g<1>{checkpoint_every}\g<3>",
            content,
            count=1,
        )
        content = re.sub(
            r"(<logCheckpoint\b[^>]*\bcheckpointFinal=\")(\d+)(\")",
            rf"\g<1>{chain_length}\g<3>",
            content,
            count=1,
        )
    content = re.sub(
        r"(<log\s+[^>]*?)overwrite=\"false\"([^>]*>)",
        r"\1\2",
        content,
    )
    content = re.sub(
        r"(<logTree\s+[^>]*?)overwrite=\"false\"([^>]*>)",
        r"\1\2",
        content,
    )
    return content


def update_xml(
    xml_path: str,
    chain_length: int,
    log_every: int,
    checkpoint_every: int | None,
) -> tuple[bool, dict]:
    with open(xml_path, encoding="utf-8") as handle:
        before = handle.read()
    before_p = read_xml_params(before)
    after = apply_params(before, chain_length, log_every, checkpoint_every)
    after_p = read_xml_params(after)
    changed = after != before
    if changed:
        with open(xml_path, "w", encoding="utf-8") as handle:
            handle.write(after)
    summary = {
        "xml": xml_path,
        "before": before_p,
        "after": after_p,
        "target": {
            "chainLength": chain_length,
            "logEvery": log_every,
            "checkpointEvery": checkpoint_every or max(1_000_000, chain_length // 50),
            "checkpointFinal": chain_length,
        },
        "changed": changed,
    }
    return changed, summary


def write_stamp(path: str, summary: dict) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", nargs="?", help="Model key in config beast.models")
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stamp",
        default=None,
        help="Write JSON stamp with applied parameters (for Snakemake output)",
    )
    parser.add_argument(
        "--all-gss",
        action="store_true",
        help="Update all template_beast/*.xml that have beast.models entries",
    )
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Config file not found: {args.config}", file=sys.stderr)
        sys.exit(1)

    config = load_config(args.config)
    models = config.get("beast", {}).get("models", {})
    if not models:
        print("No beast.models section in config.", file=sys.stderr)
        sys.exit(1)

    scenarios = list(models.keys()) if args.all_gss else [args.scenario]
    if not args.all_gss and not args.scenario:
        print("Usage: update_beast_xmls.py <scenario> OR --all-gss", file=sys.stderr)
        sys.exit(1)

    any_change = False
    for scenario in scenarios:
        params = models.get(scenario)
        if not params:
            print(f"Skip {scenario}: no params in config.", file=sys.stderr)
            continue
        chain_length = params.get("chainLength")
        log_every = params.get("logEvery")
        if chain_length is None or log_every is None:
            print(f"Skip {scenario}: missing chainLength/logEvery.", file=sys.stderr)
            continue
        xml_path = f"template_beast/{scenario}.xml"
        if not os.path.exists(xml_path):
            print(f"Skip {scenario}: {xml_path} not found.", file=sys.stderr)
            continue
        changed, summary = update_xml(
            xml_path,
            int(chain_length),
            int(log_every),
            None,
        )
        any_change = any_change or changed
        status = "updated" if changed else "unchanged"
        after = summary["after"]
        print(
            f"{status} {xml_path}: "
            f"chainLength {summary['before'].get('chainLength')} -> {after.get('chainLength')}, "
            f"logEvery {summary['before'].get('logEvery')} -> {after.get('logEvery')}, "
            f"checkpointFinal {summary['before'].get('checkpointFinal')} -> {after.get('checkpointFinal')}"
        )
        if args.stamp and scenario == (args.scenario or scenario):
            write_stamp(args.stamp, summary)

    failed = False
    for scenario in scenarios:
        params = models.get(scenario)
        if not params or not os.path.exists(f"template_beast/{scenario}.xml"):
            continue
        content = open(f"template_beast/{scenario}.xml", encoding="utf-8").read()
        current = read_xml_params(content)
        want_cl = int(params["chainLength"])
        want_le = int(params["logEvery"])
        if current.get("chainLength") != want_cl or current.get("logEvery") != want_le:
            print(
                f"ERROR {scenario}.xml still has chainLength={current.get('chainLength')} "
                f"logEvery={current.get('logEvery')} (wanted {want_cl}/{want_le})",
                file=sys.stderr,
            )
            failed = True
        if current.get("logEvery_all") and len(set(current["logEvery_all"])) > 1:
            print(f"ERROR {scenario}.xml has inconsistent logEvery values", file=sys.stderr)
            failed = True
    if failed:
        sys.exit(1)
    if not any_change:
        print("All requested XML files already match config.", file=sys.stderr)


if __name__ == "__main__":
    main()
