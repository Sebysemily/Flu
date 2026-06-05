#!/usr/bin/env python3
"""Sync BEAST XML MCMC settings from config/config.yml (beast.models).

Only the main <mcmc> block is updated (chainLength, logEvery, logTree, logCheckpoint).
Marginal-likelihood / GSS settings after </mcmc> are left unchanged.
"""
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


def split_mcmc_block(content: str) -> tuple[str, str]:
    """Return (mcmc_block_including_close_tag, rest_after_mcmc)."""
    end = content.find("</mcmc>")
    if end == -1:
        raise ValueError("No </mcmc> closing tag found in BEAST XML.")
    close_end = end + len("</mcmc>")
    return content[:close_end], content[close_end:]


def read_gss_params(content: str) -> dict[str, int | None]:
    tail = split_mcmc_block(content)[1]
    mle = re.search(
        r"<marginalLikelihoodEstimator\b[^>]*\bchainLength=\"(\d+)\"", tail
    )
    mle_log = re.search(
        r"<log\s+id=\"MLELog\"[^>]*\blogEvery=\"(\d+)\"", tail
    )
    return {
        "gssChainLength": int(mle.group(1)) if mle else None,
        "gssLogEvery": int(mle_log.group(1)) if mle_log else None,
    }


def read_xml_params(content: str) -> dict[str, int | None]:
    mcmc_block, _ = split_mcmc_block(content)
    mcmc = re.search(r"<mcmc\b[^>]*\bchainLength=\"(\d+)\"", mcmc_block)
    logs = re.findall(r"<log\b[^>]*\blogEvery=\"(\d+)\"", mcmc_block)
    chk = re.search(
        r"<logCheckpoint\b[^>]*\bcheckpointFinal=\"(\d+)\"", mcmc_block
    )
    chk_every = re.search(
        r"<logCheckpoint\b[^>]*\bcheckpointEvery=\"(\d+)\"", mcmc_block
    )
    out = {
        "chainLength": int(mcmc.group(1)) if mcmc else None,
        "logEvery": int(logs[0]) if logs else None,
        "logEvery_all": [int(x) for x in logs],
        "checkpointFinal": int(chk.group(1)) if chk else None,
        "checkpointEvery": int(chk_every.group(1)) if chk_every else None,
    }
    out.update(read_gss_params(content))
    return out


def apply_mcmc_params(
    mcmc_block: str,
    chain_length: int,
    log_every: int,
    checkpoint_every: int | None = None,
) -> str:
    mcmc_block = re.sub(
        r"(<mcmc\b[^>]*\bchainLength=\")(\d+)(\")",
        rf"\g<1>{chain_length}\g<3>",
        mcmc_block,
        count=1,
    )
    if mcmc_block.count(f'chainLength="{chain_length}"') != 1:
        raise ValueError("Expected exactly one mcmc chainLength attribute after update.")

    mcmc_block = re.sub(
        r"(<log\b[^>]*\blogEvery=\")(\d+)(\")",
        rf"\g<1>{log_every}\g<3>",
        mcmc_block,
    )
    mcmc_block = re.sub(
        r"(<logTree\b[^>]*\blogEvery=\")(\d+)(\")",
        rf"\g<1>{log_every}\g<3>",
        mcmc_block,
    )
    if checkpoint_every is None:
        checkpoint_every = max(1_000_000, chain_length // 50)
    if re.search(r"<logCheckpoint\b", mcmc_block):
        mcmc_block = re.sub(
            r"(<logCheckpoint\b[^>]*\bcheckpointEvery=\")(\d+)(\")",
            rf"\g<1>{checkpoint_every}\g<3>",
            mcmc_block,
            count=1,
        )
        mcmc_block = re.sub(
            r"(<logCheckpoint\b[^>]*\bcheckpointFinal=\")(\d+)(\")",
            rf"\g<1>{chain_length}\g<3>",
            mcmc_block,
            count=1,
        )
    mcmc_block = re.sub(
        r"(<log\s+[^>]*?)overwrite=\"false\"([^>]*>)",
        r"\1\2",
        mcmc_block,
    )
    mcmc_block = re.sub(
        r"(<logTree\s+[^>]*?)overwrite=\"false\"([^>]*>)",
        r"\1\2",
        mcmc_block,
    )
    return mcmc_block


def apply_params(
    content: str,
    chain_length: int,
    log_every: int,
    checkpoint_every: int | None = None,
) -> str:
    gss_before = read_gss_params(content)
    mcmc_block, tail = split_mcmc_block(content)
    mcmc_block = apply_mcmc_params(mcmc_block, chain_length, log_every, checkpoint_every)
    updated = mcmc_block + tail
    gss_after = read_gss_params(updated)
    if gss_before != gss_after:
        raise ValueError(
            f"GSS block changed unexpectedly: before={gss_before} after={gss_after}"
        )
    return updated


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


def model_gss_flag(models: dict, scenario: str) -> bool | None:
    entry = models.get(scenario, {})
    if "gss" not in entry:
        return None
    return bool(entry["gss"])


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
        gss = summary["before"].get("gssChainLength")
        gss_note = f", GSS chainLength={gss} (unchanged)" if gss is not None else ""
        print(
            f"{status} {xml_path}: "
            f"chainLength {summary['before'].get('chainLength')} -> {after.get('chainLength')}, "
            f"logEvery {summary['before'].get('logEvery')} -> {after.get('logEvery')}, "
            f"checkpointFinal {summary['before'].get('checkpointFinal')} -> {after.get('checkpointFinal')}"
            f"{gss_note}"
        )
        if args.stamp and scenario == (args.scenario or scenario):
            gss_cfg = model_gss_flag(models, scenario)
            xml_has_gss = read_gss_params(open(xml_path, encoding="utf-8").read()).get(
                "gssChainLength"
            )
            summary["gss_enabled"] = gss_cfg
            summary["template_has_gss"] = xml_has_gss is not None
            if gss_cfg and xml_has_gss is None:
                print(
                    f"WARNING {scenario}: config gss=true but {xml_path} has no "
                    f"marginalLikelihoodEstimator block.",
                    file=sys.stderr,
                )
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
