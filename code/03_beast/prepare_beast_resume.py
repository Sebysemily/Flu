#!/usr/bin/env python3
"""Build local.xml for BEAST resume (-load_state).

By default strips the GSS block (MCMC extension only). With --with-gss, keeps the
current template GSS section so BEAST can finish remaining MCMC then run GSS.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_beast_traces import (
    DEFAULT_LOG,
    checkpoint_state,
    ensure_combined_log,
    max_logged_state,
)

GSS_START = "<!-- START Marginal Likelihood Estimator"
GSS_END = "<!-- END Marginal Likelihood Estimator"


def extract_gss_block(content: str) -> str:
    start = content.find(GSS_START)
    if start == -1:
        return ""
    end = content.find(GSS_END, start)
    if end == -1:
        raise ValueError("GSS start marker found but end marker missing.")
    end = content.find("\n", end)
    if end == -1:
        return content[start:]
    return content[start : end + 1]


def append_gss_block(body: str, gss: str) -> str:
    if not gss:
        return body
    gss = gss.strip() + "\n"
    # Snapshot MCMC XML ends with <report>…</beast>; GSS must stay inside <beast>.
    if "<report>" in body:
        return body.replace("<report>", gss + "\n\t<report>", 1)
    if "</beast>" in body:
        return body.replace("</beast>", gss + "\n</beast>", 1)
    return body.rstrip() + "\n\n" + gss + "\n"


def append_gss_from_template(body: str, template_xml: Path) -> str:
    return append_gss_block(body, extract_gss_block(template_xml.read_text(encoding="utf-8")))


# Legacy BEAUTi exports log a shared 4D "frequencies" parameter as frequencies1..4
# and CP1+2.nu / CP3.nu instead of allNus / CP1+2.frequencies / CP3.frequencies.
_LEGACY_NU_PRIORS = """\
					<logitTransformedNormalReferencePrior fileName="H5N1_HA_panel_postQC.log" parameterColumn="CP1+2.nu" burnin="1000000">
						<parameter idref="CP1+2.nu"/>
					</logitTransformedNormalReferencePrior>
					<logitTransformedNormalReferencePrior fileName="H5N1_HA_panel_postQC.log" parameterColumn="CP3.nu" burnin="1000000">
						<parameter idref="CP3.nu"/>
					</logitTransformedNormalReferencePrior>"""

# parameterColumn must be the base name "frequencies" (BEAST appends 1..dimension);
# "frequencies1" + dimension=4 incorrectly looks for frequencies11..frequencies14.
_LEGACY_FREQ_PRIOR = """\
					<logitTransformedNormalReferencePrior fileName="H5N1_HA_panel_postQC.log" parameterColumn="frequencies" dimension="4" burnin="1000000">
						<parameter idref="frequencies"/>
					</logitTransformedNormalReferencePrior>"""


def adapt_gss_for_legacy_log(gss: str) -> str:
    gss = re.sub(
        r"\s*<logitTransformedNormalReferencePrior fileName=\"H5N1_HA_panel_postQC\.log\""
        r" parameterColumn=\"allNus\"[^>]*>.*?</logitTransformedNormalReferencePrior>\s*",
        "\n" + _LEGACY_NU_PRIORS + "\n",
        gss,
        count=1,
        flags=re.DOTALL,
    )
    gss = re.sub(
        r"\s*<logitTransformedNormalReferencePrior fileName=\"H5N1_HA_panel_postQC\.log\""
        r" parameterColumn=\"CP1\+2\.frequencies\"[^>]*>.*?</logitTransformedNormalReferencePrior>\s*"
        r"\s*<logitTransformedNormalReferencePrior fileName=\"H5N1_HA_panel_postQC\.log\""
        r" parameterColumn=\"CP3\.frequencies\"[^>]*>.*?</logitTransformedNormalReferencePrior>\s*",
        "\n" + _LEGACY_FREQ_PRIOR + "\n",
        gss,
        count=1,
        flags=re.DOTALL,
    )
    return gss


def append_gss_for_run(
    body: str,
    template_xml: Path,
    *,
    use_legacy_log_columns: bool,
) -> str:
    gss = extract_gss_block(template_xml.read_text(encoding="utf-8"))
    if use_legacy_log_columns:
        gss = adapt_gss_for_legacy_log(gss)
        print(
            "GSS reference priors adapted for legacy log columns "
            "(frequencies1-4, CP1+2.nu, CP3.nu).",
            file=sys.stderr,
        )
    return append_gss_block(body, gss)


def disable_mcmc_file_logs(content: str) -> str:
    """GSS-only: do not let the 1-step MCMC pass overwrite H5N1_HA_panel_postQC.log/trees."""
    content = re.sub(
        r"\s*<log\s+id=\"fileLog\"[^>]*>.*?</log>\s*",
        "\n",
        content,
        count=1,
        flags=re.DOTALL,
    )
    content = re.sub(
        r"\s*<logTree\s+id=\"treeFileLog\"[^>]*>.*?</logTree>\s*",
        "\n",
        content,
        count=1,
        flags=re.DOTALL,
    )
    content = re.sub(
        r"\s*<logCheckpoint\b[^>]*>.*?</logCheckpoint>\s*",
        "\n",
        content,
        count=1,
        flags=re.DOTALL,
    )
    return content


def count_mcmc_log_rows(path: Path) -> int:
    n = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if re.match(r"^[0-9]+\t", line):
            n += 1
    return n


def best_mcmc_log_path(
    run_dir: Path,
    min_rows: int = 100,
    target_state: int | None = None,
) -> Path | None:
    """Best MCMC log for GSS reference priors (prefer combined cwd log, then legacy parts)."""
    scored: list[tuple[int, int, Path]] = []
    candidates: list[Path] = [run_dir / DEFAULT_LOG]
    logs_dir = run_dir / "logs"
    if logs_dir.is_dir():
        candidates.extend(logs_dir.glob("*.log.part_*"))
        candidates.extend(logs_dir.glob("H5N1_HA_panel_postQC.log.part_*"))
        candidates.extend(logs_dir.glob("*.log"))
    for path in candidates:
        if not path.is_file():
            continue
        rows = 0
        last = 0
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if re.match(r"^[0-9]+\t", line):
                rows += 1
                last = max(last, int(line.split("\t", 1)[0]))
        if rows < min_rows:
            continue
        scored.append((rows, last, path))
    if not scored:
        return None
    if target_state is not None:
        at_target = [s for s in scored if s[1] >= target_state]
        if at_target:
            scored = at_target
    scored.sort(key=lambda t: (t[0], t[1]), reverse=True)
    return scored[0][2]


def strip_gss_block(content: str) -> str:
    start = content.find(GSS_START)
    if start == -1:
        return content
    end = content.find(GSS_END, start)
    if end == -1:
        raise ValueError("GSS start marker found but end marker missing.")
    end = content.find("\n", end)
    if end == -1:
        end = len(content)
    else:
        end += 1
    return content[:start] + content[end:]


def set_mcmc_chain_length(content: str, chain_length: int) -> str:
    updated, n = re.subn(
        r"(<mcmc\b[^>]*\bchainLength=\")(\d+)(\")",
        rf"\g<1>{chain_length}\g<3>",
        content,
        count=1,
    )
    if n != 1:
        raise ValueError("Expected exactly one <mcmc chainLength=...> to update.")
    return updated


def set_checkpoint_final(content: str, final: int) -> str:
    updated, n = re.subn(
        r"(<logCheckpoint\b[^>]*\bcheckpointFinal=\")(\d+)(\")",
        rf"\g<1>{final}\g<3>",
        content,
        count=1,
    )
    if n != 1:
        raise ValueError("Expected exactly one logCheckpoint checkpointFinal=... to update.")
    return updated


def mcmc_chain_length_for_resume(
    run_dir: Path,
    target: int,
    *,
    gss_only: bool,
) -> tuple[int, int]:
    """Return (mcmc chainLength, logCheckpoint checkpointFinal).

    With -load_state, BEAST treats <mcmc chainLength> as *additional* states on top of
    the checkpoint state, not the absolute final state number. Use remaining = target - chkpt.
    """
    logged = max_logged_state(run_dir)
    chkpt = checkpoint_state(run_dir) or 0
    if gss_only:
        # -load_state treats chainLength as *additional* MCMC steps, so use 1 to skip MCMC.
        end = max(logged, chkpt, target)
        print(
            f"GSS-only (logged={logged}, checkpoint={chkpt}, target={target}): "
            f"mcmc chainLength=1, then marginalLikelihoodEstimator.",
            file=sys.stderr,
        )
        return 1, end
    if has_checkpoint(run_dir) and chkpt > 0:
        if chkpt >= target:
            print(
                f"Checkpoint already at {chkpt} >= target {target}; "
                "no MCMC extension (chainLength=0).",
                file=sys.stderr,
            )
            return 0, target
        remaining = target - chkpt
        print(
            f"Resume from checkpoint {chkpt} (logged {logged}): chainLength={remaining} "
            f"(additional states; absolute target {target}).",
            file=sys.stderr,
        )
        return remaining, target
    if logged > 0 and logged < target:
        remaining = target - logged
        print(
            f"Resume from logged state {logged}: chainLength={remaining} "
            f"(additional states; absolute target {target}).",
            file=sys.stderr,
        )
        return remaining, target
    return target, target


def log_header(path: Path) -> str:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("state\t") or line.startswith("state "):
            return line
    return ""


def needs_legacy_snapshot(archived_log: Path, template_xml: Path) -> bool:
    """Log columns from an older BEAUTi export (frequencies1..4, not CP1+2.frequencies1)."""
    header = log_header(archived_log)
    if not header:
        return False
    if "frequencies1" in header and "CP1+2.frequencies1" not in header:
        return True
    return False


def template_uses_legacy_frequency_layout(template_xml: Path) -> bool:
    """Template still uses shared ``frequencies`` (matches frequencies1..4 log columns)."""
    text = template_xml.read_text(encoding="utf-8", errors="replace")
    return 'id="frequencies"' in text and 'id="CP1+2.frequencies"' not in text


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def has_checkpoint(run_dir: Path) -> bool:
    return bool(list(run_dir.glob("*.chkpt")) or list(run_dir.glob("*.state")))


def checkpoint_compat_xml(scenario: str) -> Path:
    return repo_root() / f"template_beast/snapshots/{scenario}_checkpoint_compat.xml"


def pick_source_xml(
    scenario: str,
    template_xml: Path,
    run_dir: Path,
    mcmc_log: Path | None,
) -> Path:
    snapshot = run_dir / "mcmc_source.xml"
    if snapshot.is_file():
        return snapshot
    legacy = checkpoint_compat_xml(scenario)
    # Resume -load_state must match the XML that created the .chkpt, not necessarily
    # the current template (parameter ids/dimensions may have changed since 50M).
    if has_checkpoint(run_dir) and legacy.is_file():
        print(
            f"Using checkpoint-compatible snapshot {legacy} "
            f"(existing checkpoint in {run_dir}).",
            file=sys.stderr,
        )
        return legacy
    if has_checkpoint(run_dir) and template_uses_legacy_frequency_layout(template_xml):
        print(
            f"Using {template_xml} (checkpoint; template matches legacy log layout).",
            file=sys.stderr,
        )
        return template_xml
    if mcmc_log and mcmc_log.is_file() and needs_legacy_snapshot(mcmc_log, template_xml):
        if legacy.is_file():
            print(
                f"Using checkpoint-compatible snapshot {legacy} "
                f"(log columns differ from current template).",
                file=sys.stderr,
            )
            return legacy
        if template_uses_legacy_frequency_layout(template_xml):
            print(
                f"Using {template_xml} (legacy log columns match template layout).",
                file=sys.stderr,
            )
            return template_xml
        raise SystemExit(
            f"ERROR: {scenario} checkpoint was created with an older operator layout "
            f"(log has legacy frequency columns). Add {legacy} or keep mcmc_source.xml "
            f"from the original run."
        )
    return template_xml


def latest_archived_log(run_dir: Path) -> Path | None:
    logs_dir = run_dir / "logs"
    if not logs_dir.is_dir():
        return None
    parts = sorted(logs_dir.glob("*.log.part_*"), key=lambda p: p.stat().st_mtime)
    return parts[-1] if parts else None


def mcmc_log_for_resume(run_dir: Path) -> Path | None:
    """Best available MCMC log to detect legacy column layout (archived or cwd)."""
    arch = latest_archived_log(run_dir)
    if arch and arch.is_file():
        return arch
    cwd = run_dir / "H5N1_HA_panel_postQC.log"
    return cwd if cwd.is_file() else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--template-xml", type=Path, required=True)
    parser.add_argument("--target-chain-length", type=int, required=True)
    parser.add_argument("--output", type=Path, default=Path("local.xml"))
    parser.add_argument(
        "--with-gss",
        action="store_true",
        help="Append GSS block (use with --gss-only after MCMC target is reached).",
    )
    parser.add_argument(
        "--gss-only",
        action="store_true",
        help=(
            "Freeze MCMC at the last logged state and run GSS only. "
            "Requires last state >= --target-chain-length."
        ),
    )
    args = parser.parse_args()

    if args.gss_only and not args.with_gss:
        raise SystemExit("ERROR: --gss-only requires --with-gss.")

    template_text = args.template_xml.read_text(encoding="utf-8")
    if args.with_gss and "marginalLikelihoodEstimator" not in template_text:
        raise SystemExit(
            f"ERROR: --with-gss requested but {args.template_xml} has no "
            "marginalLikelihoodEstimator block. Export GSS from BEAUTi or set gss: false in config."
        )

    mcmc_log = mcmc_log_for_resume(args.run_dir)
    source = pick_source_xml(args.scenario, args.template_xml, args.run_dir, mcmc_log)
    from_compat = "checkpoint_compat" in source.name
    legacy_log = from_compat or bool(
        mcmc_log
        and mcmc_log.is_file()
        and needs_legacy_snapshot(mcmc_log, args.template_xml)
    )
    content = source.read_text(encoding="utf-8")
    content = strip_gss_block(content)

    mcmc_length, checkpoint_final = mcmc_chain_length_for_resume(
        args.run_dir,
        args.target_chain_length,
        gss_only=args.gss_only,
    )

    if args.gss_only:
        logged = max_logged_state(args.run_dir)
        chkpt = checkpoint_state(args.run_dir) or 0
        effective = max(logged, chkpt)
        if effective < args.target_chain_length:
            raise SystemExit(
                f"ERROR: MCMC not at target (logged={logged}, checkpoint={chkpt} "
                f"< {args.target_chain_length}). Extend MCMC without --with-gss first."
            )
        last_state = logged
        dest = args.run_dir / DEFAULT_LOG
        combined = ensure_combined_log(args.run_dir)
        if combined is None:
            ref_log = best_mcmc_log_path(
                args.run_dir,
                target_state=args.target_chain_length,
            )
            if ref_log is None:
                raise SystemExit(
                    "ERROR: No MCMC log with enough samples for GSS reference priors. "
                    "Run MCMC extension with log merge or place a combined log in cwd."
                )
            if ref_log.resolve() != dest.resolve():
                shutil.copy2(ref_log, dest)
                print(
                    f"Restored {dest.name} from {ref_log} "
                    f"({count_mcmc_log_rows(ref_log)} MCMC rows for GSS).",
                    file=sys.stderr,
                )
        else:
            rows = count_mcmc_log_rows(dest)
            print(
                f"Using combined MCMC log {dest.name} ({rows} rows) for GSS.",
                file=sys.stderr,
            )
        content = set_mcmc_chain_length(content, mcmc_length)
        content = set_checkpoint_final(content, checkpoint_final)
        content = disable_mcmc_file_logs(content)
        content = append_gss_for_run(
            content,
            args.template_xml,
            use_legacy_log_columns=legacy_log,
        )
        mode = f"GSS-only (mcmc chainLength={mcmc_length}, logged state {last_state})"
    else:
        content = set_mcmc_chain_length(content, mcmc_length)
        content = set_checkpoint_final(content, checkpoint_final)
        if args.with_gss:
            print(
                "WARNING: --with-gss without --gss-only is discouraged on resume; "
                "BEAST may not stop MCMC at target and GSS may never start. "
                "Use MCMC-only extension, then --gss-only.",
                file=sys.stderr,
            )
            content = append_gss_for_run(
                content,
                args.template_xml,
                use_legacy_log_columns=legacy_log,
            )
        mode = "MCMC+GSS" if args.with_gss else "MCMC only"

    args.output.write_text(content, encoding="utf-8")
    print(f"Wrote {args.output} from {source} ({mode})")


if __name__ == "__main__":
    main()
