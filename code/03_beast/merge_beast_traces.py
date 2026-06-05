#!/usr/bin/env python3
"""Merge BEAST 1 MCMC log/trees segments after resume (LogCombiner).

Keeps a single canonical H5N1_HA_panel_postQC.log and .trees in the run directory.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_LOG = "H5N1_HA_panel_postQC.log"
DEFAULT_TREES = "H5N1_HA_panel_postQC.trees"
PRIOR_LOG = f"{DEFAULT_LOG}.prior"
PRIOR_TREES = f"{DEFAULT_TREES}.prior"
NEW_LOG = f"{DEFAULT_LOG}.new"
NEW_TREES = f"{DEFAULT_TREES}.new"
MIN_MCMC_ROWS = 2


def mcmc_log_header(path: Path) -> str | None:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("state\t") or line.startswith("state "):
            return line
    return None


def mcmc_state_bounds(path: Path) -> tuple[int | None, int | None, int]:
    """Return (first_state, last_state, row_count) for MCMC sample lines."""
    first: int | None = None
    last: int | None = None
    rows = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"^([0-9]+)\t", line)
        if not m:
            continue
        state = int(m.group(1))
        rows += 1
        if first is None:
            first = state
        last = state
    return first, last, rows


def has_mcmc_samples(path: Path, min_rows: int = MIN_MCMC_ROWS) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    _, _, rows = mcmc_state_bounds(path)
    return rows >= min_rows


def trees_has_samples(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    text = path.read_text(encoding="utf-8", errors="replace")
    return "tree STATE_" in text or "tree " in text


def validate_log_merge_pair(prior: Path, new: Path) -> None:
    h_prior = mcmc_log_header(prior)
    h_new = mcmc_log_header(new)
    if not h_prior or not h_new:
        raise SystemExit(
            f"ERROR: Missing MCMC header in log merge ({prior.name} / {new.name})."
        )
    if h_prior != h_new:
        raise SystemExit(
            "ERROR: Log headers differ; cannot LogCombiner resume segments.\n"
            f"  prior: {h_prior[:80]}...\n"
            f"  new:   {h_new[:80]}..."
        )
    _, last_prior, _ = mcmc_state_bounds(prior)
    first_new, _, _ = mcmc_state_bounds(new)
    if last_prior is None or first_new is None:
        raise SystemExit("ERROR: No MCMC state rows in prior or new log.")
    if first_new < last_prior:
        raise SystemExit(
            f"ERROR: New log starts at state {first_new} but prior ends at {last_prior} "
            "(overlapping or restarted numbering). Check -load_state / -overwrite resume."
        )


def run_logcombiner(
    inputs: list[Path],
    output: Path,
    *,
    trees: bool = False,
    logcombiner: str = "logcombiner",
) -> None:
    if len(inputs) < 2:
        raise ValueError("LogCombiner needs at least two input files.")
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    cmd = [logcombiner]
    if trees:
        cmd.append("-trees")
    cmd.extend(str(p) for p in inputs)
    cmd.append(str(tmp))
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)
    if output.exists():
        output.unlink()
    tmp.rename(output)


def log_part_candidates(run_dir: Path) -> list[Path]:
    logs_dir = run_dir / "logs"
    if not logs_dir.is_dir():
        return []
    parts = list(logs_dir.glob("*.log.part_*"))
    parts.extend(logs_dir.glob("H5N1_HA_panel_postQC.log.part_*"))
    parts.extend(logs_dir.glob("*.log"))
    seen: set[Path] = set()
    unique: list[Path] = []
    for p in parts:
        rp = p.resolve()
        if rp in seen or not p.is_file():
            continue
        seen.add(rp)
        unique.append(p)
    return unique


def sort_log_paths_by_first_state(paths: list[Path]) -> list[Path]:
    scored: list[tuple[int, int, Path]] = []
    for path in paths:
        if not has_mcmc_samples(path, min_rows=1):
            continue
        first, last, rows = mcmc_state_bounds(path)
        scored.append((first or 0, last or 0, path))
    scored.sort(key=lambda t: (t[0], t[1]))
    return [p for _, _, p in scored]


def consolidate_logs_to(path: Path, sources: list[Path], logcombiner: str) -> None:
    ordered = sort_log_paths_by_first_state(sources)
    if not ordered:
        raise SystemExit("ERROR: No MCMC log parts to consolidate.")
    if len(ordered) == 1:
        shutil.copy2(ordered[0], path)
        print(f"Copied single log segment to {path}", file=sys.stderr)
        return
    run_logcombiner(ordered, path, trees=False, logcombiner=logcombiner)


def consolidate_trees_to(path: Path, sources: list[Path], logcombiner: str) -> None:
    ordered = sorted(sources, key=lambda p: p.stat().st_mtime)
    ordered = [p for p in ordered if trees_has_samples(p)]
    if not ordered:
        return
    if len(ordered) == 1:
        shutil.copy2(ordered[0], path)
        print(f"Copied single trees segment to {path}", file=sys.stderr)
        return
    run_logcombiner(ordered, path, trees=True, logcombiner=logcombiner)


def tree_part_candidates(run_dir: Path) -> list[Path]:
    logs_dir = run_dir / "logs"
    if not logs_dir.is_dir():
        return []
    parts = list(logs_dir.glob("*.trees.part_*"))
    parts.extend(logs_dir.glob(f"{DEFAULT_TREES}.part_*"))
    return [p for p in parts if p.is_file()]


def snapshot_prior(run_dir: Path, logcombiner: str = "logcombiner") -> None:
    """Save pre-extension log/trees as .prior (from cwd or legacy logs/)."""
    run_dir = run_dir.resolve()
    prior_log = run_dir / PRIOR_LOG
    prior_trees = run_dir / PRIOR_TREES
    cwd_log = run_dir / DEFAULT_LOG
    cwd_trees = run_dir / DEFAULT_TREES

    if prior_log.is_file() and has_mcmc_samples(prior_log):
        print(f"Prior log already present: {prior_log}", file=sys.stderr)
    elif has_mcmc_samples(cwd_log):
        shutil.copy2(cwd_log, prior_log)
        print(f"Snapshot prior log: {cwd_log} -> {prior_log}", file=sys.stderr)
    else:
        parts = log_part_candidates(run_dir)
        if parts:
            print(f"Consolidating {len(parts)} legacy log part(s) -> {prior_log}", file=sys.stderr)
            consolidate_logs_to(prior_log, parts, logcombiner)

    if prior_trees.is_file() and trees_has_samples(prior_trees):
        print(f"Prior trees already present: {prior_trees}", file=sys.stderr)
    elif trees_has_samples(cwd_trees):
        shutil.copy2(cwd_trees, prior_trees)
        print(f"Snapshot prior trees: {cwd_trees} -> {prior_trees}", file=sys.stderr)
    else:
        parts = tree_part_candidates(run_dir)
        if parts:
            print(f"Consolidating {len(parts)} legacy trees part(s) -> {prior_trees}", file=sys.stderr)
            consolidate_trees_to(prior_trees, parts, logcombiner)


def merge_mcmc_traces(run_dir: Path, logcombiner: str = "logcombiner") -> None:
    """Merge .prior with post-BEAST cwd log/trees into canonical outputs."""
    run_dir = run_dir.resolve()
    prior_log = run_dir / PRIOR_LOG
    prior_trees = run_dir / PRIOR_TREES
    out_log = run_dir / DEFAULT_LOG
    out_trees = run_dir / DEFAULT_TREES

    if not prior_log.is_file():
        print("No prior log snapshot; skip log merge.", file=sys.stderr)
    elif not has_mcmc_samples(out_log, min_rows=1):
        print(
            f"WARNING: Post-MCMC log {out_log} empty or missing; keeping {prior_log}.",
            file=sys.stderr,
        )
    else:
        if prior_log.resolve() == out_log.resolve():
            print("Prior and output log are the same path; skip merge.", file=sys.stderr)
        else:
            validate_log_merge_pair(prior_log, out_log)
            merged = run_dir / NEW_LOG
            run_logcombiner([prior_log, out_log], merged, trees=False, logcombiner=logcombiner)
            shutil.move(str(merged), str(out_log))
            prior_log.unlink()
            first, last, rows = mcmc_state_bounds(out_log)
            print(
                f"Merged MCMC log -> {out_log} (states {first}..{last}, {rows} rows)",
                file=sys.stderr,
            )

    if not prior_trees.is_file():
        print("No prior trees snapshot; skip trees merge.", file=sys.stderr)
    elif not trees_has_samples(out_trees):
        print(
            f"WARNING: Post-MCMC trees {out_trees} empty; keeping {prior_trees}.",
            file=sys.stderr,
        )
    else:
        merged_t = run_dir / NEW_TREES
        run_logcombiner([prior_trees, out_trees], merged_t, trees=True, logcombiner=logcombiner)
        shutil.move(str(merged_t), str(out_trees))
        prior_trees.unlink()
        print(f"Merged trees -> {out_trees}", file=sys.stderr)


def recover_incomplete_merge(run_dir: Path, logcombiner: str = "logcombiner") -> bool:
    """If a prior snapshot exists and cwd has a short new segment, try merging."""
    run_dir = run_dir.resolve()
    prior_log = run_dir / PRIOR_LOG
    out_log = run_dir / DEFAULT_LOG
    if not prior_log.is_file() or not has_mcmc_samples(out_log, min_rows=1):
        return False
    _, last_prior, rows_prior = mcmc_state_bounds(prior_log)
    first_new, last_new, rows_new = mcmc_state_bounds(out_log)
    if last_prior is not None and first_new is not None and first_new >= last_prior:
        print("Recovering incomplete merge from prior + partial new log.", file=sys.stderr)
        merge_mcmc_traces(run_dir, logcombiner=logcombiner)
        return True
    if rows_new < rows_prior:
        print("Recovering: replacing short cwd log with prior-only copy.", file=sys.stderr)
        shutil.copy2(prior_log, out_log)
        prior_log.unlink()
        return True
    return False


def checkpoint_state(run_dir: Path) -> int | None:
    """Read absolute chain state from BEAST .chkpt / .state (first ``state`` line)."""
    run_dir = run_dir.resolve()
    for pattern in ("*.chkpt", "*.state"):
        for path in sorted(run_dir.glob(pattern)):
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                m = re.match(r"^state[\t ](\d+)", line)
                if m:
                    return int(m.group(1))
    return None


def chain_state_summary(run_dir: Path) -> tuple[int, int, int]:
    """Return (logged, checkpoint, effective=max(logged, checkpoint))."""
    logged = max_logged_state(run_dir)
    chkpt = checkpoint_state(run_dir) or 0
    return logged, chkpt, max(logged, chkpt)


def max_logged_state(run_dir: Path) -> int:
    """Highest MCMC state in canonical cwd log, else legacy logs/."""
    run_dir = run_dir.resolve()
    cwd_log = run_dir / DEFAULT_LOG
    if has_mcmc_samples(cwd_log, min_rows=1):
        _, last, _ = mcmc_state_bounds(cwd_log)
        return last or 0
    best = 0
    for path in log_part_candidates(run_dir):
        if not path.is_file():
            continue
        _, last, _ = mcmc_state_bounds(path)
        if last is not None:
            best = max(best, last)
    return best


def _supplemental_log_candidates(run_dir: Path, after_state: int) -> list[Path]:
    """Log files that may continue a truncated cwd log (parts, prior, cwd)."""
    run_dir = run_dir.resolve()
    candidates: list[Path] = []
    for path in log_part_candidates(run_dir):
        if path.is_file() and path.name != DEFAULT_LOG:
            candidates.append(path)
    for name in (PRIOR_LOG, DEFAULT_LOG):
        path = run_dir / name
        if path.is_file():
            candidates.append(path)
    scored: list[tuple[int, int, Path]] = []
    for path in candidates:
        first, last, rows = mcmc_state_bounds(path)
        if rows < 1 or last is None or last <= after_state:
            continue
        if first is not None and first <= after_state:
            continue
        scored.append((first or 0, last, path))
    scored.sort(key=lambda t: (t[0], t[1]))
    return [p for _, _, p in scored]


def reconcile_checkpoint_log(
    run_dir: Path,
    target: int,
    logcombiner: str = "logcombiner",
) -> bool:
    """Best-effort merge when checkpoint reached *target* but cwd log is short."""
    run_dir = run_dir.resolve()
    logged, chkpt, effective = chain_state_summary(run_dir)
    if effective < target:
        print(
            f"Reconcile skipped: effective state {effective} < target {target}.",
            file=sys.stderr,
        )
        return False
    if logged >= target:
        print(f"Log already at {logged} >= target {target}; nothing to reconcile.", file=sys.stderr)
        return True
    if chkpt < target:
        print(
            f"Reconcile skipped: checkpoint {chkpt} < target {target}.",
            file=sys.stderr,
        )
        return False

    print(
        f"Reconciling log (logged={logged}, checkpoint={chkpt}, target={target})...",
        file=sys.stderr,
    )
    cwd_log = run_dir / DEFAULT_LOG
    if not has_mcmc_samples(cwd_log, min_rows=1):
        parts = log_part_candidates(run_dir)
        if parts:
            consolidate_logs_to(cwd_log, parts, logcombiner)
            logged = max_logged_state(run_dir)
            if logged >= target:
                return True

    supplements = _supplemental_log_candidates(run_dir, logged)
    for segment in supplements:
        try:
            validate_log_merge_pair(cwd_log, segment)
        except SystemExit as exc:
            print(f"Skipping segment {segment.name}: {exc}", file=sys.stderr)
            continue
        merged = run_dir / NEW_LOG
        run_logcombiner([cwd_log, segment], merged, trees=False, logcombiner=logcombiner)
        shutil.move(str(merged), str(cwd_log))
        _, logged, _ = mcmc_state_bounds(cwd_log)
        print(f"Merged {segment.name} into {cwd_log}; log now ends at {logged}.", file=sys.stderr)
        if logged >= chkpt:
            break

    logged = max_logged_state(run_dir)
    if logged < chkpt:
        print(
            f"WARNING: Log ends at {logged} but checkpoint is {chkpt}. "
            "MCMC is complete; Tracer/GSS use the available log samples.",
            file=sys.stderr,
        )
    return logged >= target or chkpt >= target


def ensure_combined_log(run_dir: Path, logcombiner: str = "logcombiner") -> Path | None:
    """Ensure cwd has a usable combined MCMC log (for GSS reference priors)."""
    run_dir = run_dir.resolve()
    cwd_log = run_dir / DEFAULT_LOG

    if has_mcmc_samples(cwd_log, min_rows=100):
        return cwd_log

    prior_log = run_dir / PRIOR_LOG
    if has_mcmc_samples(prior_log, min_rows=100):
        shutil.copy2(prior_log, cwd_log)
        print(f"Restored combined log from {prior_log}", file=sys.stderr)
        return cwd_log

    parts = log_part_candidates(run_dir)
    if parts:
        print("Building combined log from legacy parts for GSS.", file=sys.stderr)
        consolidate_logs_to(cwd_log, parts, logcombiner)
        return cwd_log if has_mcmc_samples(cwd_log) else None

    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--logcombiner",
        default="logcombiner",
        help="LogCombiner executable (on PATH in conda env).",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--snapshot-prior",
        action="store_true",
        help="Before MCMC extension: save/consolidate prior log/trees.",
    )
    group.add_argument(
        "--merge-mcmc",
        action="store_true",
        help="After successful MCMC extension: merge prior + new segment.",
    )
    group.add_argument(
        "--recover-incomplete-merge",
        action="store_true",
        help="On job start: finish merge if prior + partial new log exist.",
    )
    group.add_argument(
        "--ensure-combined-log",
        action="store_true",
        help="Ensure canonical MCMC log exists (GSS / Tracer).",
    )
    group.add_argument(
        "--query-state",
        action="store_true",
        help="Print logged checkpoint effective states (for Snakemake shell).",
    )
    group.add_argument(
        "--reconcile-checkpoint",
        action="store_true",
        help="Merge supplemental log segments when chkpt >= target but log is short.",
    )
    parser.add_argument(
        "--target",
        type=int,
        default=0,
        help="Target chain length (required with --reconcile-checkpoint).",
    )
    parser.add_argument(
        "--shell",
        action="store_true",
        help="With --query-state, emit LOGGED=… CHKPT=… EFFECTIVE=… for eval.",
    )
    args = parser.parse_args()
    run_dir = args.run_dir

    if args.query_state:
        logged, chkpt, effective = chain_state_summary(run_dir)
        if args.shell:
            print(f"LOGGED={logged} CHKPT={chkpt} EFFECTIVE={effective}")
        else:
            print(f"{logged} {chkpt} {effective}")
        return

    if args.reconcile_checkpoint:
        if args.target <= 0:
            raise SystemExit("ERROR: --target is required with --reconcile-checkpoint.")
        reconcile_checkpoint_log(run_dir, args.target, logcombiner=args.logcombiner)
        return

    if args.recover_incomplete_merge:
        recover_incomplete_merge(run_dir, logcombiner=args.logcombiner)
        return

    if args.snapshot_prior:
        snapshot_prior(run_dir, logcombiner=args.logcombiner)
        return

    if args.merge_mcmc:
        merge_mcmc_traces(run_dir, logcombiner=args.logcombiner)
        return

    if args.ensure_combined_log:
        path = ensure_combined_log(run_dir, logcombiner=args.logcombiner)
        if path is None:
            raise SystemExit(
                "ERROR: No combined MCMC log available. "
                "Run MCMC extension or place logs under logs/*.log.part_*."
            )
        return

    raise SystemExit("No action selected.")


if __name__ == "__main__":
    main()
