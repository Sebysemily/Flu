#!/usr/bin/env python3
"""Tests for merge_beast_traces (no logcombiner required for validation helpers)."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_beast_traces import (  # noqa: E402
    DEFAULT_LOG,
    chain_state_summary,
    checkpoint_state,
    mcmc_state_bounds,
    merge_mcmc_traces,
    snapshot_prior,
    validate_log_merge_pair,
)

HEADER = "state\tposterior\tlikelihood\tprior\n"


def write_log(path: Path, states: list[int]) -> None:
    lines = [HEADER]
    for s in states:
        lines.append(f"{s}\t-100\t-100\t-100\n")
    path.write_text("".join(lines), encoding="utf-8")


class TestLogMergeValidation(unittest.TestCase):
    def test_validate_sequential_segments(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            prior = d / "a.log"
            new = d / "b.log"
            write_log(prior, [0, 50, 100])
            write_log(new, [100, 150, 200])
            validate_log_merge_pair(prior, new)

    def test_validate_rejects_overlap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            prior = d / "a.log"
            new = d / "b.log"
            write_log(prior, [0, 100])
            write_log(new, [50, 200])
            with self.assertRaises(SystemExit):
                validate_log_merge_pair(prior, new)

    def test_state_bounds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "x.log"
            write_log(p, [10, 20, 30])
            first, last, rows = mcmc_state_bounds(p)
            self.assertEqual((first, last, rows), (10, 30, 3))


class TestMergeMcmcTraces(unittest.TestCase):
    def test_merge_calls_logcombiner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            prior = d / f"{DEFAULT_LOG}.prior"
            new = d / DEFAULT_LOG
            write_log(prior, [0, 100])
            write_log(new, [100, 200])

            def fake_logcombiner(inputs, output, *, trees=False, logcombiner="logcombiner"):
                output.write_text(prior.read_text() + new.read_text(), encoding="utf-8")

            with mock.patch(
                "merge_beast_traces.run_logcombiner", side_effect=fake_logcombiner
            ):
                merge_mcmc_traces(d)
            self.assertFalse(prior.exists())
            _, last, rows = mcmc_state_bounds(d / DEFAULT_LOG)
            self.assertEqual(last, 200)
            self.assertEqual(rows, 4)


class TestSnapshotPrior(unittest.TestCase):
    def test_snapshot_copies_cwd_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_log(d / DEFAULT_LOG, [0, 1000])
            snapshot_prior(d)
            prior = d / f"{DEFAULT_LOG}.prior"
            self.assertTrue(prior.is_file())
            self.assertEqual(mcmc_state_bounds(prior)[1], 1000)


class TestCheckpointState(unittest.TestCase):
    def test_reads_chkpt_state_line(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "run.chkpt").write_text("state\t75000000\nrng\t123\n", encoding="utf-8")
            write_log(d / DEFAULT_LOG, [5000, 10000])
            self.assertEqual(checkpoint_state(d), 75_000_000)
            logged, chkpt, effective = chain_state_summary(d)
            self.assertEqual((logged, chkpt, effective), (10_000, 75_000_000, 75_000_000))


if __name__ == "__main__":
    unittest.main()
