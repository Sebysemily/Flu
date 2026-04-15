#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a BETS BEAST XML and capture logs.")
    parser.add_argument("--xml", required=True, help="Path to the BETS XML to execute.")
    parser.add_argument("--package-dir", required=True, help="Directory for BEAST user packages.")
    parser.add_argument("--threads", type=int, default=1, help="Threads passed to BEAST.")
    parser.add_argument("--stdout", required=True, help="Path to captured stdout log.")
    parser.add_argument("--stderr", required=True, help="Path to captured stderr log.")
    parser.add_argument("--done", required=True, help="Completion sentinel.")
    return parser.parse_args()


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def append_java_tool_option(env: dict[str, str], option: str) -> None:
    current = env.get("JAVA_TOOL_OPTIONS", "").strip()
    tokens = current.split()
    if option in tokens:
        env["JAVA_TOOL_OPTIONS"] = current
        return
    env["JAVA_TOOL_OPTIONS"] = f"{current} {option}".strip()


def main() -> int:
    args = parse_args()

    xml_path = Path(args.xml).resolve()
    package_dir = Path(args.package_dir).resolve()
    stdout_path = Path(args.stdout)
    stderr_path = Path(args.stderr)
    done_path = Path(args.done)

    if not xml_path.exists():
        raise FileNotFoundError(f"BETS XML not found: {xml_path}")

    ensure_parent(stdout_path)
    ensure_parent(stderr_path)
    ensure_parent(done_path)
    package_dir.mkdir(parents=True, exist_ok=True)

    beast_bin = shutil.which("beast")
    if not beast_bin:
        raise RuntimeError("Could not find 'beast' on PATH inside the active environment.")

    env = os.environ.copy()
    append_java_tool_option(env, f"-Dbeast.user.package.dir={package_dir}")

    cmd = [
        beast_bin,
        "-overwrite",
        "-threads",
        str(args.threads),
        str(xml_path),
    ]

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        proc = subprocess.run(
            cmd,
            env=env,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
        )

    if proc.returncode != 0:
        raise RuntimeError(f"BEAST exited with code {proc.returncode} for {xml_path}")

    done_path.write_text(
        "\n".join(
            [
                f"xml\t{xml_path}",
                f"stdout\t{stdout_path.resolve()}",
                f"stderr\t{stderr_path.resolve()}",
                f"package_dir\t{package_dir}",
                f"threads\t{args.threads}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
