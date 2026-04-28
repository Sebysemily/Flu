#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_SUFFIXES = [".log", ".trees", ".chkpt", ".ops", ".Lugar.rates.log"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one BEAST replicate and isolate its generated files.")
    parser.add_argument("--xml", required=True)
    parser.add_argument("--beast-binary", default="")
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--replicate-dir", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--beagle-mode", default="off")
    parser.add_argument("--beagle-resource", default="auto")
    parser.add_argument("--beagle-vendor", default="any")
    parser.add_argument("--beagle-platform", default="auto")
    parser.add_argument("--beagle-precision", default="auto")
    parser.add_argument("--beagle-scaling", default="default")
    parser.add_argument("--beagle-info", default="false")
    parser.add_argument("--beagle-threads", default="")
    parser.add_argument("--beagle-fallback-to-cpu", default="false")
    parser.add_argument("--status", required=True)
    parser.add_argument("--stdout", required=True)
    parser.add_argument("--stderr", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--heartbeat-seconds", type=int, default=60)
    return parser.parse_args()


def parse_bool(raw: str) -> bool:
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def normalize_choice(raw: str, allowed: set[str], default: str) -> str:
    value = str(raw).strip().lower()
    if not value:
        return default
    if value not in allowed:
        raise ValueError(f"Valor no soportado: {raw!r}. Opciones validas: {sorted(allowed)}")
    return value


def detect_beagle_version_text(proc: subprocess.CompletedProcess[str]) -> str:
    return "\n".join(part for part in [proc.stdout.strip(), proc.stderr.strip()] if part).strip()


def query_beagle_resources(beast_binary: str) -> tuple[list[dict[str, object]], str]:
    try:
        proc = subprocess.run(
            [beast_binary, "-beagle_info"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError:
        return [], ""

    text = detect_beagle_version_text(proc)
    resources: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    for raw_line in text.splitlines():
        line = raw_line.strip()
        match = re.match(r"^(\d+)\s*:\s*(.+)$", line)
        if match:
            label = match.group(2).strip()
            upper_label = label.upper()
            current = {
                "index": int(match.group(1)),
                "label": label,
                "is_gpu": "GPU" in upper_label,
                "is_cpu": "CPU" in upper_label,
                "has_opencl": False,
                "has_cuda": False,
            }
            resources.append(current)
            continue

        if current is None or not line.startswith("Flags:"):
            continue

        flags = set(line.replace("Flags:", "", 1).strip().split())
        current["flags"] = sorted(flags)
        current["has_opencl"] = "FRAMEWORK_OPENCL" in flags
        current["has_cuda"] = "FRAMEWORK_CUDA" in flags

    return resources, text


def choose_beagle_platform(args: argparse.Namespace) -> str:
    vendor = normalize_choice(args.beagle_vendor, {"any", "amd", "nvidia"}, "any")
    platform = normalize_choice(args.beagle_platform, {"auto", "cuda", "opencl"}, "auto")
    if platform != "auto":
        return platform
    if vendor == "amd":
        return "opencl"
    if vendor == "nvidia":
        return "cuda"
    return "auto"


def resource_matches(
    resource: dict[str, object],
    desired_resource: str,
    desired_platform: str,
) -> bool:
    is_gpu = bool(resource.get("is_gpu", False))
    is_cpu = bool(resource.get("is_cpu", False))

    if desired_resource == "gpu" and not is_gpu:
        return False
    if desired_resource == "cpu" and not is_cpu:
        return False

    if desired_platform == "opencl" and not bool(resource.get("has_opencl", False)):
        return False
    if desired_platform == "cuda" and not bool(resource.get("has_cuda", False)):
        return False

    return True


def build_beagle_args(args: argparse.Namespace, beast_binary: str) -> tuple[list[str], str, str]:
    mode = normalize_choice(args.beagle_mode, {"off", "auto", "force"}, "off")
    if mode == "off":
        return [], "off", ""

    desired_resource = normalize_choice(args.beagle_resource, {"auto", "cpu", "gpu"}, "auto")
    desired_platform = choose_beagle_platform(args)
    precision = normalize_choice(args.beagle_precision, {"auto", "single", "double"}, "auto")
    scaling = normalize_choice(
        args.beagle_scaling,
        {"default", "dynamic", "delayed", "always", "none"},
        "default",
    )
    fallback_to_cpu = parse_bool(args.beagle_fallback_to_cpu)

    resources, beagle_info_text = query_beagle_resources(beast_binary)
    selected_resource = desired_resource
    if mode == "auto":
        if desired_resource == "auto":
            gpu_match_exists = any(
                resource_matches(resource, "gpu", desired_platform) for resource in resources
            )
            cpu_match_exists = any(
                resource_matches(resource, "cpu", desired_platform if desired_platform == "auto" else "auto")
                for resource in resources
            )
            if gpu_match_exists:
                selected_resource = "gpu"
            elif fallback_to_cpu and cpu_match_exists:
                selected_resource = "cpu"
            else:
                return [], "auto_no_match", beagle_info_text
        else:
            match_exists = any(
                resource_matches(resource, desired_resource, desired_platform) for resource in resources
            )
            if not match_exists:
                if desired_resource == "gpu" and fallback_to_cpu:
                    cpu_match_exists = any(
                        resource_matches(
                            resource,
                            "cpu",
                            desired_platform if desired_platform == "auto" else "auto",
                        )
                        for resource in resources
                    )
                    if cpu_match_exists:
                        selected_resource = "cpu"
                    else:
                        return [], "auto_no_match", beagle_info_text
                else:
                    return [], "auto_no_match", beagle_info_text

    beagle_args = ["-beagle"]
    resource = normalize_choice(args.beagle_resource, {"auto", "cpu", "gpu"}, "auto")
    if mode == "force" and resource == "auto":
        beagle_args.append("-beagle_auto")
    elif selected_resource == "auto":
        beagle_args.append("-beagle_auto")
    elif selected_resource == "cpu":
        beagle_args.append("-beagle_CPU")
    elif selected_resource == "gpu":
        beagle_args.append("-beagle_GPU")

    if desired_platform == "cuda":
        beagle_args.append("-beagle_cuda")
    elif desired_platform == "opencl":
        beagle_args.append("-beagle_opencl")

    if precision == "single":
        beagle_args.append("-beagle_single")
    elif precision == "double":
        beagle_args.append("-beagle_double")

    if scaling != "default":
        beagle_args.extend(["-beagle_scaling", scaling])

    beagle_threads = str(args.beagle_threads).strip()
    if beagle_threads:
        beagle_args.extend(["-beagle_threads", beagle_threads])

    if parse_bool(args.beagle_info):
        beagle_args.append("-beagle_info")

    return beagle_args, mode, beagle_info_text


def resolve_beast_binary(configured: str) -> str:
    candidate = (configured or "").strip()
    if candidate:
        return candidate

    discovered = shutil.which("beast")
    if discovered:
        return discovered

    raise RuntimeError(
        "No se encontro el binario 'beast'. "
        "Instala 'beast' en envs/03_beast.yml o define config.beast.binary "
        "con la ruta al ejecutable compatible con BEAST 1."
    )


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def cleanup_prefix(prefix: Path) -> None:
    for suffix in EXPECTED_SUFFIXES:
        candidate = Path(f"{prefix}{suffix}")
        if candidate.exists():
            candidate.unlink()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def emit_progress(message: str) -> None:
    sys.stderr.write(message + "\n")
    sys.stderr.flush()


def detect_beast_version(beast_binary: str) -> str:
    commands = [
        [beast_binary, "-version"],
        [beast_binary, "--version"],
    ]
    for cmd in commands:
        try:
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        except OSError:
            continue

        text = "\n".join(part for part in [proc.stdout.strip(), proc.stderr.strip()] if part).strip()
        if text:
            return text.splitlines()[0].strip()
    return "UNKNOWN"


def main() -> None:
    args = parse_args()

    xml_path = Path(args.xml).resolve()
    output_prefix = Path(args.output_prefix)
    replicate_dir = Path(args.replicate_dir)
    status_path = Path(args.status)
    stdout_path = Path(args.stdout)
    stderr_path = Path(args.stderr)
    done_path = Path(args.done)
    run_label = f"{args.scenario}/{args.replicate}"

    if not xml_path.exists():
        raise FileNotFoundError(f"No existe el XML preparado: {xml_path}")

    beast_binary = resolve_beast_binary(args.beast_binary)
    beast_version = detect_beast_version(beast_binary)
    beagle_args, beagle_mode, beagle_info_text = build_beagle_args(args, beast_binary)

    replicate_dir.mkdir(parents=True, exist_ok=True)
    ensure_parent(status_path)
    ensure_parent(stdout_path)
    ensure_parent(stderr_path)
    ensure_parent(done_path)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    cleanup_prefix(output_prefix)

    cmd = [
        beast_binary,
        "-threads",
        str(args.threads),
        *beagle_args,
        "-seed",
        str(args.seed),
        str(xml_path),
    ]

    heartbeat_seconds = max(5, args.heartbeat_seconds)

    with status_path.open("w", encoding="utf-8") as status_handle, stdout_path.open(
        "w", encoding="utf-8"
    ) as stdout_handle, stderr_path.open("w", encoding="utf-8") as stderr_handle:
        status_handle.write(f"[{utc_now()}] run={run_label} status=starting\n")
        status_handle.write(f"[{utc_now()}] xml={xml_path}\n")
        status_handle.write(f"[{utc_now()}] beast_binary={beast_binary}\n")
        status_handle.write(f"[{utc_now()}] beast_version={beast_version}\n")
        status_handle.write(f"[{utc_now()}] seed={args.seed}\n")
        status_handle.write(f"[{utc_now()}] threads={args.threads}\n")
        status_handle.write(f"[{utc_now()}] beagle_mode={beagle_mode}\n")
        status_handle.write(f"[{utc_now()}] beagle_args={' '.join(beagle_args) if beagle_args else 'NONE'}\n")
        if beagle_info_text:
            compact_info = " | ".join(line.strip() for line in beagle_info_text.splitlines() if line.strip())
            status_handle.write(f"[{utc_now()}] beagle_info={compact_info}\n")
        status_handle.write(f"[{utc_now()}] replicate_dir={replicate_dir.resolve()}\n")
        status_handle.write(f"[{utc_now()}] command={' '.join(cmd)}\n")
        status_handle.flush()
        emit_progress(
            f"[{utc_now()}] run={run_label} status=starting seed={args.seed} threads={args.threads} beast_version={beast_version}"
        )

        stdout_handle.write(f"[{utc_now()}] starting BEAST replicate\n")
        stdout_handle.flush()
        stderr_handle.write(f"[{utc_now()}] starting BEAST replicate\n")
        stderr_handle.flush()

        proc = subprocess.Popen(cmd, stdout=stdout_handle, stderr=stderr_handle)
        status_handle.write(f"[{utc_now()}] run={run_label} status=running pid={proc.pid}\n")
        status_handle.flush()
        emit_progress(f"[{utc_now()}] run={run_label} status=running pid={proc.pid}")

        while True:
            returncode = proc.poll()
            if returncode is not None:
                break
            heartbeat = (
                f"[{utc_now()}] run={run_label} status=running pid={proc.pid} seed={args.seed}\n"
            )
            status_handle.write(heartbeat)
            status_handle.flush()
            stderr_handle.write(heartbeat)
            stderr_handle.flush()
            emit_progress(heartbeat.rstrip())
            time.sleep(heartbeat_seconds)

        stdout_handle.write(f"[{utc_now()}] BEAST finished with exit_code={returncode}\n")
        stdout_handle.flush()
        stderr_handle.write(f"[{utc_now()}] BEAST finished with exit_code={returncode}\n")
        stderr_handle.flush()
        status_handle.write(f"[{utc_now()}] run={run_label} status=finished exit_code={returncode}\n")
        status_handle.flush()
        emit_progress(f"[{utc_now()}] run={run_label} status=finished exit_code={returncode}")

    if returncode != 0:
        raise RuntimeError(f"BEAST termino con codigo {returncode} para {xml_path}")

    moved = []
    missing = []
    for suffix in EXPECTED_SUFFIXES:
        source = Path(f"{output_prefix}{suffix}")
        target = replicate_dir / f"{output_prefix.name}{suffix}"
        if source.exists():
            if target.exists():
                target.unlink()
            shutil.move(str(source), str(target))
            moved.append(target.resolve())
        else:
            missing.append(suffix)

    with done_path.open("w", encoding="utf-8") as handle:
        handle.write(f"xml\t{xml_path}\n")
        handle.write(f"scenario\t{args.scenario}\n")
        handle.write(f"replicate\t{args.replicate}\n")
        handle.write(f"beast_binary\t{beast_binary}\n")
        handle.write(f"beast_version\t{beast_version}\n")
        handle.write(f"seed\t{args.seed}\n")
        handle.write(f"threads\t{args.threads}\n")
        handle.write(f"beagle_mode\t{beagle_mode}\n")
        handle.write(f"beagle_args\t{' '.join(beagle_args) if beagle_args else 'NONE'}\n")
        handle.write(f"status\t{status_path.resolve()}\n")
        handle.write(f"stdout\t{stdout_path.resolve()}\n")
        handle.write(f"stderr\t{stderr_path.resolve()}\n")
        handle.write(f"replicate_dir\t{replicate_dir.resolve()}\n")
        handle.write(f"output_prefix\t{output_prefix}\n")
        handle.write(f"moved_outputs\t{';'.join(str(path) for path in moved)}\n")
        handle.write(f"missing_outputs\t{';'.join(missing)}\n")


if __name__ == "__main__":
    main()
