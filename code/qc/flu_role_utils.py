"""Helpers to attach expected_role (flu_*) to QC discard reports."""
from __future__ import annotations

import csv
import os
from typing import Any


def load_role_map(
    metadata_path: str | None,
    id_column: str = "file_name",
    role_column: str = "expected_role",
) -> dict[str, str]:
    if not metadata_path or not os.path.exists(metadata_path):
        return {}
    roles: dict[str, str] = {}
    with open(metadata_path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            taxon = (row.get(id_column) or "").strip()
            role = (row.get(role_column) or "").strip()
            if taxon and role:
                roles[taxon] = role
    return roles


def is_flu_role(role: str) -> bool:
    return role.startswith("flu_")


def write_discarded_rows(
    path: str,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_flu_subset(
    path: str,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    flu_rows = [r for r in rows if is_flu_role(str(r.get("expected_role", "")))]
    write_discarded_rows(path, flu_rows, fieldnames)
