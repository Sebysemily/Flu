"""Build GISAID context panel rows from the raw context FASTA."""

from __future__ import annotations

from build_inputs.process_raw_to_segments import (
    dedupe_metadata_rows,
    filter_complete_context_isolates,
    parse_context_isolates,
)

PANEL_COLUMNS = ["file_name", "host", "collection_date", "country", "expected_role"]


def is_ecuador_country(country: str) -> bool:
    key = str(country or "").strip().lower().replace("_", "").replace(" ", "")
    return key == "ecuador"


def context_expected_role(country: str, default_role: str) -> str:
    """Ecuadorian GISAID context isolates are coastal Ecuador panel samples."""
    if is_ecuador_country(country):
        return "flu_costa"
    return default_role


def build_gisaid_context_rows(context_fasta: str, local_epi_isls: set[str]) -> list[dict[str, str]]:
    """One metadata row per GISAID EPI_ISL in the context FASTA (non-local isolates)."""
    isolates_data = parse_context_isolates(context_fasta)
    complete_context, context_dates, context_places, context_types = filter_complete_context_isolates(
        isolates_data, local_epi_isls
    )

    rows: list[dict[str, str]] = []
    for isolate, segs in sorted(complete_context.items()):
        country = context_places[isolate]
        date_value = context_dates[isolate]
        role = context_expected_role(country, context_types[isolate])
        parts = isolate.split("/")
        host = parts[1] if (len(parts) > 1 and parts[0] == "A") else "unknown"
        for _seg, (epi_isl, _seq, _hdr) in segs.items():
            rows.append(
                {
                    "file_name": epi_isl,
                    "host": host,
                    "collection_date": date_value,
                    "country": country,
                    "expected_role": role,
                }
            )

    return dedupe_metadata_rows(rows)
