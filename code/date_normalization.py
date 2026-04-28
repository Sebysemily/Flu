import re
from datetime import datetime
from typing import Mapping, Optional, Sequence

import pandas as pd


YEAR_RE = re.compile(r"(19|20)\d{2}")


def parse_collection_date(value: str) -> Optional[str]:
	if value is None:
		return None
	text = str(value).strip()
	if not text or text.lower() in {"nan", "none", "na", "n/a"}:
		return None

	if re.match(r"^\d{4}-\d{2}-\d{2}$", text):
		return text
	if re.match(r"^\d{4}-\d{2}$", text):
		return f"{text}-01"
	if re.match(r"^\d{4}$", text):
		return f"{text}-07-01"

	for fmt in ("%d-%b-%Y", "%b-%Y"):
		try:
			parsed = datetime.strptime(text, fmt)
			if fmt == "%b-%Y":
				return parsed.strftime("%Y-%m-01")
			return parsed.strftime("%Y-%m-%d")
		except ValueError:
			pass

	for dayfirst in (False, True):
		parsed = pd.to_datetime(text, errors="coerce", dayfirst=dayfirst)
		if pd.notna(parsed):
			return parsed.date().isoformat()

	match = YEAR_RE.search(text)
	if match:
		return f"{match.group(0)}-07-01"
	return None


def extract_year(value: str) -> Optional[str]:
	parsed = parse_collection_date(value)
	if not parsed:
		return None
	return parsed[:4]


def extract_header_date(header: str) -> Optional[str]:
	if header is None:
		return None
	text = str(header).strip()
	if not text:
		return None
	return parse_collection_date(text.rsplit("/", 1)[-1].strip())


def validate_no_missing_dates(
	rows: Sequence[Mapping[str, object]],
	label_key: str,
	date_key: str,
	context: str,
	max_examples: int = 8,
) -> None:
	missing = []
	for row in rows:
		date_value = parse_collection_date(row.get(date_key))
		if date_value:
			continue
		label = str(row.get(label_key, "")).strip() or "UNKNOWN"
		raw_date = str(row.get(date_key, "")).strip()
		missing.append((label, raw_date))

	if not missing:
		return

	examples = ", ".join(
		f"{label} ({raw_date or 'sin_fecha'})"
		for label, raw_date in missing[:max_examples]
	)
	raise ValueError(
		f"Se encontraron {len(missing)} registros sin fecha valida para {context}: {examples}"
	)


def pick_ecuador_date(row: Mapping[str, object], source: str = "reception") -> Optional[str]:
	source_value = (source or "reception").strip().lower()
	if source_value == "collection":
		return parse_collection_date(
			row.get("Fecha colección") or row.get("Fecha coleccion")
		)
	if source_value == "reception":
		return parse_collection_date(
			row.get("Fecha recepción") or row.get("Fecha recepcion")
		)
	raise ValueError(
		f"ecuador_date_source invalido: {source}. Usa 'collection' o 'reception'."
	)
