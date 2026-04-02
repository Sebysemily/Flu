import re
from datetime import datetime
from typing import Mapping, Optional

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


def pick_ecuador_date(row: Mapping[str, object], source: str = "reception") -> Optional[str]:
	source_value = (source or "reception").strip().lower()
	collection_date = parse_collection_date(
		row.get("Fecha colección") or row.get("Fecha coleccion")
	)
	reception_date = parse_collection_date(
		row.get("Fecha recepción") or row.get("Fecha recepcion")
	)

	if source_value == "collection":
		return collection_date or reception_date
	return reception_date or collection_date