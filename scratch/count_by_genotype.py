import csv

complete = 0
ha_only = 0
others = 0
empty = 0

with open("metadata/H5N1_context.csv", newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("expected_role") == "regional_context":
            genotype = row.get("genotype", "").strip()
            if not genotype:
                empty += 1
                others += 1
            elif "total 8 segments" in genotype:
                complete += 1
            elif "total 1 segments" in genotype:
                ha_only += 1
            elif "Not assigned" not in genotype:
                # e.g., B3.2, B1.3 etc
                complete += 1
            else:
                others += 1

print(f"Genomas completos (8 segmentos): {complete}")
print(f"Solo la región HA (1 segmento): {ha_only}")
print(f"Otros / Vacíos: {others}")
