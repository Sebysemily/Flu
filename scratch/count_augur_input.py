import csv
import os

ha_strains = set()
with open('data/phylogeny/by_segment/H5N1_HA.fasta') as f:
    for line in f:
        if line.startswith('>'):
            ha_strains.add(line.strip()[1:])

flu_strains = []
with open('metadata/H5N1_context.csv', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("expected_role", "").startswith("flu_"):
            flu_strains.append(row.get("file_name"))

flu_ha = [s for s in flu_strains if s in ha_strains]
final_include = set(['EPI_ISL_18133029'] + flu_ha)

rc = []
aa = []

with open('metadata/H5N1_context.csv', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        fname = row.get("file_name")
        if fname not in final_include:
            role = row.get("expected_role")
            if role == "regional_context":
                rc.append(fname)
            elif role == "american_anchor":
                aa.append(fname)

rc_ha = [s for s in rc if s in ha_strains]
aa_ha = [s for s in aa if s in ha_strains]

print(f"Total HA sequences available overall: {len(ha_strains)}")
print(f"Protected flu strains with HA (bypassing filter): {len(flu_ha)}")
print(f"Total metadata records entering augur filter: {len(rc) + len(aa)}")
print(f"  - Regional Context entering filter: {len(rc)} (of which {len(rc_ha)} have HA sequences)")
print(f"  - American Anchor entering filter: {len(aa)} (of which {len(aa_ha)} have HA sequences)")
print(f"  - Total samples WITH HA sequences entering augur filter: {len(rc_ha) + len(aa_ha)}")
