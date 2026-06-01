import os
from Bio import SeqIO

main_fasta = "data/context/gisaid_epiflu_sequence.fasta"
extra_fastas = [
    "data/context/Ecuador_extra.fasta",
    "data/context/ecuador_extra2.fasta",
    "data/context/extra_usa2.fasta"
]

# Get existing headers
existing_ids = set()
for record in SeqIO.parse(main_fasta, "fasta"):
    existing_ids.add(record.id)

print(f"Original sequences in main file: {len(existing_ids)}")

# Find new sequences
new_records = []
for ef in extra_fastas:
    for record in SeqIO.parse(ef, "fasta"):
        if record.id not in existing_ids:
            new_records.append(record)
            existing_ids.add(record.id) # Add to set to avoid duplicates among extra files

print(f"Found {len(new_records)} completely NEW sequences to add.")

# Append to main file
if new_records:
    with open(main_fasta, "a") as f:
        SeqIO.write(new_records, f, "fasta")
    print(f"Successfully appended {len(new_records)} sequences to {main_fasta}.")
    
    # Recalculate totals
    total = 0
    for _ in SeqIO.parse(main_fasta, "fasta"):
        total += 1
    print(f"New total sequences in main file: {total}")
else:
    print("No new sequences to add.")

# Delete extra files
for ef in extra_fastas:
    os.remove(ef)
    print(f"Deleted {ef}")

