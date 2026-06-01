import pandas as pd

df = pd.read_csv("metadata/H5N1_context.csv")
print("Total samples in panel (H5N1_context.csv):", len(df))
print("Context samples in panel:", len(df[~df["expected_role"].str.startswith("flu_")]))
print("Ecuador samples in panel:", len(df[df["expected_role"].str.startswith("flu_")]))

ecuador = df[df["expected_role"].str.startswith("flu_")]
segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
print("\nEcuador incomplete samples:")
incomplete_count = 0
for idx, row in ecuador.iterrows():
    present = []
    missing = []
    for s in segments:
        # Some metadata versions use '{segment}_seq' or just check fasta headers
        col = f"{s}_length"
        if col in df.columns:
            if pd.notna(row[col]) and row[col] > 0:
                present.append(s)
            else:
                missing.append(s)
    if missing:
        incomplete_count += 1
        print(f"- {row['file_name']}: Missing {missing}, Present {present}")

if incomplete_count == 0:
    print("All Ecuador samples have 8 segments complete in H5N1_context.csv")

if "HA_length" not in df.columns:
    print("Columns in metadata:", df.columns.tolist())

# Check data/context fasta for raw context counts
import subprocess
try:
    context_count = subprocess.check_output("grep -c '^>' data/context/gisaid_epiflu_sequence.fasta", shell=True).decode().strip()
    print("\nTotal context sequences in data/context/gisaid_epiflu_sequence.fasta:", context_count)
except:
    print("\nCould not read data/context/gisaid_epiflu_sequence.fasta")
