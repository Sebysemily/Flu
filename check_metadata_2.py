import pandas as pd

df = pd.read_csv("metadata/H5N1_context.csv")
ecuador = df[df["expected_role"].str.startswith("flu_")]
segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
print("\nEcuador incomplete samples:")
incomplete_count = 0
for idx, row in ecuador.iterrows():
    present = []
    missing = []
    for s in segments:
        val = row[s]
        if pd.notna(val) and str(val).strip() != "" and str(val).strip() != "nan":
            present.append(s)
        else:
            missing.append(s)
    if missing:
        incomplete_count += 1
        print(f"- {row['file_name']}: Missing {missing}, Present {present}")

if incomplete_count == 0:
    print("All Ecuador samples have 8 segments complete in H5N1_context.csv")
else:
    print(f"Total incomplete Ecuador samples: {incomplete_count}")

