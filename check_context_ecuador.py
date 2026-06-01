import pandas as pd

df = pd.read_csv("metadata/H5N1_context.csv")

# Identify context samples (not starting with flu_)
context_df = df[~df["expected_role"].str.startswith("flu_")]

# Within context, find the ones from Ecuador
# The country column might have 'Ecuador'
ecuador_context = context_df[context_df["country"].str.contains("Ecuador", case=False, na=False)]

print(f"Total context sequences in panel: {len(context_df)}")
print(f"Ecuador sequences within context: {len(ecuador_context)}")

segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
if len(ecuador_context) > 0:
    print("\nEcuador context incomplete samples:")
    incomplete_count = 0
    for idx, row in ecuador_context.iterrows():
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
        print("All Ecuador context samples have 8 segments complete.")
    else:
        print(f"Total incomplete Ecuador context samples: {incomplete_count}")
else:
    print("\nNo Ecuador sequences found within the context data.")

