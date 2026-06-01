import pandas as pd
df = pd.read_csv("metadata/context_base.csv")
ec = df[df["country"].str.contains("Ecuador", case=False, na=False)]
print(f"Total sequences in context_base.csv: {len(df)}")
print(f"Ecuador sequences in context_base.csv: {len(ec)}")

segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
if len(ec) > 0:
    incomplete = 0
    for idx, row in ec.iterrows():
        present, missing = [], []
        for s in segments:
            val = None
            if s in df.columns:
                val = row[s]
            elif f"{s}_seq" in df.columns:
                val = row[f"{s}_seq"]
            elif "Isolate_Id" in df.columns:
                # GISAID format doesn't have explicit columns sometimes, wait.
                pass
            
            if pd.notna(val) and str(val).strip() != "" and str(val).strip() != "nan":
                present.append(s)
            else:
                missing.append(s)
        if missing:
            incomplete += 1
            name = row.get("EPI_ISL", row.get("file_name", "Unknown"))
            print(f"- {name}: Missing {missing}, Present {present}")
    print(f"Incomplete Ecuador context: {incomplete}")
