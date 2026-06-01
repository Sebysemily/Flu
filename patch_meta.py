import pandas as pd
import shutil

df = pd.read_csv("metadata/context_base.csv")

patch_data = {
    "EPI_ISL_18137671": "Ecuador / Provincia del Guayas",
    "EPI_ISL_18137626": "Ecuador / Provincia del Guayas",
    "EPI_ISL_17973458": "Ecuador / Provincia de Manabi",
    "EPI_ISL_17973443": "Ecuador / Provincia de Manabi",
    "EPI_ISL_16161675": "Ecuador / Provincia de Cotopaxi",
    "EPI_ISL_16161673": "Ecuador / Provincia de Cotopaxi"
}

updated = 0
for idx, row in df.iterrows():
    epi = str(row["file_name"]).strip()
    if epi in patch_data:
        df.at[idx, "country"] = patch_data[epi]
        # Also assign expected_role if needed, but let's just do country as requested
        if "Guayas" in patch_data[epi] or "Manabi" in patch_data[epi]:
            df.at[idx, "expected_role"] = "flu_costa"
        elif "Cotopaxi" in patch_data[epi]:
            df.at[idx, "expected_role"] = "flu_sierra"
        updated += 1

# Check if the other 3 are even in the file (they were from the extra fastas)
for epi, loc in patch_data.items():
    if not (df["file_name"] == epi).any():
        # Add them
        role = "flu_costa" if "Guayas" in loc or "Manabi" in loc else "flu_sierra"
        new_row = {"file_name": epi, "collection_date": "2023-01-01", "country": loc, "expected_role": role}
        df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
        updated += 1
        print(f"Added new row for {epi}")

print(f"Updated/Added {updated} records in context_base.csv")
df.to_csv("metadata/context_base.csv", index=False)
