import os
import re
import glob
from collections import defaultdict, Counter
import pandas as pd
import yaml
import shutil

CONFIG_FILE = "config/config.yml"
OUTDIR = os.path.join("data", "assembled")
AMENDED_DIR = os.path.join("data", "all_amended_fasta")

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

with open(CONFIG_FILE) as fh:
    config = yaml.safe_load(fh)

FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")

MIRA_BASE = config.get("mira_base_dir", "..")

# Load valid samples from flu_filtrado.csv
filtrado_df = pd.read_csv(FILTRADO_CSV)
valid_samples = set(filtrado_df["Código USFQ"].unique())
print(f"Muestras cargadas de {FILTRADO_CSV}: {len(valid_samples)}")

def norm_sample(x):
    x = str(x).strip()
    m = re.match(r"^(Flu-\d+)", x)
    return m.group(1) if m else x

def clean_seq(seq):
    return seq.upper().replace(" ", "").replace("\n", "").replace("-", "")

def parse_fasta(path):
    records = []
    header = None
    seq_chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records

def mira_header_to_sample_seg(header):
    if "|" not in header:
        return None, None

    sample_part, seg_part = header.split("|", 1)
    sample = norm_sample(sample_part)

    parts = seg_part.split("_")
    if len(parts) < 2:
        return sample, None

    seg = parts[1].upper()
    seg_map = {
        "PB2": "PB2",
        "PB1": "PB1",
        "PA": "PA",
        "HA": "HA",
        "NP": "NP",
        "NA": "NA",
        "MP": "MP",
        "NS": "NS",
    }

    return sample, seg_map.get(seg)

def best_expected_length(lengths):
    if not lengths:
        return None

    counts = Counter(lengths)
    top_n = max(counts.values())
    modes = [k for k, v in counts.items() if v == top_n]

    if len(modes) == 1:
        return modes[0]

    return max(lengths)

def wrap_seq(seq, width=80):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

# Create output directories automatically
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(os.path.join(OUTDIR, "ecuador_intermediate_per_sample"), exist_ok=True)
print(f"Output directory ready: {OUTDIR}")

mira_fastas = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
mira_fastas += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
mira_fastas = sorted(set(mira_fastas))

print("FASTAs detectados:")
for f in mira_fastas:
    print(" -", f)

# Create amended FASTA directory and copy files
os.makedirs(AMENDED_DIR, exist_ok=True)
print(f"Copiando FASTA a {AMENDED_DIR}...")

copied_fastas = []
for fasta in mira_fastas:
    run_name = os.path.basename(os.path.dirname(fasta))
    dest = os.path.join(AMENDED_DIR, f"{run_name}_amended_consensus.fasta")
    shutil.copy2(fasta, dest)
    copied_fastas.append(dest)
    print(f"  {fasta} -> {dest}")

# Use copied FASTA for processing
mira_fastas = copied_fastas

best_records = defaultdict(dict)
sample_runs = defaultdict(set)
lengths_by_seg = defaultdict(list)
raw_rows = []

for fasta in mira_fastas:
    run = os.path.basename(os.path.dirname(fasta))

    for header, seq in parse_fasta(fasta):
        sample, seg = mira_header_to_sample_seg(header)
        if sample is None or seg not in SEGMENTS:
            continue

        seq = clean_seq(seq)
        non_n = sum(1 for b in seq if b in "ACGT")

        raw_rows.append({
            "sample_norm": sample,
            "run": run,
            "segment": seg,
            "header": header,
            "seq_len": len(seq),
            "non_n_bases": non_n,
        })

        if non_n == 0:
            continue

        sample_runs[sample].add(run)
        lengths_by_seg[seg].append(non_n)

        prev = best_records[sample].get(seg)
        if prev is None or non_n > prev["non_n_bases"]:
            best_records[sample][seg] = {
                "run": run,
                "header": header,
                "sequence": seq,
                "non_n_bases": non_n
            }

raw_df = pd.DataFrame(raw_rows)
raw_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_raw_segments.csv"), index=False)

expected_lengths = {}
for seg in SEGMENTS:
    exp = best_expected_length(lengths_by_seg[seg])
    if exp is None:
        raise SystemExit(f"No pude estimar largo esperado para {seg}")
    expected_lengths[seg] = exp

expected_df = pd.DataFrame([
    {"segment": seg, "expected_length": expected_lengths[seg]}
    for seg in SEGMENTS
])
expected_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_expected_lengths.csv"), index=False)

summary_rows = []
fasta_audit_rows = []

# Track all found samples and which ones were processed
all_found_samples = set(best_records.keys())
processed_samples = set()
not_in_filtrado_samples = []
discrepancy_rows = []

for sample in sorted(best_records.keys()):
    # Filter: only process samples in flu_filtrado.csv
    if sample not in valid_samples:
        not_in_filtrado_samples.append(sample)
        continue
    
    processed_samples.add(sample)
    segs_found = set(best_records[sample].keys())
    runs = sorted(sample_runs[sample])

    row = {
        "Código USFQ": sample,
        "run": ",".join(runs),
    }

    n_regions = 0
    out_fa = os.path.join(OUTDIR, "ecuador_intermediate_per_sample", f"{sample}.fasta")

    # Get expected segments from flu_filtrado.csv for this sample
    sample_filtrado = filtrado_df[filtrado_df["Código USFQ"] == sample].iloc[0]
    expected_segs_in_filtrado = {seg for seg in SEGMENTS if str(sample_filtrado.get(seg, "")).strip().upper() == "SI"}

    with open(out_fa, "w") as out:
        for seg in SEGMENTS:
            if seg in segs_found:
                row[seg] = "SI"
                n_regions += 1
                seq = best_records[sample][seg]["sequence"]
                source_run = best_records[sample][seg]["run"]
                status = "assembled"
            else:
                row[seg] = ""
                seq = "N" * expected_lengths[seg]
                source_run = ""
                status = "filled_with_N"

            out.write(f">{sample}|A_{seg}\n")
            out.write(wrap_seq(seq) + "\n")

            fasta_audit_rows.append({
                "Código USFQ": sample,
                "segment": seg,
                "status": status,
                "source_run": source_run,
                "sequence_length_written": len(seq)
            })

    row["n_regions_assembled"] = n_regions
    row["has_8_regions"] = "SI" if n_regions == 8 else ""
    summary_rows.append(row)
    
    # Check for discrepancies: segments expected in filtrado but not found in MIRA
    missing_in_mira = expected_segs_in_filtrado - segs_found
    if missing_in_mira:
        discrepancy_rows.append({
            "Código USFQ": sample,
            "reason": f"Segments marked SI in flu_filtrado but NOT found in MIRA: {','.join(sorted(missing_in_mira))}",
            "found_in_mira": "PARTIAL",
            "missing_segments": ",".join(sorted(missing_in_mira)),
            "filled_with_Ns": "YES"
        })

summary_df = pd.DataFrame(summary_rows)
summary_df = summary_df[
    ["Código USFQ", "run"] + SEGMENTS + ["n_regions_assembled", "has_8_regions"]
]
summary_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_summary.csv"), index=False)

audit_df = pd.DataFrame(fasta_audit_rows)
audit_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_audit.csv"), index=False)

audit_df = pd.DataFrame(fasta_audit_rows)
audit_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_audit.csv"), index=False)

# Generate report of samples NOT processed and discrepancies
not_processed_rows = []

# Add samples found in MIRA but not in flu_filtrado.csv
for sample in not_in_filtrado_samples:
    not_processed_rows.append({
        "Código USFQ": sample,
        "status": "SKIPPED",
        "reason": "Not in config/flu_filtrado.csv",
        "found_in_mira": "YES",
        "missing_segments": "",
        "filled_with_Ns": ""
    })

# Add samples in flu_filtrado.csv but not found in MIRA
for sample in sorted(valid_samples):
    if sample not in all_found_samples:
        not_processed_rows.append({
            "Código USFQ": sample,
            "status": "NOT_FOUND",
            "reason": "Not found in MIRA amended_consensus.fasta files",
            "found_in_mira": "NO",
            "missing_segments": "ALL",
            "filled_with_Ns": "NO"
        })

# Add discrepancies (processed but with missing expected segments)
for disc in discrepancy_rows:
    not_processed_rows.append({
        "Código USFQ": disc["Código USFQ"],
        "status": "PROCESSED_WITH_ISSUES",
        "reason": disc["reason"],
        "found_in_mira": disc["found_in_mira"],
        "missing_segments": disc["missing_segments"],
        "filled_with_Ns": disc["filled_with_Ns"]
    })

not_processed_df = pd.DataFrame(not_processed_rows)
not_processed_df.to_csv(os.path.join(OUTDIR, "ecuador_intermediate_issues.csv"), index=False)

combined_fa = os.path.join(OUTDIR, "ecuador_intermediate_sequences.fasta")
with open(combined_fa, "w") as out_all:
    for fa in sorted(glob.glob(os.path.join(OUTDIR, "ecuador_intermediate_per_sample", "*.fasta"))):
        with open(fa) as fh:
            out_all.write(fh.read().rstrip("\n") + "\n")

print("Listo.")
print("Archivos generados:")
print(" - data/all_amended_fasta/  [COPIA DE INPUTS]")
print(" - data/assembled/ecuador_intermediate_raw_segments.csv")
print(" - data/assembled/ecuador_intermediate_expected_lengths.csv")
print(" - data/assembled/ecuador_intermediate_summary.csv")
print(" - data/assembled/ecuador_intermediate_audit.csv")
print(" - data/assembled/ecuador_intermediate_issues.csv  [VALIDACIONES]")
print(" - data/assembled/ecuador_intermediate_sequences.fasta")
print(" - data/assembled/ecuador_intermediate_per_sample/")
print()
print("=" * 60)
print("RESUMEN DE VALIDACIÓN:")
print("=" * 60)
print(f"✓ Muestras procesadas exitosamente: {len(summary_df) - len(discrepancy_rows)}")
print(f"✓ Con 8 regiones ensambladas: {(summary_df['has_8_regions'] == 'SI').sum()}")
print(f"⚠ Muestras procesadas CON DISCREPANCIAS: {len(discrepancy_rows)}")
print(f"✗ Muestras NO encontradas en MIRA: {len([r for r in not_processed_rows if r['status'] == 'NOT_FOUND'])}")
print(f"✗ Muestras DESCARTADAS (no en flu_filtrado): {len([r for r in not_processed_rows if r['status'] == 'SKIPPED'])}")
print(f"📊 Ver detalles en: data/assembled/ecuador_intermediate_issues.csv")
