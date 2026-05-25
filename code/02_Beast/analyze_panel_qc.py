#!/usr/bin/env python3
import csv
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def parse_fasta(path):
    records = {}
    for record in SeqIO.parse(path, "fasta"):
        records[record.id] = str(record.seq)
    return records

def get_gap_stats(seq):
    # Count leading, trailing, and internal gaps
    stripped = seq.strip("-")
    if not stripped:
        return len(seq), 0, 0, ""
    
    leading_gaps = len(seq) - len(seq.lstrip("-"))
    trailing_gaps = len(seq) - len(seq.rstrip("-"))
    internal_gaps = stripped.count("-")
    
    return leading_gaps, trailing_gaps, internal_gaps, stripped

def analyze_translation(stripped_seq):
    # Remove all gaps to translate the raw nucleotide sequence
    raw_nucs = stripped_seq.replace("-", "")
    
    # Check if length is multiple of 3
    is_multiple_of_3 = (len(raw_nucs) % 3 == 0)
    
    # Translate (using standard genetic code)
    try:
        translation = str(Seq(raw_nucs).translate(to_stop=False))
    except Exception:
        translation = ""
        
    if not translation:
        return 0, False, False, 0
    
    # Internal stops: stop codons (*) before the very last position
    internal_stops = translation[:-1].count("*")
    has_stop_at_end = (translation[-1] == "*")
    total_stops = translation.count("*")
    
    return len(raw_nucs), is_multiple_of_3, has_stop_at_end, internal_stops

def main():
    panel_tsv = "data/pre_beast/panels/panel_main_taxa.filtered.tsv"
    processed_dir = "data/processed_alignments/codon_aware"
    out_csv = "results/qc_metrics/beast_panel_qc_report.csv"
    
    if not os.path.exists(panel_tsv):
        print(f"Error: Panel taxa file not found at {panel_tsv}")
        sys.exit(1)
        
    print(f"Loading BEAST panel taxa from {panel_tsv}")
    df_panel = pd.read_csv(panel_tsv, sep="\t")
    taxa_list = df_panel["taxon"].tolist()
    taxon_to_role = dict(zip(df_panel["taxon"], df_panel["role"]))
    
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    
    # Load all processed alignments
    alignments = {}
    for seg in segments:
        aln_path = os.path.join(processed_dir, f"H5N1_{seg}.mafft")
        if os.path.exists(aln_path):
            print(f"Loading alignment for segment {seg}...")
            alignments[seg] = parse_fasta(aln_path)
        else:
            print(f"Warning: Alignment file for segment {seg} not found at {aln_path}")
            alignments[seg] = {}
            
    rows = []
    
    for taxon in taxa_list:
        role = taxon_to_role[taxon]
        
        for seg in segments:
            seq = alignments[seg].get(taxon)
            
            if seq is None:
                # Taxon is missing this segment
                rows.append({
                    "taxon": taxon,
                    "role": role,
                    "segment": seg,
                    "present": "NO",
                    "aligned_length": 0,
                    "nuc_length": 0,
                    "leading_gaps": 0,
                    "trailing_gaps": 0,
                    "internal_gaps": 0,
                    "starts_at_cp1": "N/A",
                    "frame_preserved": "N/A",
                    "stop_at_end": "N/A",
                    "internal_stops": 0,
                    "status": "MISSING"
                })
                continue
                
            leading, trailing, internal, stripped = get_gap_stats(seq)
            nuc_len, is_mult3, stop_end, int_stops = analyze_translation(stripped)
            
            starts_cp1 = "YES" if leading == 0 else "NO"
            frame_preserved = "YES" if is_mult3 else "NO"
            stop_at_end = "YES" if stop_end else "NO"
            
            # Status determination
            if int_stops > 0 or not is_mult3:
                status = "ORF_CHANGED"
            elif internal > 0:
                status = "GAPPY_OK"
            else:
                status = "CLEAN"
                
            rows.append({
                "taxon": taxon,
                "role": role,
                "segment": seg,
                "present": "YES",
                "aligned_length": len(seq),
                "nuc_length": nuc_len,
                "leading_gaps": leading,
                "trailing_gaps": trailing,
                "internal_gaps": internal,
                "starts_at_cp1": starts_cp1,
                "frame_preserved": frame_preserved,
                "stop_at_end": stop_end,
                "internal_stops": int_stops,
                "status": status
            })
            
    # Write output CSV
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)
    print(f"Wrote detailed QC report to {out_csv}")
    
    # Print summary statistics
    print("\n" + "="*50)
    print("QC SUMMARY FOR 200 CONTEXT + 32 CORE TAXA")
    print("="*50)
    print(f"Total Taxa Analyzed: {len(taxa_list)}")
    print(f"Total Segments Checked: {len(rows)}")
    
    print("\n1. ORF status counts:")
    print(df_out["status"].value_counts().to_string())
    
    print("\n2. Position CP1 Start counts:")
    print(df_out["starts_at_cp1"].value_counts().to_string())
    
    print("\n3. Frame Preservation (length multiple of 3):")
    print(df_out["frame_preserved"].value_counts().to_string())
    
    print("\n4. Internal gaps (blank spaces) count distribution (excluding 0):")
    gappy = df_out[df_out["internal_gaps"] > 0]
    if not gappy.empty:
        print(gappy["internal_gaps"].describe().to_string())
    else:
        print("No internal gaps found in any sequence.")
        
    print("\n5. Details on ORF_CHANGED (internal stop codons or frameshifts) by segment:")
    orf_chg = df_out[df_out["status"] == "ORF_CHANGED"]
    if not orf_chg.empty:
        print(orf_chg.groupby("segment").size().to_string())
        print("\nExamples of ORF_CHANGED sequences:")
        print(orf_chg[["taxon", "segment", "nuc_length", "frame_preserved", "internal_stops"]].head(10).to_string())
    else:
        print("No ORF_CHANGED sequences found! (All sequences have intact open reading frames)")

if __name__ == "__main__":
    main()
