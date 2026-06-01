import argparse
import sys
import numpy as np
from collections import Counter

def read_fasta(path):
    records = []
    header = None
    chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records

def write_fasta(path, records):
    with open(path, "w") as out:
        for h, s in records:
            out.write(f">{h}\n{s}\n")

def get_consensus_and_trim_limits(records):
    # Convert to matrix
    seqs = [list(s.upper()) for h, s in records]
    mat = np.array(seqs)
    
    n_seqs, aln_len = mat.shape
    consensus_chars = []
    
    for i in range(aln_len):
        col = mat[:, i]
        # Most common character, ignoring gaps if possible, but if gap is overwhelming, keep it?
        # Actually, just count all
        counts = Counter(col)
        # remove Ns and gaps if we just want to find the consensus nucleotide
        # but if it's an insertion in only 1 sequence, gap will be majority.
        most_common = counts.most_common(1)[0][0]
        consensus_chars.append(most_common)
        
    cons_array = np.array(consensus_chars)
    
    # Map from ungapped consensus to alignment column
    ungapped_cons = []
    col_map = []
    for i, c in enumerate(cons_array):
        if c not in ('-', 'N'):
            ungapped_cons.append(c)
            col_map.append(i)
            
    ungapped_str = "".join(ungapped_cons)
    
    # Find longest ORF in ungapped consensus
    best_start = -1
    best_end = -1
    max_len = -1
    
    # Simple longest ORF finder (forward strand)
    for i in range(len(ungapped_str) - 2):
        if ungapped_str[i:i+3] == "ATG":
            # find stop
            for j in range(i+3, len(ungapped_str) - 2, 3):
                codon = ungapped_str[j:j+3]
                if codon in ("TAA", "TAG", "TGA"):
                    orf_len = j + 3 - i
                    if orf_len > max_len:
                        max_len = orf_len
                        best_start = i
                        best_end = j + 3
                    break
                    
    if best_start == -1:
        # Fallback if no ORF found
        print("Warning: No clear ORF found. Falling back to whole alignment.")
        return 0, aln_len, cons_array
        
    start_col = col_map[best_start]
    end_col = col_map[best_end - 1] + 1
    
    return start_col, end_col, cons_array

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-divergence", type=float, default=0.10, help="Max fraction of gaps or Ns (ambiguities) in CDS")
    args = parser.parse_args()
    
    records = read_fasta(args.input)
    if not records:
        print("Empty input.")
        write_fasta(args.output, [])
        return
        
    start_col, end_col, consensus = get_consensus_and_trim_limits(records)
    print(f"Trimming alignment from column {start_col} to {end_col}")
    
    trimmed_cons = consensus[start_col:end_col]
    
    filtered_records = []
    dropped = 0
    
    for h, s in records:
        trimmed_seq = s[start_col:end_col].upper()
        
        # Calculate fraction of gaps / Ns in the trimmed region
        diffs = sum(1 for sc in trimmed_seq if sc == '-' or sc == 'N')
        valid_len = len(trimmed_cons)
        
        div = diffs / valid_len if valid_len > 0 else 0
        
        if div > args.max_divergence:
            dropped += 1
        else:
            filtered_records.append((h, trimmed_seq))
            
    print(f"Dropped {dropped} sequences due to >{args.max_divergence*100}% gaps or Ns.")
    write_fasta(args.output, filtered_records)

if __name__ == "__main__":
    main()
