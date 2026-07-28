#!/usr/bin/env python3
import sys
import os
import argparse

def deduplicate_robust(input_file, output_file):
    seen_segments = set()
    duplicates_removed = 0
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = None
        seq = []
        
        def write_seq():
            nonlocal duplicates_removed
            if header is not None:
                # Extract Sample Name
                parts = header.split("|")
                if len(parts) > 1 and parts[0].startswith(">"):
                    name = parts[0].split(">")[1].strip()
                else:
                    name = header.strip()
                
                # Extract segment
                seg = "UNKNOWN"
                if "|PB2" in header or header.endswith("PB2"): seg = "PB2"
                elif "|PB1" in header or header.endswith("PB1"): seg = "PB1"
                elif "|PA" in header or header.endswith("PA"): seg = "PA"
                elif "|HA" in header or header.endswith("HA"): seg = "HA"
                elif "|NP" in header or header.endswith("NP"): seg = "NP"
                elif "|NA" in header or header.endswith("NA"): seg = "NA"
                elif "|MP" in header or header.endswith("MP"): seg = "MP"
                elif "|NS" in header or header.endswith("NS"): seg = "NS"
                
                identifier = f"{name}_{seg}"
                
                if identifier in seen_segments:
                    duplicates_removed += 1
                else:
                    seen_segments.add(identifier)
                    f_out.write(header + "\n")
                    for line in seq:
                        f_out.write(line + "\n")
                        
        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                write_seq()
                header = line
                seq = []
            else:
                seq.append(line)
        write_seq()
        
    print(f"[{input_file}] Removed {duplicates_removed} duplicate segments based on sample name and segment type.")

def main():
    parser = argparse.ArgumentParser(description="Deduplicate GISAID FASTA files by isolate name and segment.")
    parser.add_argument("input_fasta", help="Path to input FASTA file")
    parser.add_argument("-o", "--output", help="Path to output FASTA (modifies in place if not provided)", default=None)
    
    args = parser.parse_args()
    
    in_file = args.input_fasta
    out_file = args.output if args.output else in_file + ".tmp"
    
    deduplicate_robust(in_file, out_file)
    
    if not args.output:
        os.rename(out_file, in_file)

if __name__ == "__main__":
    main()
