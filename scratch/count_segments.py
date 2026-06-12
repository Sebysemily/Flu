import csv

def process_file(filename):
    with open(filename, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
        # Check if segment columns exist
        if not set(segments).issubset(set(reader.fieldnames)):
            print(f"Segments columns not found in {filename}")
            return
            
        complete = 0
        ha_only = 0
        others = 0
        total = 0
        
        for row in reader:
            if row.get("expected_role") == "regional_context":
                total += 1
                present_count = 0
                has_ha = False
                for seg in segments:
                    val = row.get(seg, "").strip()
                    if val and val.lower() != "nan":
                        present_count += 1
                        if seg == "HA":
                            has_ha = True
                
                if present_count == 8:
                    complete += 1
                elif present_count == 1 and has_ha:
                    ha_only += 1
                else:
                    others += 1
                    
        print(f"\nUsing {filename}")
        print(f"Total regional_context: {total}")
        print(f"Complete genomes (8 segments): {complete}")
        print(f"HA region only (1 segment, HA): {ha_only}")
        print(f"Other combinations: {others}")

process_file("metadata/H5N1_context.csv")
process_file("metadata/context_base.csv")
