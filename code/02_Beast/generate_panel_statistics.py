#!/usr/bin/env python3
import os
import sys
import re
import csv
import pandas as pd

def parse_host_from_id(id_str):
    id_lower = id_str.lower()
    if 'chicken' in id_lower or 'gallus' in id_lower:
        return 'Chicken'
    elif 'turkey' in id_lower:
        return 'Turkey'
    elif 'duck' in id_lower or 'pato' in id_lower or 'wigeon' in id_lower or 'anas' in id_lower:
        return 'Duck'
    elif 'pelican' in id_lower or 'pelecanus' in id_lower:
        return 'Pelican'
    elif 'cormorant' in id_lower or 'guanay' in id_lower or 'phalacrocorax' in id_lower:
        return 'Cormorant'
    elif 'gull' in id_lower or 'larus' in id_lower:
        return 'Gull'
    elif 'tern' in id_lower or 'thalasseus' in id_lower or 'sternula' in id_lower:
        return 'Tern'
    elif 'booby' in id_lower or 'piquero' in id_lower or 'sula' in id_lower:
        return 'Booby'
    elif 'swan' in id_lower:
        return 'Swan'
    elif 'vulture' in id_lower:
        return 'Vulture'
    elif 'goose' in id_lower or 'anser' in id_lower or 'chloephaga' in id_lower:
        return 'Goose'
    elif 'sanderling' in id_lower or 'calidris' in id_lower:
        return 'Sanderling'
    elif 'skimmer' in id_lower or 'black_skimmer' in id_lower:
        return 'Skimmer'
    elif 'oystercatcher' in id_lower or 'haematopus' in id_lower:
        return 'Oystercatcher'
    elif 'grebe' in id_lower or 'grabe' in id_lower:
        return 'Grebe'
    elif 'penguin' in id_lower or 'spheniscus' in id_lower:
        return 'Penguin'
    elif 'egret' in id_lower or 'ardea' in id_lower:
        return 'Egret'
    elif 'falcon' in id_lower or 'falco' in id_lower:
        return 'Falcon'
    elif 'owl' in id_lower or 'strix' in id_lower or 'bubo' in id_lower:
        return 'Owl'
    elif any(k in id_lower for k in ['sea_lion', 'sealion', 'sea-lion', 'lobos', 'otaria', 'phocarctos', 'fur_seal', 'arctocephalus']):
        return 'Sea Lion (Mammal)'
    elif 'human' in id_lower:
        return 'Human'
    elif 'avian' in id_lower:
        return 'Avian (unspecified)'
    elif 'wild' in id_lower:
        return 'Wild Bird'
    else:
        # Try to extract capitalized name
        match = re.match(r'^A([A-Z][a-z_]+)', id_str)
        if match:
            return match.group(1).replace('_', ' ')
        return 'Other'

def map_place_to_country(place, source_type):
    if source_type in ['USFQ Core', 'Galapagos Core (MAATE)']:
        return 'Ecuador'
    
    place_clean = place.strip()
    if place_clean in ['Argentina', 'Uruguay', 'Bolivia', 'Peru', 'Chile', 'Ecuador', 'Brazil', 'Colombia', 'Venezuela', 'Paraguay']:
        return place_clean
    
    # Peru provinces
    if place_clean in ['Ica', 'Lima', 'Piura', 'Lambayeque', 'Ancash', 'Arequipa', 'LaLibertad', 'Tumbes']:
        return 'Peru'
        
    # Chile provinces
    if place_clean in ['Antofagasta', 'Araucania', 'Arica', 'AricaYParinacota', 'Arica_y_Parinacota', 'Atacama', 
                       'Aysen', 'BioBio', 'Biobio', 'Coquimbo', 'Maule', 'Metropolitana', 'Nuble', 'Ohiggins', 
                       'OHiggins', 'Tarapaca', 'Valparaiso', 'LosLagos', 'LosRios', 'Santiago']:
        return 'Chile'
        
    # Bolivia provinces
    if place_clean in ['Cochabamba', 'Potosi']:
        return 'Bolivia'
        
    # Argentina provinces
    if place_clean in ['Neuquen', 'BuenosAires', 'Chubut', 'RioNegro', 'SantaCruz', 'TierraDelFuego']:
        return 'Argentina'
        
    return place_clean

def map_local_species(specie_str):
    if not isinstance(specie_str, str):
        return 'Other/Unknown'
    specie_lower = specie_str.lower()
    if 'corral' in specie_lower:
        return 'Chicken'
    elif 'fragata' in specie_lower:
        return 'Frigatebird'
    elif 'piquero' in specie_lower:
        return 'Booby'
    elif 'pardela' in specie_lower or 'puffinus' in specie_lower or 'ardenna' in specie_lower:
        return 'Shearwater'
    elif 'elegans' in specie_lower or 'tern' in specie_lower:
        return 'Tern'
    elif 'patos' in specie_lower or 'pato' in specie_lower:
        return 'Duck'
    elif 'pelicano' in specie_lower:
        return 'Pelican'
    elif 'gaviota' in specie_lower:
        return 'Gull'
    return specie_str

def main():
    panel_tsv = sys.argv[1] if len(sys.argv) > 1 else 'data/beast/panel_main_taxa.final.tsv'
    output_csv = sys.argv[2] if len(sys.argv) > 2 else 'results/qc_metrics/phylogeny/main_analysis_panel.csv'
    
    default_metadata = 'config/flu_filtrado.csv'
    for filename in ['config/paths.yml', 'config/config.yml']:
        if os.path.exists(filename):
            try:
                import yaml
                with open(filename) as fh:
                    cfg = yaml.safe_load(fh) or {}
                    if 'flu_filtrado' in cfg:
                        default_metadata = cfg['flu_filtrado']
            except Exception:
                pass
                
    metadata_csv = sys.argv[3] if len(sys.argv) > 3 else default_metadata
    
    if not os.path.exists(panel_tsv):
        print(f"Error: Panel file not found at {panel_tsv}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Reading panel from {panel_tsv}")
    df_panel = pd.read_csv(panel_tsv, sep='\t')
    
    # Load metadata lookup if available
    metadata_lookup = {}
    if os.path.exists(metadata_csv):
        print(f"Loading local metadata from {metadata_csv}")
        df_meta = pd.read_csv(metadata_csv)
        for _, row in df_meta.iterrows():
            code = str(row.get('Código USFQ', '')).strip()
            if code:
                metadata_lookup[code] = {
                    'species': row.get('Especie', 'Unknown'),
                    'accession': row.get('EPI_ISL', '')
                }
                
    maate_hosts = {
        'EPI_ISL_17973443': ('Fregata magnificens (Frigatebird)', 'Frigatebird'),
        'EPI_ISL_17973458': ('Fregata magnificens (Frigatebird)', 'Frigatebird'),
        'EPI_ISL_18137671': ('Sula nebouxii (Blue-footed booby)', 'Booby'),
    }
    
    parsed_rows = []
    
    for _, row in df_panel.iterrows():
        taxon = row['taxon']
        role = row['role']
        lineage = row.get('lineage', '')
        dist = row.get('distance_to_seed', '')
        
        parts = taxon.split('/')
        if len(parts) < 3:
            print(f"Warning: Unexpected taxon format: {taxon}", file=sys.stderr)
            continue
            
        id_part = parts[0]
        place = parts[1]
        date_str = parts[2]
        
        # Parse year/month
        year = date_str[:4] if date_str and len(date_str) >= 4 else 'UNKNOWN'
        month = date_str[:7] if date_str and len(date_str) >= 7 else 'UNKNOWN'
        
        # Determine source, host and accession
        accession = ''
        source_type = 'GISAID Context'
        host_raw = 'Unknown'
        host_clean = 'Unknown'
        
        if id_part.startswith('Flu-EPI_ISL_'):
            # MAATE Galapagos core
            source_type = 'Galapagos Core (MAATE)'
            accession = id_part.replace('Flu-', '')
            if accession in maate_hosts:
                host_raw, host_clean = maate_hosts[accession]
            else:
                host_raw = 'Wild Bird'
                host_clean = 'Wild Bird'
        elif id_part.startswith('Flu-'):
            # USFQ Core
            source_type = 'USFQ Core'
            if id_part in metadata_lookup:
                host_raw = metadata_lookup[id_part]['species']
                host_clean = map_local_species(host_raw)
                accession = metadata_lookup[id_part]['accession']
            else:
                host_raw = 'Chicken'
                host_clean = 'Chicken'
        else:
            # GISAID Context
            source_type = 'GISAID Context'
            host_clean = parse_host_from_id(id_part)
            host_raw = id_part.split('_')[0] if '_' in id_part else id_part
            
            # Extract accession
            match = re.search(r'EPI_ISL_\d+', id_part)
            if match:
                accession = match.group(0)
                
        country = map_place_to_country(place, source_type)
        
        parsed_rows.append({
            'taxon': taxon,
            'role': role,
            'source_type': source_type,
            'accession': accession,
            'place': place,
            'country': country,
            'date': date_str,
            'year': year,
            'month': month,
            'host_raw': host_raw,
            'host_clean': host_clean,
            'lineage': lineage,
            'distance_to_seed': dist
        })
        
    df_out = pd.DataFrame(parsed_rows)
    
    # Save CSV
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_out.to_csv(output_csv, index=False)
    print(f"Wrote main analysis panel CSV to {output_csv}")
    
    # Generate a stdout summary report
    print("\n" + "="*50)
    print("             MAIN ANALYSIS PANEL SUMMARY STATISTICS")
    print("="*50)
    print(f"Total Taxa in Panel: {len(df_out)}")
    print("\n--- Composition by Role ---")
    print(df_out['role'].value_counts().to_string())
    
    print("\n--- Composition by Source Type ---")
    print(df_out['source_type'].value_counts().to_string())
    
    print("\n--- Composition by Country ---")
    print(df_out['country'].value_counts().to_string())
    
    print("\n--- Composition by Host (Cleaned) ---")
    print(df_out['host_clean'].value_counts().to_string())
    
    print("\n--- Composition by Year ---")
    print(df_out['year'].value_counts().to_string())
    
    # Also write a summary text file
    summary_txt = output_csv.replace('.csv', '_summary.txt')
    with open(summary_txt, 'w') as f:
        f.write("MAIN ANALYSIS PANEL SUMMARY STATISTICS\n")
        f.write("======================================\n")
        f.write(f"Total Taxa in Panel: {len(df_out)}\n\n")
        f.write("--- Composition by Role ---\n")
        f.write(df_out['role'].value_counts().to_string() + "\n\n")
        f.write("--- Composition by Source Type ---\n")
        f.write(df_out['source_type'].value_counts().to_string() + "\n\n")
        f.write("--- Composition by Country ---\n")
        f.write(df_out['country'].value_counts().to_string() + "\n\n")
        f.write("--- Composition by Host (Cleaned) ---\n")
        f.write(df_out['host_clean'].value_counts().to_string() + "\n\n")
        f.write("--- Composition by Year ---\n")
        f.write(df_out['year'].value_counts().to_string() + "\n")
        
    print(f"Wrote summary text to {summary_txt}")

if __name__ == '__main__':
    main()
