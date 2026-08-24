#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

# La lista VIP basada en tu requerimiento epidemiológico
CRITICAL_MARKERS = [
    # Adaptación a Mamíferos
    "PB2:E627K", "PB2:D701N", "PB2:Q591K", "PB2:Q591R",
    # Virulencia y Evasión Inmune
    "PB1-F2:N66S", "NS-1:P42S", "HA1-5:T156A", "M1:N30D",
    # Resistencia a Antivirales
    "NA-1:H275Y", "M2:S31N"
]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations_tsv", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out_png", required=True)
    parser.add_argument("--out_list", required=True)
    args = parser.parse_args()

    # 1. Cargar matriz de FluMut y metadata
    df_muts = pd.read_csv(args.mutations_tsv, sep="\t")
    df_muts.set_index("Sample", inplace=True)
    
    meta = pd.read_csv(args.metadata, sep="\t" if args.metadata.endswith(".tsv") else ",")
    meta.set_index("file_name", inplace=True)
    
    valid_samples = [s for s in df_muts.index if s in meta.index]
    
    # 2. Solo graficaremos Ecuador
    ecuador_ids = [s for s in valid_samples if meta.loc[s, 'country'] == 'Ecuador']
    all_ids = ecuador_ids
    
    if not all_ids:
        print("Error: No se encontraron muestras de Ecuador.")
        return
    
    # 3. Filtrar solo los marcadores que existen en el dataframe de FluMut
    available_markers = [m for m in CRITICAL_MARKERS if m in df_muts.columns]
    
    if not available_markers:
        print("Error: No se encontraron los marcadores críticos en el archivo de FluMut.")
        return

    # 4. Crear las matrices para colores y para texto
    state_matrix = pd.DataFrame(0, index=all_ids, columns=available_markers)
    val_matrix = df_muts.loc[all_ids, available_markers].copy()

    for col in available_markers:
        # col format: "PB2:E627K"
        gene, mut = col.split(":")
        ref_aa, alt_aa = mut[0], mut[-1]
        
        for sample in all_ids:
            aa = str(val_matrix.loc[sample, col]).upper().strip()
            if aa == alt_aa:
                state_matrix.loc[sample, col] = 1 # Mutado (Marcador de Alerta)
            elif aa in ['?', '-', 'X', 'NAN', '']:
                state_matrix.loc[sample, col] = 2 # Missing Data
                val_matrix.loc[sample, col] = '?'
            else:
                state_matrix.loc[sample, col] = 0 # Silvestre / Baseline

    # 5. DIBUJAR HEATMAP
    n_rows, n_cols = state_matrix.shape
    fig_w = max(10, 0.65 * n_cols + 4.5)
    fig_h = max(8, 0.24 * n_rows + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=220)

    # Colores: #E8EEF2 (Gris-Aviar), #E11D48 (Rojo-Alerta), #FFFFFF (Blanco-Missing)
    cmap = ListedColormap(['#E8EEF2', '#E11D48', '#FFFFFF'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    ax.imshow(state_matrix.values, aspect='auto', cmap=cmap, norm=norm)

    for i in range(n_rows):
        for j in range(n_cols):
            state = state_matrix.iat[i, j]
            aa = val_matrix.iat[i, j]
            color = 'white' if state == 1 else '#1F2933'
            ax.text(j, i, aa, ha='center', va='center', fontsize=7.5, color=color, fontweight='bold')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([c.replace(":", "\n") for c in state_matrix.columns], fontsize=9, rotation=0)
    ax.xaxis.tick_top()
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(state_matrix.index, fontsize=5.5)
    ax.tick_params(length=0)

    # Cuadrícula
    ax.set_xticks([x - 0.5 for x in range(1, n_cols)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_rows)], minor=True)
    ax.grid(which='minor', color='white', linewidth=0.7)

    # Quitar el divisor de contexto (ya no hay contexto)
    ax.text(n_cols + 0.15, len(ecuador_ids) / 2, 'Ecuador Clade', ha='left', va='center', fontsize=10, fontweight='bold')

    # Leyenda
    legend_handles = [
        mpatches.Patch(color='#E8EEF2', label='Baseline / Avian phenotype'),
        mpatches.Patch(color='#E11D48', label='Critical Marker Detected'),
        mpatches.Patch(facecolor='#FFFFFF', edgecolor='#808080', label='Missing data (partial seq)')
    ]
    ax.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False, fontsize=9)

    ax.set_title("Key Molecular Signatures for Mammalian Adaptation and Virulence", fontsize=15, fontweight='bold', pad=40)

    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches='tight')

    # 6. Escribir la lista de marcadores al archivo de salida
    with open(args.out_list, "w") as f:
        f.write("Critical Molecular Signatures Assessed:\n")
        f.write("=======================================\n")
        for marker in CRITICAL_MARKERS:
            f.write(f"- {marker}\n")

if __name__ == '__main__':
    main()
