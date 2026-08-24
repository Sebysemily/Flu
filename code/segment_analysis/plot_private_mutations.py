#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nextclade_ha", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out_png", required=True)
    args = parser.parse_args()

    # 1. Cargar Metadata y definir grupos
    meta = pd.read_csv(args.metadata, sep="\t" if args.metadata.endswith(".tsv") else ",")
    meta.set_index("file_name", inplace=True)
    
    ecuador_ids = meta[meta['country'] == 'Ecuador'].index.tolist()
    context_ids = meta[meta['country'] != 'Ecuador'].index.tolist()
    all_ids = context_ids + ecuador_ids

    # Recolectar todas las mutaciones por muestra
    sample_muts = {sample: [] for sample in all_ids}
    
    nx_path = args.nextclade_ha
    if os.path.exists(nx_path):
        nx_df = pd.read_csv(nx_path, sep=";") 
        if 'seqName' in nx_df.columns:
            nx_df.set_index('seqName', inplace=True)
            
        for sample in all_ids:
            if sample in nx_df.index:
                subs = str(nx_df.loc[sample, "aaSubstitutions"])
                if subs and subs != 'nan':
                    mut_list = subs.split(",")
                    sample_muts[sample].extend([f"HA:{m}" for m in mut_list])

    # 2. Calcular Frecuencias
    mut_counts_ecuador = {}
    mut_counts_context = {}
    
    for sample in ecuador_ids:
        for mut in sample_muts[sample]:
            mut_counts_ecuador[mut] = mut_counts_ecuador.get(mut, 0) + 1
            
    for sample in context_ids:
        for mut in sample_muts[sample]:
            mut_counts_context[mut] = mut_counts_context.get(mut, 0) + 1

    total_ecuador = len([s for s in ecuador_ids if s in nx_df.index]) if os.path.exists(nx_path) else len(ecuador_ids)
    total_context = len([s for s in context_ids if s in nx_df.index]) if os.path.exists(nx_path) else len(context_ids)
    
    # FILTRO: Mutaciones presentes en >=80% de Ecuador (Mutaciones Consenso)
    consensus_mutations = []
    for mut, count in mut_counts_ecuador.items():
        freq_ec = count / total_ecuador
        if freq_ec >= 0.80:
            consensus_mutations.append(mut)

    # Ordenar por posición (asumiendo formato HA:A156S o HA:HA1:A156S)
    consensus_mutations = sorted(consensus_mutations, key=lambda x: int(''.join(filter(str.isdigit, x.split(":")[-1]))))

    if not consensus_mutations:
        print("No se encontraron mutaciones consenso (>80% Ec). Generando imagen vacía...")
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, 'No consensus mutations found\n(>80% Ecuador)', ha='center', va='center', fontsize=12)
        ax.axis('off')
        fig.savefig(args.out_png, bbox_inches='tight')
        return

    # 3. Construir Matriz de Estados para Graficar (Solo Ecuador con HA)
    ecuador_with_ha = [s for s in ecuador_ids if s in nx_df.index] if os.path.exists(nx_path) else ecuador_ids
    state_matrix = pd.DataFrame(0, index=ecuador_with_ha, columns=consensus_mutations)
    val_matrix = pd.DataFrame("-", index=ecuador_with_ha, columns=consensus_mutations)

    for col in consensus_mutations:
        mut = col.split(":")[-1]
        ref_aa = ''.join(filter(str.isalpha, mut))[0]
        alt_aa = ''.join(filter(str.isalpha, mut))[-1]
        
        for sample in ecuador_with_ha:
            if col in sample_muts[sample]:
                state_matrix.loc[sample, col] = 1 # Tiene la mutación privada
                val_matrix.loc[sample, col] = alt_aa
            else:
                state_matrix.loc[sample, col] = 0 # Baseline
                val_matrix.loc[sample, col] = ref_aa

    # 4. GRAFICAR
    n_rows, n_cols = state_matrix.shape
    fig_w = max(11, 0.55 * n_cols + 5.0)
    fig_h = max(9, 0.24 * n_rows + 2.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=220)

    cmap = ListedColormap(['#E8EEF2', '#F28E2B'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
    ax.imshow(state_matrix.values, aspect='auto', cmap=cmap, norm=norm)

    for i in range(n_rows):
        for j in range(n_cols):
            state = state_matrix.iat[i, j]
            aa = val_matrix.iat[i, j]
            color = 'white' if state == 1 else '#1F2933'
            ax.text(j, i, aa, ha='center', va='center', fontsize=6.5, color=color, fontweight='bold')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([c.replace(":", "\n") for c in state_matrix.columns], fontsize=7, rotation=0)
    ax.xaxis.tick_top()
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(state_matrix.index, fontsize=5.5)
    ax.tick_params(length=0)

    ax.set_xticks([x - 0.5 for x in range(1, n_cols)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_rows)], minor=True)
    ax.grid(which='minor', color='white', linewidth=0.7)

    ax.text(n_cols + 0.25, len(ecuador_with_ha) / 2, 'Ecuador Clade', ha='left', va='center', fontsize=8, fontweight='bold')

    legend_handles = [
        mpatches.Patch(color='#F28E2B', label='Ecuador Clade consensus mutation (>80%)')
    ]
    ax.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.17), ncol=1, frameon=False, fontsize=8)

    ax.set_title("Ecuador Clade 1: Consensus HA Mutations Profile", fontsize=13, fontweight='bold', pad=35)

    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches='tight')

if __name__ == '__main__':
    main()
