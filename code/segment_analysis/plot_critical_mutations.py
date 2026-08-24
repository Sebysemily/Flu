#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.lines as mlines

CRITICAL_MARKERS = [
    "PB2:E627K", "PB2:D701N", "PB2:Q591K", "PB2:Q591R",
    "PA:N383D", "PB1-F2:N66S", "HA1-5:S223R", "HA1-5:T156A",
    "NP:A184K", "NS-1:P42S", "M1:N30D", "NA-1:H275Y", "M2:S31N"
]

GEO_COLOR_MAP = {
    "flu_costa": "#FF0000",   
    "flu_sierra": "#00008B",
    "flu_amazonia": "#008000",
    "Peru": "#A34700",
    "Chile": "#1E88E5",
    "default_accent": "#E11D48"
}

HOST_COLOR_MAP = {
    "domesticated bird": "#E6AB02", 
    "wild bird": "#7570B3",         
    "mammal": "#00BCD4",            
    "unknown": "#999999"            
}

def get_geo_color(meta_row):
    country = meta_row.get('country', '')
    role = meta_row.get('expected_role', '')
    if country == 'Ecuador':
        if role in GEO_COLOR_MAP:
            return GEO_COLOR_MAP[role]
        return GEO_COLOR_MAP["default_accent"]
    elif country in GEO_COLOR_MAP:
        return GEO_COLOR_MAP[country]
    else:
        return GEO_COLOR_MAP["default_accent"]

def get_host_color(meta_row):
    htype = str(meta_row.get('host_type', '')).lower()
    if 'domesticated bird' in htype:
        return HOST_COLOR_MAP['domesticated bird']
    elif 'wild bird' in htype:
        return HOST_COLOR_MAP['wild bird']
    elif 'mammal' in htype:
        return HOST_COLOR_MAP['mammal']
    else:
        return HOST_COLOR_MAP['unknown']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations_tsv", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out_png", required=True)
    parser.add_argument("--out_list", required=True)
    args = parser.parse_args()

    df_muts = pd.read_csv(args.mutations_tsv, sep="\t")
    df_muts.set_index("Sample", inplace=True)
    
    meta = pd.read_csv(args.metadata, sep="\t" if args.metadata.endswith(".tsv") else ",")
    meta.set_index("file_name", inplace=True)
    
    valid_samples = [s for s in df_muts.index if s in meta.index]
    
    meta_valid = meta.loc[valid_samples]
    meta_ecuador = meta_valid[meta_valid['country'] == 'Ecuador']
    
    ecuador_sierra = meta_ecuador[meta_ecuador['expected_role'] == 'flu_sierra'].index.tolist()
    ecuador_costa = meta_ecuador[meta_ecuador['expected_role'] == 'flu_costa'].index.tolist()
    ecuador_amazonia = meta_ecuador[meta_ecuador['expected_role'] == 'flu_amazonia'].index.tolist()
    ecuador_other = meta_ecuador[~meta_ecuador['expected_role'].isin(['flu_sierra', 'flu_costa', 'flu_amazonia'])].index.tolist()
    
    ecuador_ids = ecuador_sierra + ecuador_costa + ecuador_amazonia + ecuador_other
    
    avian_ancestors = ["EPI_ISL_20359491", "EPI_ISL_17660072", "EPI_ISL_17777528", "EPI_ISL_18054500", "EPI_ISL_19781426"]
    mammal_ancestor = ["EPI_ISL_18054502", "EPI_ISL_17777531"]
    other_context = ["EPI_ISL_18054509", "EPI_ISL_18265430", "EPI_ISL_18777129", "EPI_ISL_19391462"]
    
    jump_ids = avian_ancestors + mammal_ancestor + other_context
    
    seen = set()
    jump_ids_uniq = [x for x in jump_ids if not (x in seen or seen.add(x))]
    
    all_ids = ecuador_ids + [j for j in jump_ids_uniq if j in meta.index]
    
    if not all_ids:
        return
        
    available_markers = [m for m in CRITICAL_MARKERS if m in df_muts.columns]

    cmap_colors = ['#F8FAFC', '#CBD5E1', '#1D4ED8'] 
    
    state_matrix = pd.DataFrame(0, index=all_ids, columns=available_markers)
    val_matrix = df_muts.loc[all_ids, available_markers].copy()
    
    geo_colors = []
    host_colors = []

    for sample in all_ids:
        geo_colors.append(get_geo_color(meta.loc[sample]))
        host_colors.append(get_host_color(meta.loc[sample]))
        
        for col in available_markers:
            gene, mut = col.split(":")
            ref_aa, alt_aa = mut[0], mut[-1]
            aa = str(val_matrix.loc[sample, col]).upper().strip()
            
            if aa == alt_aa:
                state_matrix.loc[sample, col] = 2
            elif aa in ['?', '-', 'X', 'NAN', '']:
                state_matrix.loc[sample, col] = 1 
                val_matrix.loc[sample, col] = ''
            else:
                state_matrix.loc[sample, col] = 0 

    n_rows, n_cols = state_matrix.shape
    
    fig_w = max(60, 5.0 * n_cols + 20.0) 
    fig_h = max(42, 0.75 * n_rows + 15.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)

    cmap = ListedColormap(cmap_colors)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    
    # Razor thin left margin!
    LEFT_B = -1.8 
    
    ax.imshow(state_matrix.values, aspect='auto', cmap=cmap, norm=norm, extent=[-0.5, n_cols - 0.5, n_rows - 0.5, -0.5])
    ax.set_xlim(LEFT_B, n_cols - 0.5)
    
    # MASSIVELY INCREASED BOTTOM ROW HEIGHT to accommodate fontsize 55 (3 data units tall now!)
    ax.set_ylim(n_rows + 2.5, -2.5)

    for i in range(n_rows):
        for j in range(n_cols):
            state = state_matrix.iat[i, j]
            aa = val_matrix.iat[i, j]
            color = 'white' if state == 2 else '#1F2933'
            ax.text(j, i + 0.05, aa, ha='center', va='center', fontsize=45, color=color, fontweight='bold')

    # Mutation Names centered perfectly in their new 3-unit tall box! (y = n_rows + 1.0)
    for j, col in enumerate(state_matrix.columns):
        label = col.replace(":", "\n")
        ax.text(j, n_rows + 1.0, label, ha='center', va='center', fontsize=55, color='#1F2933', fontweight='bold')

    ax.plot([LEFT_B, LEFT_B], [-0.5, n_rows - 0.5], color='black', linewidth=0.5, clip_on=False, zorder=5) 
    # Extended right wall down to +2.5
    ax.plot([n_cols - 0.5, n_cols - 0.5], [-0.5, n_rows + 2.5], color='black', linewidth=0.5, clip_on=False, zorder=5) 

    ax.plot([LEFT_B, n_cols - 0.5], [-0.65, -0.65], color='black', linewidth=8, clip_on=False, zorder=7)
    ax.plot([LEFT_B, n_cols - 0.5], [-2.5, -2.5], color='black', linewidth=8, clip_on=False, zorder=7)

    for r in range(n_rows + 1):
        ax.plot([LEFT_B, n_cols - 0.5], [r - 0.5, r - 0.5], color='black', linewidth=0.5, clip_on=False, zorder=4)
        
    # Bottom closing line moved to +2.5
    ax.plot([-0.5, n_cols - 0.5], [n_rows + 2.5, n_rows + 2.5], color='black', linewidth=0.5, clip_on=False, zorder=4)
        
    # Vertical grid extending down to +2.5
    for c in range(n_cols):
        ax.plot([c - 0.5, c - 0.5], [-0.5, n_rows + 2.5], color='black', linewidth=0.5, clip_on=False, zorder=4)

    ax.text((LEFT_B - 0.5) / 2.0, -1.4, "Sample", ha='center', va='center', fontsize=70, fontweight='bold', color='#1F2933')
    ax.text((n_cols - 1) / 2.0, -1.4, "Mutation Panel", ha='center', va='center', fontsize=80, fontweight='bold', color='#1F2933')

    idx_split = len(ecuador_ids)
    if idx_split > 0 and idx_split < n_rows:
        ax.plot([LEFT_B, n_cols - 0.5], [idx_split - 0.5, idx_split - 0.5], color='black', linewidth=5, clip_on=False, zorder=6)

    for i, (sample, g_color, h_color) in enumerate(zip(all_ids, geo_colors, host_colors)):
        width = abs(LEFT_B) - 0.5
        rect = mpatches.Rectangle((LEFT_B, i - 0.5), width, 1.0, color=h_color, alpha=0.30, clip_on=False, zorder=0)
        ax.add_patch(rect)
        
        ax.text(-0.7, i + 0.05, sample, ha='right', va='center', fontsize=40, color=g_color, fontweight='bold', clip_on=False, zorder=3)
        
        if sample in avian_ancestors:
            ax.text(-0.4, i + 0.28, "*", ha='center', va='center', fontsize=65, color='red', fontweight='bold', clip_on=False, zorder=6)
        elif sample in mammal_ancestor:
            ax.text(-0.4, i + 0.28, "*", ha='center', va='center', fontsize=65, color='#F59E0B', fontweight='bold', clip_on=False, zorder=6)

    ax.set_xticks([])
    ax.set_yticks([]) 
    ax.axis('off') 

    geo_handles = [
        mpatches.Patch(color='none', label='Ecuador (Andean)'),
        mpatches.Patch(color='none', label='Ecuador (Coastal)'),
        mpatches.Patch(color='none', label='Ecuador (Amazonian)'),
        mpatches.Patch(color='none', label='Peru'),
        mpatches.Patch(color='none', label='Chile')
    ]
    host_handles = [
        mpatches.Patch(color=HOST_COLOR_MAP['domesticated bird'], label='Domesticated Bird'),
        mpatches.Patch(color=HOST_COLOR_MAP['wild bird'], label='Wild Bird'),
        mpatches.Patch(color=HOST_COLOR_MAP['mammal'], label='Mammal (Sea Lion)')
    ]
    ast_handles = [
        mlines.Line2D([], [], color='none', marker='*', markerfacecolor='red', markeredgecolor='none', markersize=60, label='Avian Intro. Ancestor'),
        mlines.Line2D([], [], color='none', marker='*', markerfacecolor='#F59E0B', markeredgecolor='none', markersize=60, label='Close Mammal Relative'),
        mpatches.Patch(color='#1D4ED8', label='Target Mutation'),
        mpatches.Patch(color='#F8FAFC', label='Baseline / WT'),
        mpatches.Patch(color='#CBD5E1', label='Missing Data')
    ]

    leg1 = ax.legend(handles=geo_handles, loc='upper center', bbox_to_anchor=(0.12, -0.01), ncol=2, frameon=False, fontsize=65, title="Location", title_fontproperties={'weight':'bold', 'size':80}, alignment='left', handlelength=0, handletextpad=0.1, columnspacing=0.6, labelspacing=0.05, borderpad=0.0)
    for text, color in zip(leg1.get_texts(), [GEO_COLOR_MAP['flu_sierra'], GEO_COLOR_MAP['flu_costa'], GEO_COLOR_MAP['flu_amazonia'], GEO_COLOR_MAP['Peru'], GEO_COLOR_MAP['Chile']]):
        text.set_color(color)
        text.set_fontweight('bold')
    
    leg2 = ax.legend(handles=host_handles, loc='upper center', bbox_to_anchor=(0.45, -0.01), ncol=1, frameon=False, fontsize=65, title="Host Type", title_fontproperties={'weight':'bold', 'size':80}, alignment='left', handletextpad=0.2, labelspacing=0.05, borderpad=0.0)
    leg3 = ax.legend(handles=ast_handles, loc='upper center', bbox_to_anchor=(0.82, -0.01), ncol=2, frameon=False, fontsize=65, title="Markers", title_fontproperties={'weight':'bold', 'size':80}, alignment='left', handletextpad=0.2, columnspacing=0.6, labelspacing=0.05, borderpad=0.0)
    
    ax.add_artist(leg1)
    ax.add_artist(leg2)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.15, left=0.05, top=0.95, right=0.95) 
    fig.savefig(args.out_png, bbox_inches='tight', bbox_extra_artists=(leg1, leg2, leg3), pad_inches=0.4, facecolor='white')

if __name__ == '__main__':
    main()
