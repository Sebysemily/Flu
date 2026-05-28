# Validated pipeline snapshot (repo)

Living master plan may also exist in Obsidian (`flu_metodos_vivo.md`). This file tracks what the **Git repo** actually runs.

## 1.1 Standardized input ingestion

- **Script:** `code/build_inputs/process_raw_to_segments.py`
- **Ecuador master metadata:** `metadata/flu_filtrado.csv`
- **Panel metadata (derived):** `metadata/H5N1_context.csv` — Ecuador + regional context, one row per `EPI_ISL`
- **Segment FASTAs:** `data/phylogeny/by_segment/H5N1_{segment}.fasta` with `EPI_ISL` headers

## 2.23 Panel pruning

- **Metadata (roles, dates):** `metadata/H5N1_context.csv` only — not `flu_filtrado` at prune time, not `panel_main_taxa.filtered.tsv`
- **Who is in the panel:** tip labels of `data/phylogeny/rtt/panel_main_concat.filtered.nwk` (output of `prune_panel_by_distance.py`)
- **Optional audit:** pass `--panel-main-out` to the prune script if you want a legacy TSV; Snakemake no longer produces it

## 2.5 Segment analysis (ggtree only)

- **Rule:** `segment_analysis_ggtree_fast_ha` → `results/figures/HA_fast_tree.png` (fast ML on all QC-filtered HA samples; no full concat)
- **Rule:** `segment_analysis_ggtree` → `results/figures/{segment}_tree.png` (8 segments + `subset_concat`)
- **Styling:** `code/segment_analysis/tree_aesthetics.R` sourced by `plot_ggtree.R`
- **Trees:** IQ-TREE `rep1` copied to `{segment}_final.treefile` / `subset_concat_final.treefile` (no dedicated rooting step)

### Removed from validated pipeline

- HA/PB2 tanglegram (`draw_tanglegrams.py`)
- iTOL color templates (`generate_itol_colors.py`)
- Verbose-label metadata tools (`build_main_panel_metadata.py`, `relabel_tree_tips.py`)
- Scaled RF distance matrix (`run_rf_matrix.py`)
- American-anchor rooting rule on final trees

### Backlog (manual, not in `rule all`)

- Tanglegrams (`ggtree::cophylo` or similar)
- Heatmap / annotation overlays on trees
- Publication-specific layouts per figure
