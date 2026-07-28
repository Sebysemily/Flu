H5N1 Ecuador Pipeline

Workflow for the phylogenetic analysis of H5N1 in Ecuador, associated with the analysis of doi..

Requirements
----------

- `snakemake` with `conda` or `mamba`.
- Internet connection to install the environments.

Useful Commands
---------------

```bash

# Complete main flow (rule all runs the entire analysis by default).
snakemake --cores all --use-conda

# To re-run model selection (GSS), use the GSS target...
snakemake --cores all --use-conda gss_model_selection
```

Current Flow
------------

1. **Inputs:** An `epi_set` file in `data/context` (which must have at least the columns `Virus name`, `Isolate ID`, `Collection date`, and `Segment`). For the Ecuador input: FASTA in `data/input/` (GISAID) or constructed from MIRA (`mira_base_dir`) → `data/input/H5N1_EC_gisaid_from_mira.fasta`. Master metadata in `metadata/flu_filtrado.csv`.
2. **Unified Ingestion:** `process_raw_to_segments` writes `data/phylogeny/by_segment/H5N1_{segment}.fasta` (headers `EPI_ISL`) and `metadata/H5N1_context.csv` (Ecuador + regional context, deduplicated by taxon).
3. **ML Phylogeny:** Nextclade by segment and QC across the 8 segments; fast HA tree (`iqtree_fast_ha` on `H5N1_HA.mafft` with all samples) to iterate `plot_ggtree.R`. The panel/subset still uses internal concat in `apply_panel_to_segment_alignments`.
4. **Panel:** `prune_panel_by_distance.py` uses roles from `H5N1_context.csv` (Ecuador core, American anchors, regional context) → `data/phylogeny/main_panel/`.
5. **Subset + RTT:** IQ-TREE on the concatenated panel; TreeTime with dates from `H5N1_context.csv`.
6. **Figures:** `plot_ggtree.R` + `tree_aesthetics.R` produce `results/figures/{segment}_tree.png` (8 segments + `subset_concat`) with a unified style.
7. **BEAST:** templates in `template_beast/`, runs in `results/beast/`.

Main Outputs
-------------------

- Processed segments: `data/phylogeny/by_segment/H5N1_{segment}.fasta`
- Final metadata: `metadata/H5N1_context.csv`, Ecuador master: `metadata/flu_filtrado.csv`
- ML Trees (IQ-TREE): `results/phylogeny/iq-tree/`
- ML tree figures per segment and composite: `figures/HA_PB2_lineage/` and `figures/`
- Model Selection (GSS): `results/beast/GSS/`
- BEAST MCC tree and final results: `results/beast/final_run/`
- Main BEAST figure: `figures/main_panel_HA_beast_mcc.png`
- QC: `results/qc_metrics/alignments/pruning_discarded_sequences.csv`
- RTT Panel: `data/phylogeny/rtt/panel_main_concat.filtered.nwk` (tip list = pruned panel)
- Ecological Groups and Reassortments: Classification files in `results/clasifications_ecological_groups.csv` and outlier details in `results/possible_reassortants_ignored.csv`.

Configuration
-------------

Files: `config/config.yml`, `config/paths.yml`.

- Inputs: `mira_base_dir`, `flu_filtrado` (default `metadata/flu_filtrado.csv`), `context_fasta_raw`, `ecuador_date_source`.
- `flu_filtrado` must include `EPI_ISL` per `Código USFQ`; segment columns indicate expected coverage.
- Reproducibility: `random_seed`, `max_threads`, `iqtree_replicates`, `iqtree_bootstrap`.
- Pruning: `prunning_parameters.n_closest`, `max_distance`, `protect_anchors_per_month`, `protect_regional_per_month`.
- TreeTime: `treetime_parameters.clock_filter`.

Visualization (manual backlog)
------------------------------

Tanglegrams, heatmaps, and publication figures will be developed separately in scripts that reuse `code/segment_analysis/tree_aesthetics.R`. They are not part of the current Snakemake targets.

Operational Notes
----------------

- Ecuador headers in GISAID FASTA: `virus_name|segment|EPI_ISL` or `virus_name|EPI_ISL|segment`.
- The pipeline does not generate iTOL colors or RF matrices; a single color scheme in R is used.
- `data/phylogeny/main_panel/` contains alignments of the final panel (do not confuse with removed metadata `metadata/main_panel.csv`).

Reassortments and Outliers (Genotype Analysis)
------------------------------------------------

During the process, samples acting as recombinants or endemic strains are identified and excluded from the main H5N1 visualization due to massive divergence in certain segments (e.g., NS or MP).

The metadata extracted from GenoFLU for these key samples (condors, boobies, otters, etc.), where the discordance of their genotypes against typical HPAI lineages is observed, is consolidated in the `results/ignored_in_per_segment.csv` file.

On the other hand, the classification of these lineages according to their ecological groups (such as Andes, Coast, Amazon) and host types (wild birds, mammals, etc.) for BEAST analyses is detailed in the `results/clasifications_ecological_groups.csv` file.
