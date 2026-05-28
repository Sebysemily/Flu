Pipeline H5N1 Ecuador
=====================

Workflow Snakemake para consolidar genomas de influenza A H5N1 de Ecuador con contexto regional, construir arboles ML por segmento y concatenado, seleccionar un panel para BEAST, hacer QC/root-to-tip, preparar XML de BEAST 1 y correr/combinar escenarios exploratorios.

Requisitos
----------
- `snakemake` con `conda` o `mamba`.
- Internet si se reconstruye el contexto desde NCBI.
- BEAST 1, `logcombiner` y `treeannotator` en `PATH`, o `config.beast.binary` definido.

Comandos utiles
---------------
```bash
# Ingesta unificada: segmentos por EPI_ISL + metadata de panel.
snakemake --cores all --use-conda process_raw_to_segments

# Ingesta, alineacion por segmento, arbol HA rapido y figura ggtree.
snakemake --cores all --use-conda results/figures/HA_fast_tree.png

# Panel por distancia, subset, RTT, figuras ggtree y pre-BEAST.
snakemake --cores all --use-conda -R rtt_and_ml_subset

# Flujo principal completo (targets en snakefile rule all).
snakemake --cores all --use-conda
```

Flujo actual
------------
1. **Entrada Ecuador:** FASTA en `data/input/` (GISAID) o construccion desde MIRA (`mira_base_dir`) → `data/input/H5N1_EC_gisaid_from_mira.fasta`. Metadata maestra en `metadata/flu_filtrado.csv`.
2. **Ingesta unificada:** `process_raw_to_segments` escribe `data/phylogeny/by_segment/H5N1_{segment}.fasta` (headers `EPI_ISL`) y `metadata/H5N1_context.csv` (Ecuador + contexto regional, deduplicado por taxon).
3. **Filogenia ML:** Nextclade por segmento y QC en los 8 segmentos; arbol rapido HA (`iqtree_fast_ha` sobre `H5N1_HA.mafft` con todas las muestras) para iterar `plot_ggtree.R`. El panel/subset sigue usando concat interno en `apply_panel_to_segment_alignments`.
4. **Panel:** `prune_panel_by_distance.py` usa roles de `H5N1_context.csv` (core Ecuador, anchors americanos, contexto regional) → `data/phylogeny/main_panel/`.
5. **Subset + RTT:** IQ-TREE en concatenado del panel; TreeTime con fechas de `H5N1_context.csv`.
6. **Figuras:** `plot_ggtree.R` + `tree_aesthetics.R` producen `results/figures/{segment}_tree.png` (8 segmentos + `subset_concat`) con estilo unificado.
7. **BEAST:** templates en `template_beast/`, corridas en `results/beast/`.

Salidas principales
-------------------
- Segmentos: `data/phylogeny/by_segment/H5N1_{PB2,PB1,PA,HA,NP,NA,MP,NS}.fasta`
- Metadata panel: `metadata/H5N1_context.csv`, maestro Ecuador: `metadata/flu_filtrado.csv`
- Arboles ML: `results/phylogeny/iq-tree/<segment>/rep<rep>.treefile`, `{segment}_final.treefile`, `subset_concat/subset_concat_final.treefile`
- Figura rapida (iteracion ggtree): `results/figures/HA_fast_tree.png` desde alineacion HA QC + `iqtree_fast_ha`
- Figuras panel: `results/figures/{segment}_tree.png`, `results/figures/subset_concat_tree.png`
- QC: `results/qc_metrics/alignments/pruning_discarded_sequences.csv`
- Panel RTT: `data/phylogeny/rtt/panel_main_concat.filtered.nwk` (tip list = pruned panel), metadata en `metadata/H5N1_context.csv`, alineaciones en `data/phylogeny/main_panel/`

Configuracion
-------------
Archivos: `config/config.yml`, `config/paths.yml`.

- Entradas: `mira_base_dir`, `flu_filtrado` (default `metadata/flu_filtrado.csv`), `context_fasta_raw`, `ecuador_date_source`.
- `flu_filtrado` debe incluir `EPI_ISL` por `Código USFQ`; columnas de segmentos indican cobertura esperada.
- Reproducibilidad: `random_seed`, `max_threads`, `iqtree_replicates`, `iqtree_bootstrap`.
- Pruning: `prunning_parameters.n_closest`, `max_distance`, `protect_anchors_per_month`, `protect_regional_per_month`.
- TreeTime: `treetime_parameters.clock_filter`.

Visualizacion (backlog manual)
------------------------------
Tanglegrams, heatmaps y figuras de publicacion se desarrollaran aparte en scripts que reutilicen `code/segment_analysis/tree_aesthetics.R`. No forman parte de los targets Snakemake actuales.

Notas operativas
----------------
- Headers Ecuador en FASTA GISAID: `virus_name|segment|EPI_ISL` o `virus_name|EPI_ISL|segment`.
- El pipeline no genera colores iTOL ni matrices RF; un solo esquema de color en R.
- `data/phylogeny/main_panel/` son alineaciones del panel final (no confundir con metadata eliminada `metadata/main_panel.csv`).
