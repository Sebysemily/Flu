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
# Construir/normalizar inputs de Ecuador y luego headers con fechas.
snakemake --cores all --use-conda build_gisaid

# Descargar/mezclar contexto, alinear, arboles ML por segmento y concatenado.
snakemake --cores all --use-conda main_phylogeny

# Seleccion/QC del panel, root-to-tip, FASTA finales por segmento y XML de BEAST.
snakemake --cores all --use-conda pre_beast_outputs

# Flujo principal completo hasta corridas BEAST y arboles finales combinados.
snakemake --cores all --use-conda beast_runs
```

Nota: `snakemake --cores all --use-conda` puede incluir targets historicos de sensibilidad NP/MP/NS definidos en `rule all`. Para reproducir el flujo principal actual sin ese subset, usar `beast_runs`.

Flujo actual
------------
1. Entrada Ecuador: `data/input/` queda reservado para FASTA descargados de GISAID. Puede contener uno o varios archivos con nombres libres (`*.fasta`, `*.fa`, `*.fna`, `*.fas`). Si no hay FASTA en esa carpeta, el pipeline puede construir uno equivalente desde `run*/amended_consensus.fasta` bajo `mira_base_dir` y lo escribe como `data/assembled/H5N1_EC_gisaid_from_mira.fasta`.
2. Headers de trabajo: `build_standard_date_headers.py` toma los FASTA tipo GISAID y genera `data/assembled/H5N1_EC.fasta`/`data/assembled/H5N1_EC_summary.csv` con headers `Flu-XXXX/segment/province/date`.
3. Contexto regional: usa `config/final_metadata_50_per_country_isolates.tsv`, descarga secuencias de contexto (sin filtro de redundancia por país-mes `max_per_country_month` para preservar la diversidad nativa) y produce `data/assembled/H5N1_context.fasta`, `data/assembled/H5N1_context_summary.csv` y `data/final/H5N1_final.fasta`.
4. Filogenia ML: separa los 8 segmentos (`PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `MP`, `NS`), alinea con Nextclade, construye particiones por codon para `PB2/PB1/PA/HA/NP/NA`, usa `GTR+G` simple para `MP/NS`, corre IQ-TREE por segmento y replica, y arma el concatenado completo `data/phylogeny/aligned/H5N1_full_concat.mafft` con `results/phylogeny/iq-tree/full_concat/H5N1_full_concat.treefile`.
5. Pre-BEAST: selecciona el panel mediante un filtrado por proximidad patristica (Dijkstra) en el árbol concatenado. Protege por mes los `protect_anchors_per_month` (default: 3) anchors de Norteamérica más cercanos al core de Ecuador hasta mayo 2022, y los `protect_regional_per_month` (default: 3) contextos regionales más cercanos al core hasta la fecha más reciente. Luego calcula el clock de TreeTime (filtrando outliers según `clock_filter` en `config.yml`, default: 3.0).
6. Entregables pre-BEAST: genera la filogenia y distancias del subset concatenado (`results/phylogeny/iq-tree/subset_concat/`) con soporte de bootstrap (`iqtree_bootstrap` en `config.yml`) y distancias Robinson-Foulds intra-replicate. Produce `data/phylogeny/rtt/panel_main_concat.filtered.fasta`, `data/phylogeny/rtt/panel_main_concat.filtered.nwk` y los FASTA finales por segmento para BEAST.
7. XML BEAST: parametriza templates versionados en `template_beast/` y escribe XML en `results/beast/xml/` para `strict_constant`, `strict_constant_lugar`, `ucln_constant`, `strict_exp` y `ucln_exp`.
8. Corridas BEAST: valida cada XML, corre dos replicas (`r1`, `r2`) por escenario con seeds de `config.beast.seeds` y `seed + 100000` para `r2`, guarda logs/trees/chkpt/ops en `results/beast/runs/<scenario>/<replicate>/` y resume cada escenario con `results/beast/runs/<scenario>/run.done`.
9. Finales combinados: combina las dos replicas de `strict_constant` y `strict_constant_lugar` con 10% de burn-in y anota arbol MCC con alturas medias. Salidas: `results/beast/final/time/strict_constant.*` y `results/beast/final/geography/strict_constant_lugar.*`.

Salidas principales
-------------------
- Entrada GISAID fuente: FASTA descargados por el usuario en `data/input/`.
- FASTA derivados: `data/assembled/H5N1_EC_gisaid_from_mira.fasta`, `data/assembled/H5N1_EC.fasta`, `data/assembled/H5N1_context.fasta`, `data/final/H5N1_final.fasta`.
- Arboles ML: `results/phylogeny/iq-tree/<segment>/rep<rep>.treefile`, `results/phylogeny/iq-tree/full_concat/H5N1_full_concat.treefile` y `results/phylogeny/iq-tree/subset_concat/subset_concat_final.treefile`.
- QC y RTT: `results/phylogeny/rtt/treetime_clock.done`, `results/qc_metrics/alignments/pruning_discarded_sequences.csv`.
- Panel final: `data/phylogeny/rtt/panel_main_taxa.filtered.tsv`, `data/phylogeny/rtt/panel_main_concat.filtered.fasta`, `data/phylogeny/main_panel/H5N1_{PB2,PB1,PA,HA,NP,NA,MP,NS}.fasta`.
- XML: `results/beast/xml/{strict_constant,strict_constant_lugar,ucln_constant,strict_exp,ucln_exp}.xml`.
- BEAST: `results/beast/runs/<scenario>/r{1,2}/` y `results/beast/runs/<scenario>/run.done`.
- Finales BEAST: `results/beast/final/time/strict_constant.combined.{log,trees}`, `results/beast/final/time/strict_constant.mcc.mean.tree`, `results/beast/final/geography/strict_constant_lugar.combined.{log,trees}`, `results/beast/final/geography/strict_constant_lugar.mcc.mean.tree`.

Configuracion
-------------
El archivo central es `config/config.yml`.

- Entradas: `mira_base_dir`, `flu_filtrado`, `context_metadata_tsv`, `ecuador_date_source`.
- `flu_filtrado` debe incluir `EPI_ISL` para cada `Código USFQ`; las columnas `PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `MP`, `NS` son la fuente simple de segmentos esperados.
- Recursos/reproducibilidad: `random_seed`, `max_threads`, `iqtree_replicates`, `iqtree_bootstrap`.
- Prunning y RTT: `prunning_parameters.n_closest`, `prunning_parameters.max_distance`, `prunning_parameters.protect_anchors_per_month`, `prunning_parameters.protect_regional_per_month`, `treetime_parameters.clock_filter`.
- BEAST: `beast.enabled`, `binary`, `max_hours`, `threads`, `chain_length`, `log_every`, `tree_every`, `echo_every`, `seeds` y bloque `beagle`.

Notas operativas
----------------
- Para usar descargas directas de GISAID, poner uno o varios FASTA en `data/input/`. Esa carpeta no debe usarse para archivos generados por el pipeline.
- El FASTA GISAID-standard de Ecuador usa headers `virus_name|segment|EPI_ISL`, por ejemplo `A/chicken/Ecuador/Flu-0610/2023|NP|EPI_ISL_20450104`.
- `build_standard_date_headers.py` acepta FASTA descargados de GISAID si los headers vienen como `virus_name|segment|EPI_ISL` o `virus_name|EPI_ISL|segment`. No acepta el formato viejo `Flu-XXXX|A_SEG`.
- Los templates de BEAST viven en `template_beast/`; `prepare_beast_run_xml.py` solo actualiza parametros de corrida y nombres de salida.
- `beast.max_hours` se aplica como recurso de Snakemake; no corta el XML por tiempo.
- `beast.beagle.mode` acepta `off`, `auto` o `force`. En `auto`, el runner consulta `beast -beagle_info` y usa BEAGLE solo si detecta recurso compatible.
- Para GPU, `beast.beagle.vendor: amd` favorece OpenCL con `platform: auto`; `nvidia` favorece CUDA. Si no hay recurso compatible y `fallback_to_cpu: false`, la corrida falla en vez de cambiar silenciosamente a CPU.
