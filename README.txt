Pipeline para consolidar filodinámica de influenza A H5N1 en Ecuador con contexto regional.

El flujo usa el concatenado solo como artefacto interno para selección de panel, QC y root-to-tip. Los alineamientos finales para revisión manual siguen siendo los FASTA por segmento en `data/beast/final_panel_segment/`, y ahora el workflow además prepara y corre cuatro escenarios exploratorios de BEAST a partir de templates XML versionados.

Requisitos
- `snakemake`
- `conda` o `mamba`
- conexión a internet si se quiere descargar contexto desde NCBI
- un binario compatible con BEAST 1 disponible en `PATH` o configurable con `config.beast.binary`

Ejecución
- Preparar inputs de Ecuador desde MIRA:

```bash
snakemake --cores all --use-conda build_gisaid
```

- Correr el workflow principal hasta las corridas exploratorias de BEAST:

```bash
snakemake --cores all --use-conda
```

- Correr solo el flujo pre-BEAST con XML preparados:

```bash
snakemake --cores all --use-conda pre_beast_outputs
```

- Correr solo las corridas exploratorias de BEAST:

```bash
snakemake --cores all --use-conda beast_runs
```

Resumen del flujo

Stage 1. Preparación de Ecuador desde MIRA
1. `build_ecuador_intermediate_input`
2. `build_h5n1_ec_fasta`

Stage 2. Contexto, alineación y árboles ML
3. `build_h5n1_final_fasta`
4. `split_h5n1_final_by_segment`
5. `mafft_align_per_segment`
6. `build_segment_codon_partitions`
7. `raxml_ng_tree_per_segment_codon`
8. `raxml_ng_tree_per_segment_simple`
9. `concat_aligned_segments_with_partitions`
10. `raxml_ng_tree_full_concat`

Stage 3. Pre-BEAST con panel, QC y root-to-tip
11. `build_beast_panels`
12. `observe_beast_subset_source_qc`
13. `filter_beast_panel_by_qc`
14. `subset_filtered_panel_concat_alignment_and_prune_tree`
15. `build_treetime_dates` y `run_root_to_tip`
   - Generan un unico RTT visible en `results/beast_pre/rtt/`
16. `filter_beast_panel_by_rtt_outliers`
   - Usa el RTT para excluir outliers y escribir el TSV final de fechas
17. `publish_final_panel_concat_alignment`
18. `subset_filtered_raw_segment_fasta`
19. `subset_filtered_segment_alignment`
20. `observe_segment_alignment_qc`
21. `summarize_final_segment_qc`
22. `prepare_beast_run_xml`
   - Genera:
     - `results/beast/xml/strict_constant.xml`
     - `results/beast/xml/ucln_constant.xml`
     - `results/beast/xml/strict_exp.xml`
     - `results/beast/xml/ucln_exp.xml`

Stage 4. Corridas exploratorias de BEAST
23. `validate_beast_xml`
24. `run_beast_replicate`
   - Corre cada escenario por duplicado (`r1`, `r2`)
   - Usa seed base y seed base `+100000`
25. `summarize_beast_run`
   - Genera:
     - `results/beast/runs/strict_constant/run.done`
     - `results/beast/runs/ucln_constant/run.done`
     - `results/beast/runs/strict_exp/run.done`
     - `results/beast/runs/ucln_exp/run.done`

Outputs principales

QC y root-to-tip
- `results/beast_pre/qc_validation/source_panel_qc.metrics.tsv`
- `results/beast_pre/qc_validation/source_panel_qc.summary.tsv`
- `results/beast_pre/qc_validation/final_segment_panel_qc.summary.tsv`
- `results/beast_pre/rtt/treetime_clock.done`
- `data/beast/panel_main_dates.final.tsv`

Entregables finales por segmento
- `data/beast/final_panel_segment/H5N1_PB2.fasta`
- `data/beast/final_panel_segment/H5N1_PB1.fasta`
- `data/beast/final_panel_segment/H5N1_PA.fasta`
- `data/beast/final_panel_segment/H5N1_HA.fasta`
- `data/beast/final_panel_segment/H5N1_NP.fasta`
- `data/beast/final_panel_segment/H5N1_NA.fasta`
- `data/beast/final_panel_segment/H5N1_MP.fasta`
- `data/beast/final_panel_segment/H5N1_NS.fasta`

XML y corridas exploratorias de BEAST
- `results/beast/xml/strict_constant.xml`
- `results/beast/xml/ucln_constant.xml`
- `results/beast/xml/strict_exp.xml`
- `results/beast/xml/ucln_exp.xml`
- `results/beast/runs/strict_constant/r1/run.done`
- `results/beast/runs/strict_constant/r2/run.done`
- `results/beast/runs/strict_constant/run.done`
- `results/beast/runs/ucln_constant/r1/run.done`
- `results/beast/runs/ucln_constant/r2/run.done`
- `results/beast/runs/ucln_constant/run.done`
- `results/beast/runs/strict_exp/r1/run.done`
- `results/beast/runs/strict_exp/r2/run.done`
- `results/beast/runs/strict_exp/run.done`
- `results/beast/runs/ucln_exp/r1/run.done`
- `results/beast/runs/ucln_exp/r2/run.done`
- `results/beast/runs/ucln_exp/run.done`

Configuración

El archivo central es `config/config.yml`.

Variables globales
- `mira_base_dir`
- `flu_filtrado`
- `context_metadata_tsv`
- `ecuador_date_source`
- `random_seed`
- `max_threads`
- `mafft_threads`
- `raxml_segment_threads`
- `raxml_full_concat_threads`

Variables de `beast_pre`
- `relaxed_panel_mode`

Variables de `beast`
- `enabled`
- `binary`
- `max_hours`
- `threads`
- `beagle.mode`
- `beagle.vendor`
- `beagle.resource`
- `beagle.platform`
- `beagle.precision`
- `beagle.scaling`
- `beagle.fallback_to_cpu`
- `beagle.info`
- `chain_length`
- `log_every`
- `tree_every`
- `echo_every`
- `seeds.strict_constant`
- `seeds.ucln_constant`
- `seeds.strict_exp`
- `seeds.ucln_exp`

Notas
- Los templates XML base viven en `template_beast/` y son la fuente versionada para las corridas exploratorias.
- El script `prepare_beast_run_xml.py` solo parametriza la corrida; no reescribe el contenido científico del template.
- El scheduler del workflow aplica el tope operativo de tiempo por corrida con `resources.runtime`; el XML no se corta por tiempo.
- `config.beast.beagle.mode` acepta `off`, `auto` o `force`.
- Con `mode: auto`, el runner consulta `beast -beagle_info` y solo agrega BEAGLE si detecta un recurso compatible; si no, corre sin BEAGLE.
- `config.beast.beagle.vendor` permite preferir `amd`, `nvidia` o `any`; con `platform: auto`, `amd` favorece `opencl` y `nvidia` favorece `cuda`.
- `config.beast.seeds` define la seed base por escenario; la repeticion `r2` se deriva automaticamente como `seed_base + 100000`.
- Si el binario de BEAST 1 no está en `PATH`, define `config.beast.binary`.
