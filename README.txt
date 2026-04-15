Pipeline para consolidar filodinámica de influenza A H5N1 en Ecuador con contexto regional y un flujo pre-BEAST basado en el concatenado completo.

Requisitos
- `snakemake`
- `conda` o `mamba`
- conexión a internet si se quiere descargar el contexto desde NCBI

Ejecución
- Preparar inputs de Ecuador desde MIRA:

```bash
snakemake --cores all --use-conda build_gisaid
```

- Correr el pipeline principal completo:

```bash
snakemake --cores all --use-conda
```

- Correr solo la parte pre-BEAST:

```bash
snakemake --cores all --use-conda results/beast_pre/model_test/panel_main_concat.best_model.txt
```

- Correr la validación RF/codón:

```bash
snakemake --cores all --use-conda rf_validation
```

Flujo actual

Stage 1. Preparación de Ecuador desde MIRA
1. `build_ecuador_intermediate_input`
   - Copia y normaliza `amended_consensus.fasta` desde `run*/`
   - Detecta segmentos disponibles por muestra
   - Produce FASTA por muestra y tablas intermedias en `data/assembled/`

2. `build_h5n1_ec_fasta`
   - Toma solo los segmentos ensamblados de Ecuador
   - Usa `flu_filtrado.csv` para metadata
   - Produce:
     - `data/input/H5N1_EC.fasta`
     - `data/input/H5N1_EC_summary.csv`

Stage 2. Contexto, alineación y árboles ML
3. `build_h5n1_final_fasta`
   - Descarga secuencias contextuales por accession
   - Renombra headers al mismo estilo usado en Ecuador
   - Une Ecuador + contexto en:
     - `data/input/H5N1_context.fasta`
     - `data/input/H5N1_context_summary.csv`
     - `data/final/H5N1_final.fasta`

4. `split_h5n1_final_by_segment`
   - Divide el FASTA final en 8 segmentos

5. `mafft_align_per_segment`
   - Alinea cada segmento con MAFFT

6. `build_segment_codon_partitions`
   - Genera particiones codón para PB2, PB1, PA, HA, NP y NA

7. `raxml_ng_tree_per_segment_codon`
   - Corre ML + bootstrap + TBE para PB2, PB1, PA, HA, NP y NA

8. `raxml_ng_tree_per_segment_simple`
   - Corre ML + bootstrap + TBE para NS y MP

9. `concat_aligned_segments_with_partitions`
   - Concatena los 8 segmentos y conserva las particiones

10. `raxml_ng_tree_full_concat`
   - Corre el árbol ML concatenado completo con soporte TBE

Stage 3. Pre-BEAST sobre el concatenado subseteado
11. `build_beast_panels`
   - Usa el árbol concatenado completo para elegir el panel principal
   - El dataset por defecto es:
     - `ecuador_core`
     - `regional_context`
     - `OQ968009`
     - `1` `american_anchor` extra
     - sin `usa_distal` por defecto
   - Produce:
     - `data/beast_pre/panels/panel_main_taxa.tsv`
     - `data/beast_pre/panels/panel_selection_audit.tsv`
     - `data/beast_pre/panels/panel_country_month_coverage.tsv`

12. `subset_panel_concat_alignment_and_prune_tree`
   - Extrae del concatenado completo solo los taxa del panel
   - Produce:
     - `data/beast_pre/panels/panel_main_concat.subset.fasta`
     - `data/beast_pre/panels/panel_main_concat.subset.nwk`

13. `build_treetime_dates`
   - Construye la tabla de fechas desde los headers

14. `run_root_to_tip`
   - Corre TreeTime root-to-tip sobre el subset final

15. `ensure_beast_model_selection_package`
   - Instala el paquete `MODEL_SELECTION` en un directorio local de BEAST

16. `run_bets_scenario` y `summarize_bets`
   - Ejecutan BETS a partir de XMLs externos si está habilitado
   - Resumen en:
     - `results/beast_pre/bets/bets_runs.tsv`
     - `results/beast_pre/bets/bets_bayes_factors.tsv`

17. `build_concat_iqtree_partition_nexus`
   - Convierte las particiones del concatenado a NEXUS para IQ-TREE

18. `iqtree_model_test_concat_subset`
   - Corre ModelFinder de IQ-TREE sobre el subset final
   - Produce:
     - `results/beast_pre/model_test/panel_main_concat.iqtree`
     - `results/beast_pre/model_test/panel_main_concat.best_model.txt`

Validación RF opcional
- `rf_validation` genera comparaciones RF para:
  - 6 segmentos con particiones codón
  - NS y MP con un esquema simple

Cómo funciona el modo “relajado” del subset pre-BEAST
- Primero corre una selección estricta de `regional_context`:
  - exige que el taxón comparta con `ecuador_core` un MRCA no raíz
  - ese MRCA debe cumplir `min_mrca_support`
  - además reparte la muestra por país y por mes usando `n_per_country`, `n_total` y `max_per_country_month`

- Después, solo si `relaxed_fill > 0`, corre una segunda pasada:
  - vuelve a mirar solo `regional_context`
  - usa un umbral menor: `relaxed_min_mrca_support`
  - añade hasta `relaxed_fill` taxones extra
  - sigue respetando el tope por país `n_per_country`
  - no reintroduce `usa_distal` ni anchors dentro de esa pasada

Interpretación práctica
- Si subes `relaxed_fill`, metes más contexto regional
- Si bajas `relaxed_min_mrca_support`, permites clados algo menos robustos
- Si dejas `usa_distal_quota: 0`, ese contexto extra sigue siendo regional y no distal USA

Configuración

El archivo central es `config/config.yml`.

Variables globales
- `mira_base_dir`: directorio base donde viven `run*/`
- `flu_filtrado`: metadata/filtro de Ecuador
- `context_metadata_tsv`: tabla curada del contexto regional
- `ecuador_date_source`: fecha que se usa para Ecuador en headers y resúmenes
- `random_seed`: semilla base reproducible
- `max_threads`: techo general de hilos
- `mafft_threads`: hilos para MAFFT
- `raxml_segment_threads`: hilos para RAxML-NG por segmento
- `raxml_full_concat_threads`: hilos para RAxML-NG concatenado
- `raxml_rf_threads`: hilos para validación RF
- `iqtree_modeltest_threads`: hilos para ModelFinder en pre-BEAST

Variables de `beast_pre`
- `ecuador_core_ids`: lista de muestras Ecuador que define el clado focal
- `regional_blacklist_tokens`: substrings que excluyen taxa del pool `regional_context`
- `max_cluster_dist`: parámetro heredado; hoy no controla la selección
- `n_per_country`: máximo de `regional_context` por país en la pasada estricta
- `n_total`: máximo total de `regional_context` en la pasada estricta
- `min_mrca_support`: soporte TBE mínimo para aceptar candidatos regionales en la pasada estricta
- `max_per_country_month`: máximo por combinación país-mes en la pasada estricta
- `usa_distal_quota`: cuántos `usa_distal` agregar por distancia mínima a Ecuador
- `forced_american_anchor_accession`: accession que siempre entra como `american_anchor` si está presente
- `additional_american_anchor_quota`: anchors adicionales aparte del accession forzado
- `relaxed_min_mrca_support`: soporte TBE mínimo de la pasada relajada
- `relaxed_fill`: número de taxa extra que añade la pasada relajada

Variables de `beast_pre.bets`
- `enabled`: fuerza activar o desactivar BETS
- `threads`: hilos usados por BEAST/BETS
- `package_dir`: directorio local para instalar `MODEL_SELECTION`
- `*_xml`: rutas a los XMLs ya preparados para cada escenario BETS

Perfiles útiles del panel pre-BEAST

Perfil actual por defecto
- `usa_distal_quota: 0`
- `forced_american_anchor_accession: "OQ968009"`
- `additional_american_anchor_quota: 1`
- `relaxed_fill: 0`

Más contexto regional sin abrir mucho el panel
- `relaxed_fill: 5`
- `relaxed_min_mrca_support: 50.0`

Más permisivo todavía
- `relaxed_fill: 10`
- `relaxed_min_mrca_support: 40.0`

Salidas principales
- `data/input/H5N1_EC.fasta`
- `data/final/H5N1_final.fasta`
- `data/phylogeny/aligned/H5N1_full_concat_beast.mafft`
- `results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE`
- `data/beast_pre/panels/panel_main_taxa.tsv`
- `data/beast_pre/panels/panel_main_concat.subset.fasta`
- `results/beast_pre/root_to_tip/treetime_clock.log`
- `results/beast_pre/bets/bets_runs.tsv`
- `results/beast_pre/model_test/panel_main_concat.best_model.txt`

Notas
- El FASTA final contiene primero Ecuador y luego el contexto.
- El subset pre-BEAST se construye sobre el concatenado ya alineado; no se realinea.
- `panel_selection_audit.tsv` resume cuántos taxa entraron por cada mecanismo.
- `panel_country_month_coverage.tsv` ayuda a revisar si el contexto quedó balanceado en tiempo y país.
