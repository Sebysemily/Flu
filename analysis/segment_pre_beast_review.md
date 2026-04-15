# Review: Pre-BEAST con root-to-tip, BETS e IQ-TREE sobre el concatenado subseteado

## Objetivo del cambio

Dejar `02_pre_beast.smk` en un flujo simple y congruente con lo conversado:

1. usar el árbol concatenado completo para selección de taxa,
2. construir el dataset final subseteado del concatenado sin realinear,
3. correr root-to-tip con TreeTime sobre ese mismo dataset final,
4. correr BETS desde XMLs externos ya preparados para ese mismo dataset,
5. convertir las particiones ya generadas en `01_ml_trees` a NEXUS para IQ-TREE,
6. correr un solo IQ-TREE particionado sobre ese concatenado subseteado final.

En esta versión no se generan plantillas XML para BETS: `02_pre_beast.smk` solo ejecuta los XMLs que se indiquen en `config/config.yml`.

## Qué se mantiene

- `01_ml_trees.smk` sigue haciendo el ML global y construyendo:
  - `data/phylogeny/aligned/H5N1_full_concat_beast.mafft`
  - `data/phylogeny/H5N1_full_concat_beast.partitions`
  - `results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE`
- `build_beast_panels.py` sigue usando el árbol concatenado completo.
- `subset_alignment_and_prune_tree.py` sigue preservando coordenadas del alignment y ya no realinea con MAFFT.

## Cambios implementados

### 1. `02_pre_beast.smk` reescrito al flujo concat-subset

Archivo: `rules/02_pre_beast.smk`

Ahora el flujo principal es:

- `build_beast_panels`
- `subset_panel_concat_alignment_and_prune_tree`
- `build_treetime_dates`
- `run_root_to_tip`
- `ensure_beast_model_selection_package`
- `run_bets_scenario`
- `summarize_bets`
- `build_concat_iqtree_partition_nexus`
- `iqtree_model_test_concat_subset`

### 2. Dataset final subseteado del concatenado

La nueva regla usa:

- alignment:
  - `data/phylogeny/aligned/H5N1_full_concat_beast.mafft`
- tree:
  - `results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE`
- taxa:
  - `data/beast_pre/panels/panel_main_taxa.tsv`

Y genera:

- `data/beast_pre/panels/panel_main_concat.subset.fasta`
- `data/beast_pre/panels/panel_main_concat.subset.nwk`
- `data/beast_pre/panels/panel_main_concat.subset.audit.tsv`

### 3. Particiones IQ-TREE a partir de las particiones ya hechas en `01`

No se tocaron los scripts de `01_ml_trees`.

En lugar de eso, `02_pre_beast.smk` toma:

- `data/phylogeny/H5N1_full_concat_beast.partitions`

y lo convierte a:

- `data/beast_pre/panels/panel_main_concat.subset.iqtree.nex`

La conversión:

- elimina el prefijo de modelo RAxML,
- conserva exactamente los nombres y rangos de las 14 particiones,
- escribe un archivo NEXUS con `charset ... = DNA, ...;`

Como el subset no realinea ni cambia columnas, las coordenadas siguen siendo válidas.

### 4. Root-to-tip reintroducido dentro de pre-BEAST

Se volvió a añadir:

- `build_treetime_dates`
- `run_root_to_tip`

Usan:

- `data/beast_pre/panels/panel_main_concat.subset.fasta`
- `data/beast_pre/panels/panel_main_concat.subset.nwk`

Y escriben en:

- `results/beast_pre/root_to_tip/dates_from_headers.tsv`
- `results/beast_pre/root_to_tip/treetime_clock.log`
- `results/beast_pre/root_to_tip/treetime_clock.done`

Además, la regla de IQ-TREE ahora depende explícitamente de `treetime_clock.done`, para que root-to-tip corra antes.

### 5. BETS añadido como paso ejecutable, sin generar XMLs

Se añadió una rama de BETS dentro de `02_pre_beast.smk` que:

- instala el paquete `MODEL_SELECTION` de BEAST en:
  - `results/beast_pre/bets/packages`
- ejecuta cuatro XMLs externos:
  - `strict + constant + heterochronous`
  - `strict + constant + isochronous`
  - `UCLN + constant + heterochronous`
  - `UCLN + constant + isochronous`
- resume los marginal likelihoods y los log Bayes factors en:
  - `results/beast_pre/bets/bets_runs.tsv`
  - `results/beast_pre/bets/bets_bayes_factors.tsv`

Los XMLs no se generan en el pipeline. Se esperan en estas rutas configurables:

- `config/beast_pre/bets/strict_constant_heterochronous.xml`
- `config/beast_pre/bets/strict_constant_isochronous.xml`
- `config/beast_pre/bets/ucln_constant_heterochronous.xml`
- `config/beast_pre/bets/ucln_constant_isochronous.xml`

### 6. Un solo IQ-TREE particionado sobre el dataset final

La regla `iqtree_model_test_concat_subset` corre:

```bash
iqtree2 \
  -s data/beast_pre/panels/panel_main_concat.subset.fasta \
  -p data/beast_pre/panels/panel_main_concat.subset.iqtree.nex \
  -m MFP \
  -mset HKY,GTR \
  -pre results/beast_pre/model_test/panel_main_concat
```

Salidas:

- `results/beast_pre/model_test/panel_main_concat.iqtree`
- `results/beast_pre/model_test/panel_main_concat.best_model.txt`

Además, ahora IQ-TREE depende del resumen de BETS para que el chequeo temporal formal quede aguas arriba dentro de `pre_beast`.

### 7. `snakefile` actualizado

Archivo: `snakefile`

`rule all` ahora exige:

- `results/beast_pre/root_to_tip/treetime_clock.done`
- `results/beast_pre/bets/bets_runs.tsv`
- `results/beast_pre/bets/bets_bayes_factors.tsv`
- `results/beast_pre/model_test/panel_main_concat.iqtree`
- `results/beast_pre/model_test/panel_main_concat.best_model.txt`

### 8. Config actualizada para BETS

Archivo: `config/config.yml`

- `beast_pre` mantiene los parámetros de panel selection.
- Se añade `beast_pre.bets` con:
  - `threads`
  - `package_dir`
  - las cuatro rutas XML a ejecutar

### 9. Limpieza de la rama segment-aware

Se eliminaron del flujo principal:

- `BEAST_SEGMENT`
- toda la rama por segmento dentro de `02_pre_beast.smk`
- el script `code/02_Beast/build_single_segment_iqtree_partition.py`

## Archivos realmente afectados

- `rules/02_pre_beast.smk`
- `snakefile`
- `config/config.yml`
- `analysis/segment_pre_beast_review.md`
- `code/02_Beast/run_beast_bets.py`
- `code/02_Beast/summarize_bets.py`

Y se mantiene el cambio útil en:

- `code/02_Beast/subset_alignment_and_prune_tree.py`

## Qué ya no forma parte de esta versión

- cualquier lógica segment-aware para pre-BEAST
- generación de particiones IQ-TREE en `01_ml_trees`
- generación automática de XMLs para BETS
- BEAUti / XML / particionado final de BEAST fuera de ejecutar BETS ya preparado

## Validación esperada

El dry-run debería mostrar estos jobs pre-BEAST:

- `build_beast_panels`
- `subset_panel_concat_alignment_and_prune_tree`
- `build_treetime_dates`
- `run_root_to_tip`
- `ensure_beast_model_selection_package`
- `run_bets_scenario`
- `summarize_bets`
- `build_concat_iqtree_partition_nexus`
- `iqtree_model_test_concat_subset`

Y las entradas clave de IQ-TREE deberían ser:

- alignment:
  - `data/beast_pre/panels/panel_main_concat.subset.fasta`
- partition file:
  - `data/beast_pre/panels/panel_main_concat.subset.iqtree.nex`

## Nota

Este review describe la versión final pedida ahora:

- concatenado subseteado,
- root-to-tip sobre ese mismo dataset,
- BETS ejecutado desde XMLs externos ya preparados,
- un solo IQ-TREE particionado,
- sin generación automática de plantillas XML.
