Pipeline del trabajo realizado para "" que consolida la filodinamica de influenza A en ecuador con contexto.

Requisitos
- snakemake
- Conexion a internet si se quiere descargar el contexto regional desde NCBI

Ejecucion
- Preparar inputs de Ecuador desde MIRA (stage 1, para construir desde raw):

snakemake --cores all --use-conda build_gisaid


- Correr todo el pipeline principal (stage 2, desde descarga de contexto y GISAID hasta arboles ML):

snakemake --cores all --use-conda 


- Correr la validacion codon+RF (opcional, para revision interna o reviewers):

snakemake --cores all --use-conda rf_validation

    - Validacion con particiones codon en 6 segmentos (PB2, PB1, PA, HA, NP, NA):
        - Alineacion concatenada codon-aware en results/validation/concat/
        - Arbol base y 5 replicas RF sin bootstrap en results/validation/raxml/full_concat/ y results/validation/rf/
        - Resumen de distancias RF en results/validation/rf/rf_summary/rf_summary.tsv (incluye columna run_label para trazabilidad)
        - Arboles per-segment codon sin bootstrap en results/validation/per_segment/{segment}/
    
    - Validacion RF simple (GTR+G sin particiones codon) para 2 segmentos restantes (NS, MP):
        - Arbol base y 5 replicas RF por segmento en results/validation/rf/
        - Resumen de distancias RF por segmento en results/validation/rf/rf_summary_{segment}/rf_summary.tsv

Flujo de reglas
Stage 1 (build_gisaid) — preparacion de inputs Ecuador desde MIRA

1. build_ecuador_intermediate_input
    - Copia los amended_consensus.fasta a data/all_amended_fasta
    - Selecciona el mejor segmento por muestra
    - Genera FASTA por muestra y archivos intermedios de Ecuador en data/assembled

2. build_h5n1_ec_fasta
    - Toma solo segmentos realmente ensamblados de Ecuador
    - Usa Provincia y Fecha recepcion de flu_filtrado
    - Genera el FASTA de Ecuador: data/input/H5N1_EC.fasta
    - Formato de encabezado: muestra/segmento/lugar/año
    - El lugar se normaliza a una forma canonica sin espacios ni guiones, por ejemplo SantaElena

Stage 2 (pipeline principal) — contexto, alineacion y filogenia

3. build_h5n1_final_fasta
    - Descarga secuencias contextuales por accession desde NCBI
    - Las renombra al mismo formato de influenza usado para Ecuador
    - Une Ecuador + contexto:
        - data/input/H5N1_context.fasta
        - data/input/H5N1_context_summary.csv
        - data/final/H5N1_final.fasta

4. split_h5n1_final_by_segment
    - Divide data/final/H5N1_final.fasta en 8 FASTA por segmento
    - Salida en data/phylogeny/by_segment/

5. mafft_align_per_segment
    - Alinea cada segmento con MAFFT
    - Salida por segmento en data/phylogeny/aligned/H5N1_{segment}.mafft

6. build_segment_codon_partitions
    - Genera particiones codon (cp12/cp3) para los 6 segmentos codificantes (PB2, PB1, PA, HA, NP, NA)
    - Salida en data/phylogeny/aligned/partitions/

7. raxml_ng_tree_per_segment_codon
    - Arboles ML en RAxML-NG para PB2, PB1, PA, HA, NP, NA bajo un modelo nucleotidico con
      particiones por posicion codon (posiciones 1+2 vs 3), con GTR+G asignado a cada particion
    - Busqueda ML + bootstrap, calcula soporte TBE

8. raxml_ng_tree_per_segment_simple
    - Arboles ML en RAxML-NG para NS y MP bajo un modelo GTR+G (sin particiones codon)
    - Busqueda ML + bootstrap, calcula soporte TBE

9. concat_aligned_segments_with_partitions
    - Concatena alineaciones por segmento con particiones para arbol full-concat
    - Salida: data/phylogeny/aligned/H5N1_full_concat_beast.mafft

10. raxml_ng_tree_full_concat
    - Arbol ML en RAxML-NG sobre el concatenado de los 8 segmentos bajo un modelo nucleotidico
      con particiones por posicion codon (posiciones 1+2 vs 3), con GTR+G asignado a cada particion
    - Salida: results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE

Estructura de salidas
- Intermedios de ensamblaje:
    - data/all_amended_fasta/
    - data/assembled/
    - data/assembled/ecuador_intermediate_sequences.fasta
    - data/assembled/ecuador_intermediate_summary.csv
    - data/assembled/ecuador_intermediate_audit.csv
    - data/assembled/ecuador_intermediate_issues.csv
    - data/assembled/ecuador_intermediate_per_sample/

- Build_inputs.smk:
    - inputs: 
        - data/input/H5N1_EC.fasta
        - data/input/H5N1_EC_summary.csv
        - data/input/H5N1_context.fasta
        - data/input/H5N1_context_summary.csv

    - Salida final:
        - data/final/H5N1_final.fasta

- 01_ml_trees:
    -inputs:
        - data/phylogeny/by_segment/
        - data/phylogeny/by_segment_summary.csv
        - data/phylogeny/aligned/H5N1_{segment}.mafft
        - data/phylogeny/aligned/H5N1_full_concat_beast.mafft
    -outputs:
        - results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE
        - results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE

- 00_validation_codon_rf:
    - outputs:
        - results/validation/concat/
        - results/validation/raxml/full_concat/
        - results/validation/rf/rf_summary/rf_summary.tsv (6 segmentos con particiones codon)
        - results/validation/rf/rf_summary_NS/rf_summary.tsv (NS simple RF)
        - results/validation/rf/rf_summary_MP/rf_summary.tsv (MP simple RF)
        - results/validation/per_segment/{segment}/H5N1_{segment}_codon_validation.raxml.bestTreeCollapsed

Graficos del workflow (DAG, rulegraph y filegraph) guardados en analysis/.

./scripts/gen_snakemake_graphs.sh all analysis


Notas
- El FASTA final contiene primero Ecuador y despues el contexto regional.
- Ejemplo de encabezado: Flu-0008/NS/SantaElena/2023
- En contexto, el nombre de muestra se fuerza a ser unico anexando accession cuando existe isolate (isolate_accession) y usando accession cuando isolate falta.
- La normalizacion del lugar se aplica igual a Ecuador y al contexto para evitar conteos dobles por variantes como Santa Elena, Santa_Elena o Santa-Elena.
- La descarga contextual se hace en lotes desde NCBI para reducir tiempo de ejecucion.
- opciones de cambio en config/config.yml
    mira_base_dir: ".."
    flu_filtrado: "config/flu_filtrado.csv"
    context_metadata_tsv: "config/final_metadata_50_per_country_isolates.tsv"
    ecuador_date_source: "collection"
    random_seed: 39809473
    max_threads: 20
-Preparacion de runs MIRA
    1. Instalar y levantar el contenedor de MIRA.
    2. Colocar cada run dentro de ~/MIRA_NGS/runX.
    3. Asegurar que cada run tenga samplesheet.csv.
    4. Ejecutar MIRA para cada run:

    cd ~/MIRA_NGS
    for run in run*; do
            if [ -d "$run" ] && [ -f "$run/samplesheet.csv" ]; then
                    echo "=== Reanudando $run ==="
                    docker exec -w /data mira MIRA.sh \
                            -s "$run/samplesheet.csv" \
                            -r "$run" \
                            -e Flu-ONT \
                            > "$run/mira_batch_resume.out" \
                            2> "$run/mira_batch_resume.err"
            fi
    done
    Corrido originalmente con MIRA v2.0.0.
