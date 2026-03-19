Pipeline de ensamblaje y consolidacion de H5N1 a partir de salidas de MIRA.

Requisitos
- Python 3
- pandas
- snakemake
- Conexion a internet si se quiere descargar el contexto regional desde NCBI

Preparacion de runs MIRA
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

Configuracion
- Metadata local de muestras Ecuador: config/flu_filtrado.csv
- Metadata contextual regional: config/final_metadata_50_per_country.tsv
- Ruta base de runs MIRA: config/config.yml

Ejecucion
- Correr todo el pipeline (incluye alineamientos MAFFT y arboles ML por segmento):

snakemake --cores all --use-conda all

- Reintentar solo incompletos/fallidos:

snakemake --cores all --use-conda --rerun-incomplete all

- Forzar recomputacion completa desde cero:

snakemake --cores all --use-conda --forceall all

- Construir solo el FASTA final combinado (sin filogenia):

snakemake --cores 1 --use-conda data/final/H5N1_final.fasta

Flujo de reglas
1. build_ecuador_intermediate_input
     - Copia los amended_consensus.fasta a data/all_amended_fasta
     - Selecciona el mejor segmento por muestra
    - Genera FASTA por muestra y archivos intermedios de Ecuador en data/assembled

2. build_h5n1_ec_fasta
    - Toma solo segmentos realmente ensamblados de Ecuador
     - Usa Provincia y Fecha recepcion de flu_filtrado
    - Genera el FASTA exportable de Ecuador: data/exports/H5N1_EC.fasta
    - Formato de encabezado: muestra/segmento/lugar/año
    - El lugar se normaliza a una forma canonica sin espacios ni guiones, por ejemplo SantaElena

3. build_h5n1_final_fasta
     - Descarga secuencias contextuales por accession desde NCBI
    - Las renombra al mismo formato de influenza usado para Ecuador
    - Une Ecuador + contexto en un unico FASTA final

4. split_h5n1_final_by_segment
    - Divide data/final/H5N1_final.fasta en 8 FASTA (uno por segmento)

5. mafft_align_per_segment
    - Alinea cada segmento con MAFFT
    - Salida por segmento en data/phylogeny/aligned/H5N1_{segment}.mafft

6. raxml_ng_tree_per_segment
    - Ejecuta RAxML-NG con modelo GTR+G
    - Corre busqueda ML + bootstrap y calcula soporte FBP/TBE

7. raxml_trees_all_segments
    - Agregador final que exige soporte FBP para los 8 segmentos

Estructura de salidas
- Intermedios de ensamblaje:
    - data/all_amended_fasta/
    - data/assembled/
    - data/assembled/ecuador_intermediate_sequences.fasta
    - data/assembled/ecuador_intermediate_summary.csv
    - data/assembled/ecuador_intermediate_audit.csv
    - data/assembled/ecuador_intermediate_issues.csv
    - data/assembled/ecuador_intermediate_per_sample/

- Exportables para analisis:
    - data/exports/H5N1_EC.fasta
    - data/exports/H5N1_EC_summary.csv
    - data/exports/H5N1_context.fasta
    - data/exports/H5N1_context_summary.csv

- Salida final:
    - data/final/H5N1_final.fasta

- Filogenia por segmento:
    - data/phylogeny/by_segment/
    - data/phylogeny/by_segment_summary.csv
    - data/phylogeny/aligned/H5N1_*.mafft
    - results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.bestTreeCollapsed
    - results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportFBP
    - results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE

Graficos del workflow (DAG, rulegraph y filegraph)
- Generar graficos hasta la ultima regla (arboles ML):

./scripts/gen_snakemake_graphs.sh all analysis

- Archivo clave para revisar dependencias de archivos:
    - analysis/pipeline_filegraph.svg

- Si quieres solo hasta el FASTA final (sin filogenia):

./scripts/gen_snakemake_graphs.sh data/final/H5N1_final.fasta analysis

Notas
- El FASTA final contiene primero Ecuador y despues el contexto regional.
- Ejemplo de encabezado: Flu-0008/NS/SantaElena/2023
- En contexto, el nombre de muestra se fuerza a ser unico anexando accession cuando existe isolate (isolate_accession) y usando accession cuando isolate falta.
- La normalizacion del lugar se aplica igual a Ecuador y al contexto para evitar conteos dobles por variantes como Santa Elena, Santa_Elena o Santa-Elena.
- La descarga contextual se hace en lotes desde NCBI para reducir tiempo de ejecucion.
- Si cambia el archivo contextual, actualizar config/final_metadata_50_per_country.tsv o la ruta en config/config.yml.

Corrido originalmente con MIRA v2.0.0.
