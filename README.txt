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
- Correr todo el pipeline:

snakemake -j 1

- O construir solo el FASTA final combinado:

snakemake data/final/H5N1_final.fasta -j 1

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

Notas
- El FASTA final contiene primero Ecuador y despues el contexto regional.
- Ejemplo de encabezado: Flu-0008/NS/SantaElena/2023
- La normalizacion del lugar se aplica igual a Ecuador y al contexto para evitar conteos dobles por variantes como Santa Elena, Santa_Elena o Santa-Elena.
- La descarga contextual se hace en lotes desde NCBI para reducir tiempo de ejecucion.
- Si cambia el archivo contextual, actualizar config/final_metadata_50_per_country.tsv o la ruta en config/config.yml.

Corrido originalmente con MIRA v2.0.0.
