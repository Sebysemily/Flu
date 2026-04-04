FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
ECUADOR_DATE_SOURCE = config.get("ecuador_date_source", "reception")


rule build_h5n1_final_fasta:
    input:
        ecuador_fasta="data/input/H5N1_EC.fasta",
        context_metadata_tsv=config.get("context_metadata_tsv", "config/final_metadata_50_per_country.tsv")
    output:
        context_fasta="data/input/H5N1_context.fasta",
        context_summary="data/input/H5N1_context_summary.csv",
        final_fasta="data/final/H5N1_final.fasta"
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/download_context_and_merge_denv2_fasta.py \
            --ecuador-fasta {input.ecuador_fasta} \
            --context-metadata-tsv {input.context_metadata_tsv} \
            --context-fasta-out {output.context_fasta} \
            --context-summary-out {output.context_summary} \
            --final-fasta-out {output.final_fasta}
        """


rule build_h5n1_final_beast_fasta:
    input:
        final_fasta="data/final/H5N1_final.fasta",
        metadata_csv=FILTRADO_CSV,
        context_metadata_tsv=config.get("context_metadata_tsv", "config/final_metadata_50_per_country.tsv")
    output:
        beast_fasta="data/final/H5N1_final_beast.fasta",
        beast_summary="data/final/H5N1_final_beast_summary.csv"
    params:
        ecuador_date_source=ECUADOR_DATE_SOURCE
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/build_h5n1_beast_fasta.py \
            --final-fasta {input.final_fasta} \
            --flu-filtrado-csv {input.metadata_csv} \
            --context-metadata-tsv {input.context_metadata_tsv} \
            --ecuador-date-source {params.ecuador_date_source} \
            --output-fasta {output.beast_fasta} \
            --output-summary {output.beast_summary}
        """