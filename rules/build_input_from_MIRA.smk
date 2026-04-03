import os
import glob

MIRA_BASE = config["mira_base_dir"]

MIRA_FASTAS = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
MIRA_FASTAS += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
MIRA_FASTAS = sorted(set(MIRA_FASTAS))
 
FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
ECUADOR_DATE_SOURCE = config.get("ecuador_date_source", "reception")

rule build_ecuador_intermediate_input:
    input:
        MIRA_FASTAS,
        FILTRADO_CSV
    output:
        directory("data/all_amended_fasta"),
        "data/assembled/ecuador_intermediate_raw_segments.csv",
        "data/assembled/ecuador_intermediate_expected_lengths.csv",
        "data/assembled/ecuador_intermediate_summary.csv",
        "data/assembled/ecuador_intermediate_audit.csv",
        "data/assembled/ecuador_intermediate_issues.csv",
        "data/assembled/ecuador_intermediate_sequences.fasta",
        directory("data/assembled/ecuador_intermediate_per_sample")
    shell:
        r"""
        python code/build_input_from_MIRA/build_gisaid_input_from_mira.py
        """


rule build_h5n1_ec_fasta:
    input:
        per_sample_fastas="data/assembled/ecuador_intermediate_per_sample",
        audit_csv="data/assembled/ecuador_intermediate_audit.csv",
        metadata_csv=FILTRADO_CSV
    output:
        fasta="data/exports/H5N1_EC.fasta",
        summary="data/exports/H5N1_EC_summary.csv"
    params:
        ecuador_date_source=ECUADOR_DATE_SOURCE
    shell:
        r"""
        python code/build_input_from_MIRA/build_denv2_style_fasta_from_assembled.py \
            --per-sample-dir {input.per_sample_fastas} \
            --audit-csv {input.audit_csv} \
            --metadata-csv {input.metadata_csv} \
            --ecuador-date-source {params.ecuador_date_source} \
            --output-fasta {output.fasta} \
            --summary-csv {output.summary}
        """


rule build_h5n1_final_fasta:
    input:
        ecuador_fasta="data/exports/H5N1_EC.fasta",
        context_metadata_tsv=config.get("context_metadata_tsv", "config/final_metadata_50_per_country.tsv")
    output:
        context_fasta="data/exports/H5N1_context.fasta",
        context_summary="data/exports/H5N1_context_summary.csv",
        final_fasta="data/final/H5N1_final.fasta"
    shell:
        r"""
        python code/build_input_from_MIRA/download_context_and_merge_denv2_fasta.py \
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
        python code/build_input_from_MIRA/build_h5n1_beast_fasta.py \
            --final-fasta {input.final_fasta} \
            --flu-filtrado-csv {input.metadata_csv} \
            --context-metadata-tsv {input.context_metadata_tsv} \
            --ecuador-date-source {params.ecuador_date_source} \
            --output-fasta {output.beast_fasta} \
            --output-summary {output.beast_summary}
        """
