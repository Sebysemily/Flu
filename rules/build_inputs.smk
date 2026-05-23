import glob
import os


MIRA_BASE = config["mira_base_dir"]
FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
ECUADOR_DATE_SOURCE = config.get("ecuador_date_source", "reception")
MIRA_GISAID_FASTA = "data/assembled/H5N1_EC_gisaid_from_mira.fasta"


def get_existing_mira_fastas():
    mira_fastas = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
    mira_fastas += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
    return sorted(set(mira_fastas))


MIRA_FASTAS = get_existing_mira_fastas()


def get_gisaid_input_fastas():
    patterns = ["*.fasta", "*.fa", "*.fna", "*.fas"]
    fastas = []
    for pattern in patterns:
        fastas.extend(glob.glob(os.path.join("data", "input", pattern)))
    return sorted(set(fastas))


GISAID_INPUT_FASTAS = get_gisaid_input_fastas()
STANDARD_HEADER_FASTAS = GISAID_INPUT_FASTAS if GISAID_INPUT_FASTAS else [MIRA_GISAID_FASTA]


rule build_gisaid_input_from_mira:
    input:
        MIRA_FASTAS,
        FILTRADO_CSV
    output:
        "data/assembled/ecuador_intermediate_raw_segments.csv",
        "data/assembled/ecuador_intermediate_summary.csv",
        "data/assembled/ecuador_intermediate_audit.csv",
        "data/assembled/ecuador_intermediate_issues.csv",
        MIRA_GISAID_FASTA
    shell:
        r"""
        python code/build_inputs/build_gisaid_input_from_mira.py
        """


rule build_standard_date_headers:
    input:
        gisaid_fastas=STANDARD_HEADER_FASTAS,
        metadata_csv=FILTRADO_CSV
    output:
        fasta="data/assembled/H5N1_EC.fasta",
        summary="data/assembled/H5N1_EC_summary.csv"
    params:
        ecuador_date_source=ECUADOR_DATE_SOURCE
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/build_standard_date_headers.py \
            --input-fasta {input.gisaid_fastas:q} \
            --metadata-csv {input.metadata_csv:q} \
            --ecuador-date-source {params.ecuador_date_source:q} \
            --output-fasta {output.fasta:q} \
            --summary-csv {output.summary:q}
        """


rule download_context_and_merge_denv2_fasta:
    input:
        ecuador_fasta="data/assembled/H5N1_EC.fasta",
        context_metadata_tsv=config.get("context_metadata_tsv", "config/final_metadata_50_per_country.tsv")
    output:
        context_fasta="data/assembled/H5N1_context.fasta",
        context_summary="data/assembled/H5N1_context_summary.csv",
        final_fasta="data/final/H5N1_final.fasta"
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/download_context_and_merge_denv2_fasta.py \
            --ecuador-fasta {input.ecuador_fasta:q} \
            --context-metadata-tsv {input.context_metadata_tsv:q} \
            --context-fasta-out {output.context_fasta:q} \
            --context-summary-out {output.context_summary:q} \
            --final-fasta-out {output.final_fasta:q}
        """
