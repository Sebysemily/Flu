import glob
import os

# Organize paths dynamically from config (paths.yml)
MIRA_BASE = config.get("mira_base_dir", "../../../MIRA_NGS")
FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
DATA_INPUT = config.get("data_input", "data/input")
CONTEXT_FASTA_RAW = config.get("context_fasta_raw", "data/context/gisaid_epiflu_sequence.fasta")
DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/standard_header_input_fasta")
RESULTS_QC_METRICS = config.get("results_qc_metrics", "results/qc_metrics")

MIRA_GISAID_FASTA = f"{DATA_INPUT}/H5N1_EC_gisaid_from_mira.fasta"


def get_existing_mira_fastas():
    mira_fastas = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
    mira_fastas += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
    return sorted(set(mira_fastas))


MIRA_FASTAS = get_existing_mira_fastas()


def get_gisaid_input_fastas():
    patterns = ["*.fasta", "*.fa", "*.fna", "*.fas"]
    fastas = []
    for pattern in patterns:
        fastas.extend(glob.glob(os.path.join(DATA_INPUT, pattern)))
    # Exclude the generated MIRA fasta to keep the same fallback logic
    mira_gen = os.path.abspath(MIRA_GISAID_FASTA)
    fastas = [f for f in fastas if os.path.abspath(f) != mira_gen]
    return sorted(set(fastas))


GISAID_INPUT_FASTAS = get_gisaid_input_fastas()
STANDARD_HEADER_FASTAS = GISAID_INPUT_FASTAS if GISAID_INPUT_FASTAS else [MIRA_GISAID_FASTA]


rule build_gisaid_input_from_mira:
    input:
        MIRA_FASTAS,
        FILTRADO_CSV
    output:
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
        fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_EC.fasta"
    params:
        ecuador_date_source="collection"
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/build_standard_date_headers.py \
            --input-fasta {input.gisaid_fastas:q} \
            --metadata-csv {input.metadata_csv:q} \
            --ecuador-date-source {params.ecuador_date_source:q} \
            --output-fasta {output.fasta:q}
        """


rule build_standard_headers_for_gisaid_context:
    input:
        ecuador_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_EC.fasta",
        context_fasta=CONTEXT_FASTA_RAW,
        metadata_csv=FILTRADO_CSV
    output:
        context_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_context.fasta",
        final_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final.fasta"
    params:
        max_per_country_month=config.get("max_per_country_month", 2)
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/process_gisaid_context.py \
            --ecuador-fasta {input.ecuador_fasta:q} \
            --context-fasta-in {input.context_fasta:q} \
            --context-fasta-out {output.context_fasta:q} \
            --final-fasta-out {output.final_fasta:q} \
            --metadata-csv {input.metadata_csv:q} \
            --max-per-country-month {params.max_per_country_month}
        """
