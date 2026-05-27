import glob
import os

FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/standard_header_input_fasta")

# =====================================================================
# Rule: build_gisaid_input_from_mira
# =====================================================================

MIRA_BASE = config.get("mira_base_dir", "../../../MIRA_NGS")

def get_existing_mira_fastas():
    mira_fastas = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
    mira_fastas += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
    return sorted(set(mira_fastas))

MIRA_FASTAS = get_existing_mira_fastas()
DATA_INPUT = config.get("data_input", "data/input")
MIRA_GISAID_FASTA = f"{DATA_INPUT}/H5N1_EC_gisaid_from_mira.fasta"

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

# =====================================================================
# Rule: build_standard_date_headers
# =====================================================================

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

# =====================================================================
# Rule: build_standard_headers_for_gisaid_context
# =====================================================================

CONTEXT_FASTA_RAW = config.get("context_fasta_raw", "data/context/gisaid_epiflu_sequence.fasta")

rule build_standard_headers_for_gisaid_context:
    input:
        ecuador_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_EC.fasta",
        context_fasta=CONTEXT_FASTA_RAW,
        metadata_csv=FILTRADO_CSV
    output:
        context_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_context.fasta",
        final_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final.fasta"
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/process_gisaid_context.py \
            --ecuador-fasta {input.ecuador_fasta:q} \
            --context-fasta-in {input.context_fasta:q} \
            --context-fasta-out {output.context_fasta:q} \
            --final-fasta-out {output.final_fasta:q} \
            --metadata-csv {input.metadata_csv:q}
        """
