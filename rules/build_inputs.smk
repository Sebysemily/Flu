import glob
import os


MIRA_BASE = config["mira_base_dir"]
FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
MIRA_GISAID_FASTA = "data/input/H5N1_EC_gisaid_from_mira.fasta"


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
    # Exclude the generated MIRA fasta to keep the same fallback logic
    mira_gen = os.path.abspath(os.path.join("data", "input", "H5N1_EC_gisaid_from_mira.fasta"))
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
        fasta="data/standard_header_input_fasta/H5N1_EC.fasta"
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
        ecuador_fasta="data/standard_header_input_fasta/H5N1_EC.fasta",
        context_fasta="data/context/gisaid_epiflu_sequence.fasta",
        metadata_csv=FILTRADO_CSV
    output:
        context_fasta="data/standard_header_input_fasta/H5N1_context.fasta",
        final_fasta="data/standard_header_input_fasta/H5N1_final.fasta"
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
