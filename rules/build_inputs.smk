import glob
import os

FILTRADO_CSV = config.get("flu_filtrado", "metadata/flu_filtrado.csv")
DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/standard_header_input_fasta")
DATA_PHYLOGENY = config.get("data_phylogeny", "data/phylogeny")
PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
MAIN_PANEL_METADATA = "metadata/H5N1_context.csv"

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
        ancient(FILTRADO_CSV)
    output:
        MIRA_GISAID_FASTA
    shell:
        r"""
        python code/build_inputs/build_gisaid_input_from_mira.py
        """

# =====================================================================
# Rule: process_raw_to_segments
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
CONTEXT_FASTA_RAW = config.get("context_fasta_raw", "data/context/gisaid_epiflu_sequence.fasta")

rule process_raw_to_segments:
    input:
        ecuador_fastas=STANDARD_HEADER_FASTAS,
        context_fasta=CONTEXT_FASTA_RAW,
        metadata_csv=ancient(FILTRADO_CSV)
    output:
        segment_fastas=expand(f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS),
        metadata_out=MAIN_PANEL_METADATA
    params:
        output_dir=f"{DATA_PHYLOGENY}/by_segment",
        ecuador_date_source="collection"
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/process_raw_to_segments.py \
            --ecuador-fastas {input.ecuador_fastas} \
            --context-fasta {input.context_fasta} \
            --metadata-csv {input.metadata_csv} \
            --ecuador-date-source {params.ecuador_date_source} \
            --output-dir {params.output_dir} \
            --metadata-out {output.metadata_out}
        """

# =====================================================================
# Rule: genoflu_multi
# =====================================================================
INPUT_FASTAS = sorted(glob.glob(os.path.join(DATA_INPUT, "*.fasta")))

rule genoflu_multi:
    input:
        input_fastas=INPUT_FASTAS
    output:
        results="metadata/genoflu_results.tsv"
    params:
        fasta_dir=DATA_INPUT,
        genoflu_dir="resources/GenoFLU-multi"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python {params.genoflu_dir}/bin/genoflu-multi.py -f {params.fasta_dir} -m
        cp {params.fasta_dir}/results/results.tsv {output.results}
        """

# =====================================================================
# Rule: genoflu_results_to_metadata
# =====================================================================
rule genoflu_results_to_metadata:
    input:
        results="metadata/genoflu_results.tsv"
    output:
        metadata_csv=FILTRADO_CSV
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/genoFlu_to_metadata.py \
            --genoflu-results {input.results} \
            --metadata-csv {output.metadata_csv}
        """
