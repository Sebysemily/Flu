import glob
import os

FILTRADO_CSV = config.get("flu_filtrado", "metadata/flu_filtrado.csv")
DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/standard_header_input_fasta")
DATA_PHYLOGENY = config.get("data_phylogeny", "data/phylogeny")
PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
MAIN_PANEL_METADATA = "metadata/H5N1_context.csv"
CONTEXT_BASE_METADATA = "metadata/context_base.csv"
CONTEXT_DIR = config.get("context_dir", "data/context")
CONTEXT_GENOFLU_RUN = config.get("context_genoflu_run", f"{CONTEXT_DIR}/genoflu_run")
GENOFLU_DIR = "resources/GenoFLU-multi"
GENOFLU_CONDA = "../envs/01_ml_trees.yml"

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
            --output-dir {params.output_dir}
        """

# =====================================================================
# Rule: genoflu_multi (All sequences post-QC)
# =====================================================================

rule genoflu_multi:
    input:
        input_fastas=expand(f"{DATA_PHYLOGENY}/by_segment_qc_filtered/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS)
    output:
        results="metadata/genoflu_results.tsv"
    params:
        fasta_dir=f"{DATA_PHYLOGENY}/by_segment_qc_filtered",
        genoflu_dir=GENOFLU_DIR,
    conda:
        GENOFLU_CONDA
    shell:
        r"""
        if [ -n "$(ls -A {params.fasta_dir}/*.fasta 2>/dev/null)" ]; then
            mkdir -p {params.fasta_dir}/results
            python {params.genoflu_dir}/bin/genoflu-multi.py -f {params.fasta_dir} -m --run_incomplete
            if [ -f {params.fasta_dir}/results/results.tsv ]; then
                cp {params.fasta_dir}/results/results.tsv {output.results}
            else
                touch {output.results}
            fi
        else
            touch {output.results}
        fi
        rm -rf "{params.fasta_dir}/results" "{params.fasta_dir}/temp" "{params.fasta_dir}/blast"
        """

# =====================================================================
# Rule: genoflu_results_to_metadata
# =====================================================================
rule genoflu_results_to_metadata:
    input:
        results="metadata/genoflu_results.tsv"
    output:
        done="metadata/genoflu_results.done",
    params:
        metadata_csv=FILTRADO_CSV
    conda:
        GENOFLU_CONDA
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/genoFlu_to_metadata.py \
            --genoflu-results {input.results} \
            --metadata-csv {params.metadata_csv}
        rm -rf "data/input/results" "data/input/temp" "data/input/blast"
        touch {output.done}
        """

# Removed prepare_context_genoflu_fastas and genoflu_context_multi
# Since genoflu_multi now runs on all QC-filtered sequences.

# =====================================================================
# Rule: build_h5n1_context_metadata
# =====================================================================
rule build_h5n1_context_metadata:
    input:
        flu_filtrado=FILTRADO_CSV,
        context_fasta=CONTEXT_FASTA_RAW,
        genoflu_context="metadata/genoflu_results.tsv",
        genoflu_ecuador_done="metadata/genoflu_results.done",
        human_fasta=os.path.join(CONTEXT_DIR, "Human_host.fasta"),
        avian_fasta=os.path.join(CONTEXT_DIR, "avian_only_regional_context.fasta"),
    output:
        context_base=CONTEXT_BASE_METADATA,
        metadata=MAIN_PANEL_METADATA,
    params:
        ecuador_date_source="collection",
    conda:
        GENOFLU_CONDA
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_inputs/build_h5n1_context_metadata.py \
            --context-fasta {input.context_fasta} \
            --flu-filtrado {input.flu_filtrado} \
            --genoflu-context-results {input.genoflu_context} \
            --ecuador-date-source {params.ecuador_date_source} \
            --context-base-out {output.context_base} \
            --metadata-out {output.metadata} \
            --human-fasta {input.human_fasta} \
            --avian-fasta {input.avian_fasta}
        """

# =====================================================================
# Rule: genoflu_plot
# =====================================================================
rule plot_flu_lineage:
    input:
        metadata=MAIN_PANEL_METADATA,
    output:
        plot="figures/build_inputs/flu_lineage.png",
        rds="figures/build_inputs/flu_lineage.rds",
    params:
        script="code/segment_analysis/flu_lineage_heatmap.R"
    conda:
        "../envs/r.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.plot})
        Rscript {params.script} {input.metadata} {output.plot}
        test -s {output.rds}
        """
