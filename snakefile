configfile: "config/config.yml"

include: "rules/build_gisaid_input_from_mira.smk"
include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_pre_beast.smk"
include: "rules/03_beast.smk"

BUILD_GISAID_TARGETS = [
    "data/all_amended_fasta",
    "data/assembled/ecuador_intermediate_summary.csv",
    "data/assembled/ecuador_intermediate_audit.csv",
    "data/assembled/ecuador_intermediate_issues.csv",
    "data/assembled/ecuador_intermediate_sequences.fasta",
    "data/assembled/ecuador_intermediate_per_sample",
    "data/input/H5N1_EC.fasta",
    "data/input/H5N1_EC_summary.csv",
]

INPUT_CONTEXT_FASTA = "data/input/H5N1_context.fasta"
INPUT_CONTEXT_SUMMARY = "data/input/H5N1_context_summary.csv"
FINAL_FASTA = "data/final/H5N1_final.fasta"
BUILD_INPUTS_TARGETS = [
    INPUT_CONTEXT_FASTA,
    INPUT_CONTEXT_SUMMARY,
    FINAL_FASTA,
]

BY_SEGMENT_FASTAS = expand(
    "data/phylogeny/by_segment/H5N1_{segment}.fasta",
    segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
)
BY_SEGMENT_SUMMARY = "data/phylogeny/by_segment_summary.csv"
SEGMENT_TREE_SUPPORTS = expand(
    "results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
    segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
)
MAIN_PHYLOGENY_TARGETS = [
    *BY_SEGMENT_FASTAS,
    BY_SEGMENT_SUMMARY,
    FULL_CONCAT_ALIGNMENT,
    FULL_CONCAT_PREFIX + ".raxml.supportTBE",
    *SEGMENT_TREE_SUPPORTS,
]

PRE_BEAST_TARGETS = [
    PRE_BEAST_SOURCE_QC_METRICS,
    PRE_BEAST_SOURCE_QC_OUTLIERS,
    PRE_BEAST_SOURCE_QC_SUMMARY,
    PRE_BEAST_SOURCE_QC_REPORT,
    BEAST_FILTERED_SUBSET_ALIGNMENT,
    BEAST_FILTERED_SUBSET_TREE,
    BEAST_FILTERED_SUBSET_AUDIT,
    PRE_BEAST_DATES,
    *PRE_BEAST_SEGMENT_QC_METRICS,
    *PRE_BEAST_SEGMENT_QC_OUTLIERS,
    *PRE_BEAST_SEGMENT_QC_SUMMARIES,
    *PRE_BEAST_SEGMENT_QC_REPORTS,
    BEAST_FINAL_SUBSET_ALIGNMENT,
    BEAST_FINAL_SUBSET_AUDIT,
    PRE_BEAST_FINAL_SEGMENT_QC_METRICS,
    PRE_BEAST_FINAL_SEGMENT_QC_OUTLIERS,
    PRE_BEAST_FINAL_SEGMENT_QC_SUMMARY,
    PRE_BEAST_FINAL_SEGMENT_QC_REPORT,
    PRE_BEAST_PANEL_EXCLUSIONS,
    PRE_BEAST_PANEL_EXCLUSIONS_SUMMARY,
    PRE_BEAST_RTT_EXCLUSIONS,
    PRE_BEAST_RTT_EXCLUSIONS_SUMMARY,
    *PRE_BEAST_RAW_FINAL_SEGMENT_FASTAS,
    *PRE_BEAST_FINAL_SEGMENT_FASTAS,
    PRE_BEAST_ROOT_TO_TIP_DONE,
    BEAST_FINAL_DATES,
    *PREPARED_BEAST_XMLS,
]

PAPER_FIGURE_TARGETS = []

ALL_VALIDATION_TARGETS = [
    *BUILD_INPUTS_TARGETS,
    *MAIN_PHYLOGENY_TARGETS,
    *PRE_BEAST_TARGETS,
    *BEAST_PUBLIC_TARGETS,
]


rule all:
    input:
        *BEAST_PUBLIC_TARGETS


rule build_gisaid_outputs:
    input:
        BUILD_GISAID_TARGETS


rule build_inputs_outputs:
    input:
        BUILD_INPUTS_TARGETS


rule main_phylogeny:
    input:
        MAIN_PHYLOGENY_TARGETS


rule pre_beast_outputs:
    input:
        PRE_BEAST_TARGETS


rule beast_runs:
    input:
        BEAST_PUBLIC_TARGETS


rule paper_figures:
    input:
        PAPER_FIGURE_TARGETS


rule all_validation_outputs:
    input:
        ALL_VALIDATION_TARGETS


rule build_gisaid:
    input:
        BUILD_GISAID_TARGETS
