configfile: "config/config.yml"
configfile: "config/paths.yml"

include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_pre_beast.smk"
include: "rules/03_beast.smk"

BUILD_GISAID_TARGETS = [
    f"{DATA_COMBINED_CONTEXT_EC}/H5N1_EC.fasta",
]
if not GISAID_INPUT_FASTAS:
    BUILD_GISAID_TARGETS = [
        MIRA_GISAID_FASTA,
        *BUILD_GISAID_TARGETS,
    ]

INPUT_CONTEXT_FASTA = f"{DATA_COMBINED_CONTEXT_EC}/H5N1_context.fasta"
FINAL_FASTA = f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final.fasta"
BUILD_INPUTS_TARGETS = [
    INPUT_CONTEXT_FASTA,
    FINAL_FASTA,
]

BY_SEGMENT_FASTAS = expand(
    f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta",
    segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
)
BY_SEGMENT_SUMMARY = f"{DATA_PHYLOGENY}/by_segment_summary.csv"
MAIN_PHYLOGENY_TARGETS = [
    *BY_SEGMENT_FASTAS,
    BY_SEGMENT_SUMMARY,
    FULL_CONCAT_ALIGNMENT,
    FULL_CONCAT_PREFIX + ".treefile",
    *expand(
        f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile",
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS", "concat_except_HA"],
        rep=range(1, int(config.get("iqtree_replicates", 5)) + 1)
    ),
    *expand(
        f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}_final.treefile",
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    ),
    *expand(
        f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rf_dist_matrix.rf",
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS", "concat_except_HA"]
    ),
    "results/phylogeny/concordance/full_concat.gcs",
    "results/phylogeny/concordance/full_concat.scf",
    f"{RESULTS_PHYLOGENY}/itol_styles/tree_colors.txt",
    f"{RESULTS_PHYLOGENY}/itol_styles/dataset_colorstrip.txt",
]

PRE_BEAST_TARGETS = [
    BEAST_FILTERED_SUBSET_ALIGNMENT,
    BEAST_FILTERED_SUBSET_TREE,
    BEAST_FILTERED_SUBSET_AUDIT,
    PRE_BEAST_DATES,
    FINAL_ALIGNMENT,
    FINAL_PARTITIONS,
    BAYESIAN_ALIGNMENT,
    PRE_BEAST_RTT_EXCLUSIONS,
    *PRE_BEAST_FINAL_SEGMENT_FASTAS,
    PRE_BEAST_ROOT_TO_TIP_DONE,
    BEAST_FINAL_DATES,
    *PREPARED_BEAST_XMLS,
    f"{RESULTS_PHYLOGENY}/iq-tree/subset_fast/panel_main_concat.final.treefile",
    *expand(
        f"{RESULTS_PHYLOGENY}/itol_styles/{{subdir}}/tree_colors.txt",
        subdir=PRE_BEAST_SEGMENTS + ["subset_fast"]
    ),
    *expand(
        f"{RESULTS_PHYLOGENY}/itol_styles/{{subdir}}/dataset_colorstrip.txt",
        subdir=PRE_BEAST_SEGMENTS + ["subset_fast"]
    ),
]

PAPER_FIGURE_TARGETS = []

ALL_VALIDATION_TARGETS = [
    *BUILD_INPUTS_TARGETS,
    *MAIN_PHYLOGENY_TARGETS,
    *PRE_BEAST_TARGETS,
]


rule all:
    input:
        ALL_VALIDATION_TARGETS


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
