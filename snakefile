configfile: "config/config.yml"
configfile: "config/paths.yml"

include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_pre_beast.smk"

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
MAIN_PHYLOGENY_TARGETS = [
    *BY_SEGMENT_FASTAS,
    f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat.mafft",
    f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat.treefile",
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
        f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rf_dist_matrix.rfdist",
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS", "concat_except_HA"]
    ),
    f"{RESULTS_PHYLOGENY}/segment_analysis/au_test_matrix.csv",
    f"{RESULTS_PHYLOGENY}/segment_analysis/rf_matrix.csv",
    f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram_ha_vs_na.png",
    f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram_ha_vs_concat_except_ha.png",
    f"{RESULTS_PHYLOGENY}/itol_styles/tree_colors.txt",
    f"{RESULTS_PHYLOGENY}/itol_styles/dataset_colorstrip.txt",
]

RTT_AND_ML_SUBSET_TARGETS = [
    MAIN_PANEL_FILTERED_ALIGNMENT,
    MAIN_PANEL_FILTERED_TREE,
    MAIN_PANEL_FILTERED_AUDIT,
    PRE_BEAST_DATES,
    MAIN_PANEL_ALIGNMENT,
    MAIN_PANEL_PARTITIONS,
    PRE_BEAST_ROOT_TO_TIP_DONE,
    *expand(
        f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/rep{{rep}}.treefile",
        rep=range(1, int(config.get("iqtree_replicates", 5)) + 1)
    ),
    f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/subset_concat_final.treefile",
    f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/rf_dist_matrix.rfdist",
    f"{RESULTS_PHYLOGENY}/itol_styles/subset_concat/tree_colors.txt",
    f"{RESULTS_PHYLOGENY}/itol_styles/subset_concat/dataset_colorstrip.txt",
]

ALL_VALIDATION_TARGETS = [
    *BUILD_INPUTS_TARGETS,
    *MAIN_PHYLOGENY_TARGETS,
    *RTT_AND_ML_SUBSET_TARGETS,
]


rule all:
    input:
        ALL_VALIDATION_TARGETS
