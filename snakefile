configfile: "config/config.yml"
configfile: "config/paths.yml"

include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_pre_beast.smk"

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
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
        rep=range(1, int(config.get("iqtree_replicates", 5)) + 1)
    ),
    *expand(
        f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}_final.treefile",
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    ),
    f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram_ha_vs_pb2.png",
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
    MAIN_PANEL_METADATA,
    *expand(
        "results/figures/{segment}_tree.png",
        segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS", "subset_concat"]
    ),
]

ALL_VALIDATION_TARGETS = [
    *MAIN_PHYLOGENY_TARGETS,
    *RTT_AND_ML_SUBSET_TARGETS,
]

rule all:
    input:
        ALL_VALIDATION_TARGETS
