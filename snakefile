configfile: "config/config.yml"
configfile: "config/paths.yml"

include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_pre_beast.smk"

BY_SEGMENT_FASTAS = expand(
    f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta",
    segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
)

LINEAGE_PB2_HA_TARGETS = [
    f"{PROCESSED_ALIGNMENTS_HA_PB2_LINEAGE}/H5N1_PB2.mafft",
    f"{PROCESSED_ALIGNMENTS_HA_PB2_LINEAGE}/H5N1_HA.mafft",
    f"{RESULTS_QC_METRICS}/alignments/pb2_ha_lineage_filter_audit.csv",
    f"{RESULTS_PHYLOGENY}/iq-tree/HA/lineage/H5N1_HA_fast.treefile",
    f"{RESULTS_PHYLOGENY}/iq-tree/PB2/lineage/H5N1_PB2_fast.treefile",
    *expand("figures/HA_PB2_lineage/{segment}_lineage_fast_tree.png", segment=["HA", "PB2"]),
    *expand("figures/HA_PB2_lineage/{segment}_lineage_fast_tree.rds", segment=["HA", "PB2"]),
    "figures/HA_PB2_lineage/HA_PB2_lineage_composite.png",
]

MAIN_PHYLOGENY_TARGETS = [
    *BY_SEGMENT_FASTAS,
    *LINEAGE_PB2_HA_TARGETS,
    MAIN_PANEL_METADATA,
]

ALL_VALIDATION_TARGETS = [
    *MAIN_PHYLOGENY_TARGETS,
    FILTRADO_CSV,
]

GENOFLU_METADATA_TARGETS = [
    MAIN_PANEL_METADATA,
    "metadata/genoflu_results.done",
    "metadata/genoflu_context_results.tsv",
]

rule all:
    input:
        ALL_VALIDATION_TARGETS

rule genoflu_metadata:
    input:
        GENOFLU_METADATA_TARGETS
