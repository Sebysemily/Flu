configfile: "config/config.yml"
configfile: "config/paths.yml"

include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_pre_beast.smk"
include: "rules/03_beast.smk"

BY_SEGMENT_FASTAS = expand(
    f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta",
    segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
)

LINEAGE_ALIGNMENT_TARGETS = [
    *expand(f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft", segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]),
    f"{RESULTS_QC_METRICS}/nextclade/HA_eliminated.csv",
    f"{RESULTS_QC_METRICS}/nextclade/HA/nextclade_report.csv",

    *expand("figures/HA_PB2_lineage/{segment}_lineage_fast_tree.png", segment=["HA", "PB2"]),
    *expand("figures/HA_PB2_lineage/{segment}_lineage_fast_tree.rds", segment=["HA", "PB2"]),
    "figures/HA_PB2_lineage/HA_PB2_lineage_composite.png",
]

MAIN_PANEL_TREES = expand(
    f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/H5N1_{{segment}}_fast.treefile",
    segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
)

MAIN_PHYLOGENY_TARGETS = [
    *BY_SEGMENT_FASTAS,
    *LINEAGE_ALIGNMENT_TARGETS,
    MAIN_PANEL_METADATA,
    *MAIN_PANEL_TREES,
]

ALL_VALIDATION_TARGETS = [
    *MAIN_PHYLOGENY_TARGETS,
    FILTRADO_CSV,
    "figures/main_panel_8_segments_collapsed.png",
    "results/pre_beast/rtt_HA/treetime_clock.done",
    "metadata/beast/metadata_beast.tsv",
    "results/phylogeny/flu_mut/flumut_report_markers.tsv",
    "results/phylogeny/flu_mut/flumut_report_mutations.tsv",
    *BEAST_GSS_TARGETS
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
