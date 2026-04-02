configfile: "config/config.yml"

include: "rules/build_input_from_MIRA.smk"
include: "rules/01_ml_trees.smk"
include: "rules/02_Beast.smk"

rule all:
    input:
        "data/all_amended_fasta",
        "data/assembled/ecuador_intermediate_summary.csv",
        "data/assembled/ecuador_intermediate_audit.csv",
        "data/assembled/ecuador_intermediate_issues.csv",
        "data/assembled/ecuador_intermediate_sequences.fasta",
        "data/assembled/ecuador_intermediate_per_sample",
        "data/exports/H5N1_EC.fasta",
        "data/exports/H5N1_EC_summary.csv",
        "data/exports/H5N1_context.fasta",
        "data/exports/H5N1_context_summary.csv",
        "data/final/H5N1_final.fasta",
        "data/final/H5N1_final_beast.fasta",
        "data/final/H5N1_final_beast_summary.csv",
        "data/phylogeny/aligned/H5N1_full_concat_beast.mafft",
        expand("data/phylogeny/by_segment/H5N1_{segment}.fasta", segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]),
        "data/phylogeny/by_segment_summary.csv",
        "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE",
        "tmp_rf_full_concat/rf_summary.tsv",
        "data/beast_pre/panels/panel_main_taxa.tsv",
        "data/beast_pre/panels/panel_main.subset.mafft",
        "data/beast_pre/root_to_tip/treetime_clock.done",
        expand(
            "results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
            segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
        )
