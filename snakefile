configfile: "config/config.yml"

include: "rules/build_gisaid_input_from_mira.smk"
include: "rules/build_inputs.smk"
include: "rules/01_ml_trees.smk"

include: "rules/00_validation_codon_rf.smk"

rule all:
    input:
        "data/input/H5N1_context.fasta",
        "data/input/H5N1_context_summary.csv",
        "data/final/H5N1_final.fasta",
        "data/final/H5N1_final_beast.fasta",
        "data/final/H5N1_final_beast_summary.csv",
        "data/phylogeny/aligned/H5N1_full_concat_beast.mafft",
        expand("data/phylogeny/by_segment/H5N1_{segment}.fasta", segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]),
        "data/phylogeny/by_segment_summary.csv",
        "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE",
        expand(
            "results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
            segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
        )


rule build_gisaid:
    input:
        "data/all_amended_fasta",
        "data/assembled/ecuador_intermediate_summary.csv",
        "data/assembled/ecuador_intermediate_audit.csv",
        "data/assembled/ecuador_intermediate_issues.csv",
        "data/assembled/ecuador_intermediate_sequences.fasta",
        "data/assembled/ecuador_intermediate_per_sample",
        "data/input/H5N1_EC.fasta",
        "data/input/H5N1_EC_summary.csv"


rule rf_validation:
    input:
        # Concat codon+RF summary
        "results/validation/rf/rf_summary/rf_summary.tsv",
        # Codon RF summaries for 6 segments
        expand(
            "results/validation/rf/rf_summary_{segment}_codon/rf_summary.tsv",
            segment=["PB2", "PB1", "PA", "HA", "NP", "NA"],
        ),
        # Simple RF summaries for NS and MP
        expand(
            "results/validation/rf/rf_summary_{segment}/rf_summary.tsv",
            segment=["NS", "MP"],
        )
