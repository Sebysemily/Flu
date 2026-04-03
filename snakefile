configfile: "config/config.yml"

include: "rules/build_input_from_MIRA.smk"
include: "rules/01_ml_trees.smk"

include: "rules/00_validation_codon_rf.smk"

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
        expand(
            "results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
            segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
        )


rule rf_validation:
    input:
        # Codon+RF summary for 6 segments
        "results/validation/rf/rf_summary/rf_summary.tsv",
        expand(
            "results/validation/per_segment/{segment}/H5N1_{segment}_codon_validation.raxml.bestTreeCollapsed",
            segment=["PB2", "PB1", "PA", "HA", "NP", "NA"],
        ),
        # Simple RF summaries for NS and MP
        expand(
            "results/validation/rf/rf_summary_{segment}/rf_summary.tsv",
            segment=["NS", "MP"],
        )
