configfile: "config/config.yml"

include: "rules/build_input_from_MIRA.smk"
include: "rules/ml_per_segment.smk"

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
        expand("data/phylogeny/by_segment/H5N1_{segment}.fasta", segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]),
        "data/phylogeny/by_segment_summary.csv",
        expand(
            "results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportFBP",
            segment=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
        )
