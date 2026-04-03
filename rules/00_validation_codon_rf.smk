VALIDATION_CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
VALIDATION_SIMPLE_RF_SEGMENTS = ["NS", "MP"]
VALIDATION_ALL_8_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "NS", "MP"]
VALIDATION_DIR = "results/validation"
VALIDATION_CONCAT_DIR = f"{VALIDATION_DIR}/concat"
VALIDATION_RAXML_DIR = f"{VALIDATION_DIR}/raxml"
VALIDATION_RF_DIR = f"{VALIDATION_DIR}/rf"
VALIDATION_SEGMENT_DIR = f"{VALIDATION_DIR}/per_segment"

VALIDATION_CONCAT_ALIGNMENT = f"{VALIDATION_CONCAT_DIR}/H5N1_full_concat_codon_validation.mafft"
VALIDATION_CONCAT_PARTITIONS = f"{VALIDATION_CONCAT_DIR}/H5N1_full_concat_codon_validation.partitions"
VALIDATION_CONCAT_BASE_PREFIX = f"{VALIDATION_RAXML_DIR}/full_concat/H5N1_full_concat_codon_validation"
VALIDATION_CONCAT_BASE_TREE = f"{VALIDATION_CONCAT_BASE_PREFIX}.raxml.bestTreeCollapsed"

VALIDATION_RF_REPLICATES = [1, 2, 3, 4, 5]
VALIDATION_RF_REP_PREFIXES = {
    i: f"{VALIDATION_RF_DIR}/H5N1_full_concat_codon_validation_rep{i}" for i in VALIDATION_RF_REPLICATES
}
VALIDATION_RF_REP_TREES = {
    i: f"{VALIDATION_RF_REP_PREFIXES[i]}.raxml.bestTreeCollapsed" for i in VALIDATION_RF_REPLICATES
}
VALIDATION_RF_SUMMARY_TSV = f"{VALIDATION_RF_DIR}/rf_summary/rf_summary.tsv"

# Simple RF validation for NS and MP (no codon partitions)
VALIDATION_SIMPLE_RF_REP_PREFIXES = {
    segment: {
        i: f"{VALIDATION_RF_DIR}/H5N1_{segment}_simple_rf_rep{i}" for i in VALIDATION_RF_REPLICATES
    }
    for segment in VALIDATION_SIMPLE_RF_SEGMENTS
}
VALIDATION_SIMPLE_RF_REP_TREES = {
    segment: {
        i: f"{VALIDATION_SIMPLE_RF_REP_PREFIXES[segment][i]}.raxml.bestTreeCollapsed" for i in VALIDATION_RF_REPLICATES
    }
    for segment in VALIDATION_SIMPLE_RF_SEGMENTS
}
VALIDATION_SIMPLE_RF_SUMMARY_PATHS = {
    segment: f"{VALIDATION_RF_DIR}/rf_summary_{segment}/rf_summary.tsv" for segment in VALIDATION_SIMPLE_RF_SEGMENTS
}


rule concat_codon_validation_with_partitions:
    input:
        alignments=expand("data/phylogeny/aligned/H5N1_{segment}.mafft", segment=PHYLO_SEGMENTS),
        segment_trees=expand(
            "results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
            segment=PHYLO_SEGMENTS,
        )
    output:
        aligned=VALIDATION_CONCAT_ALIGNMENT,
        partitions=VALIDATION_CONCAT_PARTITIONS
    params:
        segment_order=",".join(PHYLO_SEGMENTS),
        codon_segments=",".join(VALIDATION_CODON_SEGMENTS)
    shell:
        r"""
        python code/00_validation_codon_rf/build_codon_validation_partitions.py \
            --segment-order {params.segment_order} \
            --codon-segments {params.codon_segments} \
            --output-alignment {output.aligned} \
            --output-partitions {output.partitions} \
            {input.alignments}
        """


rule raxml_ng_tree_full_concat_codon_validation_base:
    input:
        alignment=VALIDATION_CONCAT_ALIGNMENT,
        partitions=VALIDATION_CONCAT_PARTITIONS
    output:
        best_tree=VALIDATION_CONCAT_BASE_TREE
    params:
        prefix=VALIDATION_CONCAT_BASE_PREFIX,
        extra=lambda wildcards: f"--seed {RANDOM_SEED} --redo --tree 'pars{{20}},rand{{20}}'"
    threads: FULL_RAXML_THREADS
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        mkdir -p {VALIDATION_RAXML_DIR}/full_concat
        raxml-ng \
            --search \
            --msa {input.alignment} \
            --model {input.partitions} \
            --prefix {params.prefix} \
            --threads {threads} \
            {params.extra}
        """


rule raxml_ng_tree_full_concat_codon_validation_rf_rep:
    input:
        alignment=VALIDATION_CONCAT_ALIGNMENT,
        partitions=VALIDATION_CONCAT_PARTITIONS,
        base_tree=VALIDATION_CONCAT_BASE_TREE
    output:
        rep_tree=f"{VALIDATION_RF_DIR}/H5N1_full_concat_codon_validation_rep{{rep}}.raxml.bestTreeCollapsed"
    params:
        rep_prefix=lambda wildcards: VALIDATION_RF_REP_PREFIXES[int(wildcards.rep)],
        seed=lambda wildcards: RANDOM_SEED + 100 + int(wildcards.rep),
        extra=lambda wildcards: "--redo --tree 'pars{20},rand{20}'"
    threads: RF_RAXML_THREADS
    wildcard_constraints:
        rep="[1-5]"
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        mkdir -p {VALIDATION_RF_DIR}
        raxml-ng \
            --search \
            --msa {input.alignment} \
            --model {input.partitions} \
            --prefix {params.rep_prefix} \
            --threads {threads} \
            --seed {params.seed} \
            {params.extra}
        """


rule summarize_full_concat_codon_validation_rf_instability:
    input:
        base_tree=VALIDATION_CONCAT_BASE_TREE,
        replicate_trees=expand(
            f"{VALIDATION_RF_DIR}/H5N1_full_concat_codon_validation_rep{{rep}}.raxml.bestTreeCollapsed",
            rep=VALIDATION_RF_REPLICATES,
        )
    output:
        rf_summary=VALIDATION_RF_SUMMARY_TSV
    params:
        cutoff=0.05,
        base_seed=RANDOM_SEED,
        label="full_concat_codon_partition",
        base_label="full_concat_base",
        replicate_tree_args=lambda wildcards, input: " ".join(
            f"--replicate-tree {tree}" for tree in input.replicate_trees
        ),
        replicate_label_args=" ".join(
            f"--replicate-label full_concat_rep{rep}" for rep in VALIDATION_RF_REPLICATES
        ),
        replicate_seed_args=" ".join(
            f"--replicate-seed {RANDOM_SEED + 100 + rep}" for rep in VALIDATION_RF_REPLICATES
        )
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        python code/00_validation_codon_rf/rf_tree_instability_summary.py \
            --base-tree {input.base_tree} \
            {params.replicate_tree_args} \
            --base-seed {params.base_seed} \
            --base-label {params.base_label} \
            {params.replicate_label_args} \
            {params.replicate_seed_args} \
            --cutoff {params.cutoff} \
            --label {params.label} \
            --output {output.rf_summary}
        """


rule build_segment_codon_partitions_validation:
    input:
        alignment="data/phylogeny/aligned/H5N1_{segment}.mafft"
    output:
        partitions=f"{VALIDATION_SEGMENT_DIR}/{{segment}}/H5N1_{{segment}}.codon.partitions"
    wildcard_constraints:
        segment="PB2|PB1|PA|HA|NP|NA"
    shell:
        r"""
        python code/00_validation_codon_rf/build_single_segment_codon_partition.py \
            --alignment {input.alignment} \
            --segment {wildcards.segment} \
            --output {output.partitions}
        """


rule raxml_ng_tree_per_segment_codon_validation:
    input:
        alignment="data/phylogeny/aligned/H5N1_{segment}.mafft",
        partitions=f"{VALIDATION_SEGMENT_DIR}/{{segment}}/H5N1_{{segment}}.codon.partitions"
    output:
        best_tree=f"{VALIDATION_SEGMENT_DIR}/{{segment}}/H5N1_{{segment}}_codon_validation.raxml.bestTreeCollapsed"
    params:
        prefix=lambda wildcards: f"{VALIDATION_SEGMENT_DIR}/{wildcards.segment}/H5N1_{wildcards.segment}_codon_validation",
        extra=lambda wildcards: f"--seed {RANDOM_SEED} --redo --tree 'pars{{20}},rand{{20}}'"
    threads: SEGMENT_RAXML_THREADS
    wildcard_constraints:
        segment="PB2|PB1|PA|HA|NP|NA"
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        mkdir -p {VALIDATION_SEGMENT_DIR}/{wildcards.segment}
        raxml-ng \
            --search \
            --msa {input.alignment} \
            --model {input.partitions} \
            --prefix {params.prefix} \
            --threads {threads} \
            {params.extra}
        """


rule validation_segment_codon_trees:
    input:
        expand(
            f"{VALIDATION_SEGMENT_DIR}/{{segment}}/H5N1_{{segment}}_codon_validation.raxml.bestTreeCollapsed",
            segment=VALIDATION_CODON_SEGMENTS,
        )


rule validation_codon_rf_bundle:
    input:
        VALIDATION_RF_SUMMARY_TSV,
        expand(
            f"{VALIDATION_SEGMENT_DIR}/{{segment}}/H5N1_{{segment}}_codon_validation.raxml.bestTreeCollapsed",
            segment=VALIDATION_CODON_SEGMENTS,
        )


# ============================================================================
# Simple RF Validation for NS and MP (GTR+G, no codon partitions)
# ============================================================================

rule raxml_ng_tree_segment_simple_rf_base:
    """RF validation base tree for NS/MP segments (GTR+G model)."""
    input:
        alignment="data/phylogeny/aligned/H5N1_{segment}.mafft"
    output:
        best_tree=f"{VALIDATION_RF_DIR}/H5N1_{{segment}}_simple_rf_base.raxml.bestTreeCollapsed"
    params:
        prefix=lambda wildcards: f"{VALIDATION_RF_DIR}/H5N1_{wildcards.segment}_simple_rf_base",
        extra=lambda wildcards: f"--seed {RANDOM_SEED} --redo --tree 'pars{{20}},rand{{20}}'"
    threads: SEGMENT_RAXML_THREADS
    wildcard_constraints:
        segment="NS|MP"
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        mkdir -p {VALIDATION_RF_DIR}
        raxml-ng \
            --search \
            --msa {input.alignment} \
            --model GTR+G \
            --prefix {params.prefix} \
            --threads {threads} \
            {params.extra}
        """


rule raxml_ng_tree_segment_simple_rf_rep:
    """RF validation replicate trees for NS/MP segments (GTR+G model)."""
    input:
        alignment="data/phylogeny/aligned/H5N1_{segment}.mafft",
        base_tree=f"{VALIDATION_RF_DIR}/H5N1_{{segment}}_simple_rf_base.raxml.bestTreeCollapsed"
    output:
        rep_tree=f"{VALIDATION_RF_DIR}/H5N1_{{segment}}_simple_rf_rep{{rep}}.raxml.bestTreeCollapsed"
    params:
        rep_prefix=lambda wildcards: f"{VALIDATION_RF_DIR}/H5N1_{wildcards.segment}_simple_rf_rep{wildcards.rep}",
        seed=lambda wildcards: RANDOM_SEED + 100 + int(wildcards.rep),
        extra=lambda wildcards: "--redo --tree 'pars{20},rand{20}'"
    threads: SEGMENT_RAXML_THREADS
    wildcard_constraints:
        segment="NS|MP",
        rep="[1-5]"
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        raxml-ng \
            --search \
            --msa {input.alignment} \
            --model GTR+G \
            --prefix {params.rep_prefix} \
            --threads {threads} \
            --seed {params.seed} \
            {params.extra}
        """


rule summarize_segment_simple_rf_instability:
    """Compute RF distances for NS/MP per-segment validation."""
    input:
        base_tree=f"{VALIDATION_RF_DIR}/H5N1_{{segment}}_simple_rf_base.raxml.bestTreeCollapsed",
        replicate_trees=expand(
            f"{VALIDATION_RF_DIR}/H5N1_{{{{segment}}}}_simple_rf_rep{{rep}}.raxml.bestTreeCollapsed",
            rep=VALIDATION_RF_REPLICATES,
        )
    output:
        rf_summary=f"{VALIDATION_RF_DIR}/rf_summary_{{segment}}/rf_summary.tsv"
    params:
        cutoff=0.05,
        base_seed=RANDOM_SEED,
        label=lambda wildcards: f"{wildcards.segment}_simple_rf",
        base_label=lambda wildcards: f"{wildcards.segment}_base",
        replicate_tree_args=lambda wildcards, input: " ".join(
            f"--replicate-tree {tree}" for tree in input.replicate_trees
        ),
        replicate_label_args=lambda wildcards: " ".join(
            f"--replicate-label {wildcards.segment}_rep{rep}" for rep in VALIDATION_RF_REPLICATES
        ),
        replicate_seed_args=" ".join(
            f"--replicate-seed {RANDOM_SEED + 100 + rep}" for rep in VALIDATION_RF_REPLICATES
        )
    wildcard_constraints:
        segment="NS|MP"
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.rf_summary})
        python code/00_validation_codon_rf/rf_tree_instability_summary.py \
            --base-tree {input.base_tree} \
            {params.replicate_tree_args} \
            --base-seed {params.base_seed} \
            --base-label {params.base_label} \
            {params.replicate_label_args} \
            {params.replicate_seed_args} \
            --cutoff {params.cutoff} \
            --label {params.label} \
            --output {output.rf_summary}
        """


rule validation_rf_all_8_segments:
    """Bundle all RF validation results (codon for 6 segments + simple for NS,MP)."""
    input:
        # Codon+RF validation for 6 segments
        VALIDATION_RF_SUMMARY_TSV,
        # Simple RF validation for NS and MP
        expand(
            f"{VALIDATION_RF_DIR}/rf_summary_{{segment}}/rf_summary.tsv",
            segment=VALIDATION_SIMPLE_RF_SEGMENTS,
        )
