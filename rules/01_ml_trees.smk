RANDOM_SEED = config.get("random_seed", 39809473)
MAX_THREADS = int(config.get("max_threads", 18))

PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
SIMPLE_SEGMENTS = ["NS", "MP"]
MAIN_PANEL_METADATA = "metadata/H5N1_context.csv"

# =====================================================================
# =====================================================================

# =====================================================================
# Rule: download_nextclade_db
# =====================================================================

rule download_nextclade_db:
    output:
        ref=directory("resources/nextclade_h5n1_ha")
    conda:
        "../envs/01_ml_trees.yml"
    log:
        "logs/download_db.log"
    shell:
        """
        nextclade dataset get \
          --name 'community/moncla-lab/iav-h5/ha/2.3.4.4' \
          --output-dir {output.ref} > {log} 2>&1
        """

# =====================================================================
# Rule: extract_segment_reference
# =====================================================================

rule extract_segment_reference:
    input:
        fasta=f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta"
    output:
        ref=temp(f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}_ref.fasta")
    params:
        reference_id=config.get("reference_id", "Flu-0316")
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/extract_reference.py \
            --input-fasta {input.fasta} \
            --reference-id {params.reference_id} \
            --output-fasta {output.ref}
        """

# =====================================================================
# Rule: run_nextclade_alignment_per_segment
# =====================================================================

PROCESSED_ALIGNMENTS_CODON_AWARE = config.get("processed_alignments_codon_aware", "data/processed_alignments/codon_aware")
RESULTS_QC_METRICS = config.get("results_qc_metrics", "results/qc_metrics")

rule run_nextclade_alignment_per_segment:
    input:
        fasta=f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta",
        db="resources/nextclade_h5n1_ha",
        ref_seq=f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}_ref.fasta"
    output:
        trimmed=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft",
        report=f"{RESULTS_QC_METRICS}/alignments/nextclade_{{segment}}_report.tsv"
    conda:
        "../envs/01_ml_trees.yml"
    threads: MAX_THREADS
    shell:
        r"""
        if [ "{wildcards.segment}" = "HA" ]; then
            nextclade run \
                --input-dataset {input.db} \
                --output-fasta {output.trimmed} \
                --output-tsv {output.report} \
                --jobs {threads} \
                {input.fasta}
        else
            nextclade run \
                --input-ref {input.ref_seq} \
                --output-fasta {output.trimmed} \
                --output-tsv {output.report} \
                --jobs {threads} \
                {input.fasta}
        fi
        """

# =====================================================================
# Rule: filter_alignments_by_nextclade_qc
# =====================================================================

PROCESSED_ALIGNMENTS_QC_FILTERED = config.get("processed_alignments_qc_filtered", "data/processed_alignments/qc_filtered")

rule filter_alignments_by_nextclade_qc:
    input:
        alignments=expand(f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS),
        report=f"{RESULTS_QC_METRICS}/alignments/nextclade_HA_report.tsv"
    output:
        filtered=expand(f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS),
        discarded_csv=f"{RESULTS_QC_METRICS}/alignments/qc_discarded_sequences.csv"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/filter_qc_nextclade.py \
            --input-dir {PROCESSED_ALIGNMENTS_CODON_AWARE} \
            --report {input.report} \
            --output-dir {PROCESSED_ALIGNMENTS_QC_FILTERED} \
            --discarded-csv {output.discarded_csv}
        """

# =====================================================================
# Rule: build_single_segment_codon_partition
# =====================================================================

MAIN_PANEL_DIR = "data/phylogeny/main_panel"
MAIN_PANEL_ALIGNMENT = f"{MAIN_PANEL_DIR}/main_alignment.fasta"
MAIN_PANEL_PARTITIONS = f"{MAIN_PANEL_DIR}/main_alignment.partitions"

rule build_single_segment_codon_partition:
    input:
        alignment=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta"
    output:
        partition=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.partitions"
    wildcard_constraints:
        segment="|".join(CODON_SEGMENTS)
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/build_single_segment_codon_partition.py \
            --alignment {input.alignment} \
            --segment {wildcards.segment} \
            --output {output.partition}
        """

# =====================================================================
# Rule: concat_aligned_segments_with_partitions
# =====================================================================

rule concat_aligned_segments_with_partitions:
    input:
        alignments=expand(f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS)
    output:
        aligned=f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat.mafft",
        partitions=f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat.partitions"
    params:
        segment_order=",".join(PHYLO_SEGMENTS),
        codon_segments=",".join(CODON_SEGMENTS)
    shell:
        r"""
        python code/01_ml_trees/build_concat_codon_partitions.py \
            --segment-order {params.segment_order} \
            --codon-segments {params.codon_segments} \
            --output-alignment {output.aligned} \
            --output-partitions {output.partitions} \
            {input.alignments}
        """

# =====================================================================
# Rule: run_iqtree_codon_segment_replicate
# =====================================================================

def get_iqtree_replicate_seed(wildcards):
    return RANDOM_SEED + int(wildcards.rep)

RESULTS_PHYLOGENY = config.get("results_phylogeny", "results/phylogeny")

rule run_iqtree_codon_segment_replicate:
    input:
        alignment=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta",
        partitions=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.partitions"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed,
        bootstrap=config.get("iqtree_bootstrap", 1000)
    wildcard_constraints:
        segment="|".join(CODON_SEGMENTS)
    threads: 4
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb {params.bootstrap} -bnni
        cp {params.prefix}.treefile {output.treefile}
        """

# =====================================================================
# Rule: run_iqtree_simple_segment_replicate
# =====================================================================

rule run_iqtree_simple_segment_replicate:
    input:
        alignment=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed,
        bootstrap=config.get("iqtree_bootstrap", 1000)
    wildcard_constraints:
        segment="|".join(SIMPLE_SEGMENTS)
    threads: 4
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb {params.bootstrap} -bnni
        cp {params.prefix}.treefile {output.treefile}
        """

# =====================================================================
# Rule: iqtree_fast_full_concat
# =====================================================================

rule iqtree_fast_full_concat:
    input:
        alignment=f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat.mafft",
        partitions=f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat.partitions"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/pathogen/pathogen",
        seed=RANDOM_SEED
    threads: MAX_THREADS
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} --fast -seed {params.seed} -nt {threads} 
        cp {params.prefix}.treefile {output.treefile}
        """

# =====================================================================
# Rule: prune_tree_by_distance_to_core
# =====================================================================

RTT_DIR = "results/phylogeny/rtt"
RTT_DATA_DIR = "data/phylogeny/rtt"
MAIN_PANEL_FILTERED_TAXA = f"{RTT_DATA_DIR}/panel_main_taxa.filtered.tsv"
MAIN_PANEL_FILTERED_ALIGNMENT = f"{RTT_DATA_DIR}/panel_main_concat.filtered.fasta"
MAIN_PANEL_FILTERED_TREE = f"{RTT_DATA_DIR}/panel_main_concat.filtered.nwk"
MAIN_PANEL_FILTERED_AUDIT = f"{RTT_DATA_DIR}/panel_main_concat.filtered.audit.tsv"

rule prune_tree_by_distance_to_core:
    input:
        alignment=f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat.mafft",
        tree=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat.treefile",
    output:
        alignment=MAIN_PANEL_FILTERED_ALIGNMENT,
        tree=MAIN_PANEL_FILTERED_TREE,
        audit=MAIN_PANEL_FILTERED_AUDIT,
        taxa=MAIN_PANEL_FILTERED_TAXA,
        discarded_csv=f"{RESULTS_QC_METRICS}/alignments/pruning_discarded_sequences.csv",
    conda:
        "../envs/01_ml_trees.yml"
    params:
        n_closest=config.get("prunning_parameters", {}).get("n_closest", 200),
        max_distance=config.get("prunning_parameters", {}).get("max_distance", 0.08),
        protect_anchors_per_month=config.get("prunning_parameters", {}).get("protect_anchors_per_month", 2),
        protect_regional_per_month=config.get("prunning_parameters", {}).get("protect_regional_per_month", 3),
        temp_audit=f"{RTT_DATA_DIR}/panel_distance_audit.tsv",
    shell:
        r"""
        python code/01_ml_trees/prune_panel_by_distance.py \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --panel-main-out {output.taxa} \
            --out-alignment {output.alignment} \
            --out-tree {output.tree} \
            --audit-out {output.audit} \
            --distance-audit-out {params.temp_audit} \
            --discarded-csv {output.discarded_csv} \
            --n-closest {params.n_closest} \
            --max-distance {params.max_distance} \
            --protect-anchors-per-month {params.protect_anchors_per_month} \
            --protect-regional-per-month {params.protect_regional_per_month}
        """

# =====================================================================
# Rule: generate_itol_colors_full_concat
# =====================================================================

# =====================================================================
# Rule: copy_trees_for_segment_analysis
# =====================================================================
rule copy_trees_for_segment_analysis:
    input:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep1.treefile"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/{{dir_path}}/{{segment}}_final.treefile"
    wildcard_constraints:
        dir_path=r"(iq-tree/[^/]+|segment_analysis/tanglegram)"
    shell:
        "cp {input.treefile} {output.treefile}"

# =====================================================================
# Rule: segment_analysis_tanglegram
# =====================================================================
rule segment_analysis_tanglegram:
    input:
        ha_tree = f"{RESULTS_PHYLOGENY}/iq-tree/HA/HA_final.treefile",
        pb2_tree = f"{RESULTS_PHYLOGENY}/iq-tree/PB2/PB2_final.treefile",
        metadata = MAIN_PANEL_METADATA,
        ecuador_metadata = "config/flu_filtrado.csv"
    output:
        tanglegram = f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram_ha_vs_pb2.png",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/segment_analysis/draw_tanglegrams.py \
            --tree1 {input.ha_tree} \
            --tree2 {input.pb2_tree} \
            --label1 "HA Segment" \
            --label2 "PB2 Segment" \
            --metadata {input.metadata} \
            --ecuador-metadata {input.ecuador_metadata} \
            --collapse-context \
            --hide-non-core-labels \
            --output {output.tanglegram}
        """

# =====================================================================
# Rule: apply_panel_to_segment_alignments
# =====================================================================

rule apply_panel_to_segment_alignments:
    # Applies the distance-pruned taxa list to each per-segment MAFFT alignment,
    # and also builds the concatenated alignment and partitions for the subset.
    input:
        alignments=expand("data/processed_alignments/qc_filtered/H5N1_{segment}.mafft", segment=PHYLO_SEGMENTS),
        taxa=MAIN_PANEL_FILTERED_TAXA,
    output:
        segment_alignments=expand(f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS),
        aligned=MAIN_PANEL_ALIGNMENT,
        partitions=MAIN_PANEL_PARTITIONS,
    params:
        segments=PHYLO_SEGMENTS,
        segment_order=",".join(PHYLO_SEGMENTS),
        codon_segments=",".join(CODON_SEGMENTS),
        final_panel_dir=MAIN_PANEL_DIR
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        # 1. Subset each segment alignment
        for seg in {params.segments}; do
            python code/02_Beast/subset_alignment_by_taxa.py \
                --alignment data/processed_alignments/qc_filtered/H5N1_${{seg}}.mafft \
                --taxa {input.taxa} \
                --out-alignment {params.final_panel_dir}/H5N1_${{seg}}.fasta
        done
        
        # 2. Build concatenated alignment and partitions
        python code/01_ml_trees/build_concat_codon_partitions.py \
            --segment-order {params.segment_order} \
            --codon-segments {params.codon_segments} \
            --output-alignment {output.aligned} \
            --output-partitions {output.partitions} \
            {output.segment_alignments}
        """

# =====================================================================
# =====================================================================
# Rule: run_iqtree_subset_concat_replicate
# =====================================================================

rule run_iqtree_subset_concat_replicate:
    input:
        alignment=MAIN_PANEL_ALIGNMENT,
        partitions=MAIN_PANEL_PARTITIONS
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed,
        bootstrap=config.get("iqtree_bootstrap", 1000)
    threads: MAX_THREADS
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb {params.bootstrap} -bnni
        cp {params.prefix}.treefile {output.treefile}
        """

# =====================================================================
# Rule: segment_analysis_ggtree
# =====================================================================
rule segment_analysis_ggtree:
    input:
        tree=lambda wildcards: (
            f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/subset_concat_final.treefile"
            if wildcards.segment == "subset_concat"
            else f"{RESULTS_PHYLOGENY}/iq-tree/{wildcards.segment}/{wildcards.segment}_final.treefile"
        ),
        metadata=MAIN_PANEL_METADATA,
        ecuador_metadata="config/flu_filtrado.csv"
    output:
        png="results/figures/{segment}_tree.png"
    conda:
        "../envs/ggtree.yml"
    shell:
        r"""
        Rscript code/segment_analysis/plot_ggtree.R \
            {input.tree} \
            {input.metadata} \
            {input.ecuador_metadata} \
            {output.png} \
            "{wildcards.segment} Tree"
        """
