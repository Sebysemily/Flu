RANDOM_SEED = config.get("random_seed", 39809473)
MAX_THREADS = int(config.get("max_threads", 18))

PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
SIMPLE_SEGMENTS = ["NS", "MP"]

# =====================================================================
# Rule: split_h5n1_final_by_segment
# =====================================================================

DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/standard_header_input_fasta")
DATA_PHYLOGENY = config.get("data_phylogeny", "data/phylogeny")

rule split_h5n1_final_by_segment:
    input:
        final_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final.fasta"
    output:
        expand(f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS)
    shell:
        r"""
        python code/01_ml_trees/split_final_fasta_by_segment.py \
            --input-fasta {input.final_fasta} \
            --output-dir {DATA_PHYLOGENY}/by_segment
        """

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
# Rule: concat_except_ha_with_partitions
# =====================================================================

rule concat_except_ha_with_partitions:
    input:
        alignments=expand(f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta", segment=["PB2", "PB1", "PA", "NP", "NA", "MP", "NS"])
    output:
        aligned=f"{DATA_PHYLOGENY}/aligned/H5N1_concat_except_HA.mafft",
        partitions=f"{DATA_PHYLOGENY}/aligned/H5N1_concat_except_HA.partitions"
    params:
        segment_order=",".join(["PB2", "PB1", "PA", "NP", "NA", "MP", "NS"]),
        codon_segments=",".join(["PB2", "PB1", "PA", "NP", "NA"])
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
# Rule: concat_except_HA_replicate
# =====================================================================

rule concat_except_HA_replicate:
    input:
        alignment=f"{DATA_PHYLOGENY}/aligned/H5N1_concat_except_HA.mafft",
        partitions=f"{DATA_PHYLOGENY}/aligned/H5N1_concat_except_HA.partitions"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/concat_except_HA/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/concat_except_HA/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed,
        bootstrap=config.get("iqtree_bootstrap", 1000)
    wildcard_constraints:
        segment="concat_except_HA"
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
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} --pathogen -seed {params.seed} -nt {threads} 
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

rule generate_itol_colors_full_concat:
    input:
        tree=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat.treefile"
    output:
        tree_colors=f"{RESULTS_PHYLOGENY}/itol_styles/tree_colors.txt",
        colorstrip=f"{RESULTS_PHYLOGENY}/itol_styles/dataset_colorstrip.txt"
    shell:
        r"""
        python code/01_ml_trees/generate_itol_colors.py \
            -i {input.tree} \
            --tree-colors {output.tree_colors} \
            --colorstrip {output.colorstrip}
        """

# =====================================================================
# Rule: itol_color_subsets
# =====================================================================

def get_itol_subset_tree(wildcards):
    if wildcards.subdir == "subset_concat":
        return f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/subset_concat_final.treefile"
    return f"{RESULTS_PHYLOGENY}/iq-tree/{wildcards.subdir}/{wildcards.subdir}_final.treefile"

rule itol_color_subsets:
    input:
        tree=get_itol_subset_tree
    output:
        tree_colors=f"{RESULTS_PHYLOGENY}/itol_styles/{{subdir}}/tree_colors.txt",
        colorstrip=f"{RESULTS_PHYLOGENY}/itol_styles/{{subdir}}/dataset_colorstrip.txt"
    wildcard_constraints:
        subdir="|".join(PHYLO_SEGMENTS + ["subset_concat"])
    shell:
        r"""
        python code/01_ml_trees/generate_itol_colors.py \
            -i {input.tree} \
            --tree-colors {output.tree_colors} \
            --colorstrip {output.colorstrip}
        """

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
# Rule: segment_analysis_au
# =====================================================================
rule segment_analysis_au:
    input:
        species_tree = f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/subset_concat_final.treefile",
        gene_trees = expand(f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}_final.treefile", segment=PHYLO_SEGMENTS),
        alignments = expand(f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS),
        concat_alignment = MAIN_PANEL_ALIGNMENT,
        concat_partitions = MAIN_PANEL_PARTITIONS,
    output:
        matrix = f"{RESULTS_PHYLOGENY}/segment_analysis/au_test_matrix.csv"
    params:
        segments = ",".join(PHYLO_SEGMENTS),
        alignments_dir = MAIN_PANEL_DIR,
        work_dir = f"{RESULTS_PHYLOGENY}/segment_analysis/au"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/segment_analysis/run_au_tests.py \
            --segments {params.segments} \
            --trees {input.gene_trees} \
            --concat-tree {input.species_tree} \
            --alignments-dir {params.alignments_dir} \
            --concat-alignment {input.concat_alignment} \
            --concat-partitions {input.concat_partitions} \
            --work-dir {params.work_dir} \
            --output {output.matrix}
        """

# =====================================================================
# Rule: segment_analysis_rf
# =====================================================================
rule segment_analysis_rf:
    input:
        species_tree = f"{RESULTS_PHYLOGENY}/iq-tree/subset_concat/subset_concat_final.treefile",
        gene_trees = expand(f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}_final.treefile", segment=PHYLO_SEGMENTS),
    output:
        matrix = f"{RESULTS_PHYLOGENY}/segment_analysis/rf_matrix.csv"
    params:
        segments = ",".join(PHYLO_SEGMENTS),
        work_dir = f"{RESULTS_PHYLOGENY}/segment_analysis/rf"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/segment_analysis/run_rf_matrix.py \
            --segments {params.segments} \
            --trees {input.gene_trees} \
            --concat-tree {input.species_tree} \
            --work-dir {params.work_dir} \
            --output {output.matrix}
        """

# =====================================================================
# Rule: segment_analysis_tanglegram
# =====================================================================
rule segment_analysis_tanglegram:
    input:
        ha_tree = f"{RESULTS_PHYLOGENY}/iq-tree/HA/HA_final.treefile",
        na_tree = f"{RESULTS_PHYLOGENY}/iq-tree/NA/NA_final.treefile",
        concat_except_ha_tree = f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram/concat_except_HA_final.treefile",
    output:
        ha_vs_na = f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram_ha_vs_na.png",
        ha_vs_except_ha = f"{RESULTS_PHYLOGENY}/segment_analysis/tanglegram_ha_vs_concat_except_ha.png",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        # Tanglegram 1: HA vs NA
        python code/segment_analysis/draw_tanglegrams.py \
            --tree1 {input.ha_tree} \
            --tree2 {input.na_tree} \
            --label1 "HA Segment" \
            --label2 "NA Segment" \
            --output {output.ha_vs_na}
            
        # Tanglegram 2: HA vs concat_except_HA
        python code/segment_analysis/draw_tanglegrams.py \
            --tree1 {input.ha_tree} \
            --tree2 {input.concat_except_ha_tree} \
            --label1 "HA Segment" \
            --label2 "All Segments except HA (Concat)" \
            --output {output.ha_vs_except_ha}
        """

# =====================================================================
# Rule: calculate_rf_distances
# =====================================================================

rule calculate_rf_distances:
    input:
        trees=lambda wildcards: expand(
            f"{RESULTS_PHYLOGENY}/iq-tree/{wildcards.segment}/rep{{rep}}.treefile",
            rep=range(1, int(config.get("iqtree_replicates", 5)) + 1)
        )
    output:
        matrix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rf_dist_matrix.rfdist"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rf_dist_matrix",
        list_file=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/5_replicas.trees"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        if [ $(echo {input.trees} | wc -w) -le 1 ]; then
            echo "Only one tree replica provided, writing dummy RF distance matrix."
            mkdir -p $(dirname {output.matrix})
            echo "Only one tree replica. RF distance calculation skipped." > {output.matrix}
        else
            cat {input.trees} > {params.list_file}
            iqtree -rf_all {params.list_file} -pre {params.prefix}
            rm -f {params.list_file}
        fi
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
