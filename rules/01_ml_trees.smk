PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
SIMPLE_SEGMENTS = ["NS", "MP"]
NP_MP_NS_SEGMENTS = ["NP", "MP", "NS"]
NP_MP_NS_CODON_SEGMENTS = ["NP"]
RANDOM_SEED = config.get("random_seed", 39809473)
MAX_THREADS = int(config.get("max_threads", 18))

DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/standard_header_input_fasta")
DATA_PHYLOGENY = config.get("data_phylogeny", "data/phylogeny")
RESULTS_PHYLOGENY = config.get("results_phylogeny", "results/phylogeny")
DATA_PRE_BEAST = config.get("data_pre_beast", "data/pre_beast")
DATA_BEAST = config.get("data_beast", "data/beast")
PROCESSED_ALIGNMENTS_CODON_AWARE = config.get("processed_alignments_codon_aware", "data/processed_alignments/codon_aware")
PROCESSED_ALIGNMENTS_QC_FILTERED = config.get("processed_alignments_qc_filtered", "data/processed_alignments/qc_filtered")
RESULTS_QC_METRICS = config.get("results_qc_metrics", "results/qc_metrics")
RESULTS_BEAST_PRE = config.get("results_beast_pre", "results/beast_pre")

QC_METRICS_DIR = f"{RESULTS_QC_METRICS}/processed_alignments/codon_aware"

FULL_CONCAT_FASTA = f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final_beast.fasta"
FULL_CONCAT_ALIGNMENT = f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat_beast.mafft"
FULL_CONCAT_PARTITIONS = f"{DATA_PHYLOGENY}/H5N1_full_concat_beast.partitions"
CONCAT_EXCEPT_HA_ALIGNMENT = f"{DATA_PHYLOGENY}/aligned/H5N1_concat_except_HA.mafft"
CONCAT_EXCEPT_HA_PARTITIONS = f"{DATA_PHYLOGENY}/H5N1_concat_except_HA.partitions"
FULL_CONCAT_PREFIX = f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat_beast"
NP_MP_NS_ALIGNMENT = f"{DATA_PHYLOGENY}/aligned/H5N1_NP_MP_NS.mafft"
NP_MP_NS_PARTITIONS = f"{DATA_PHYLOGENY}/H5N1_NP_MP_NS.partitions"
NP_MP_NS_PREFIX = f"{RESULTS_PHYLOGENY}/iq-tree/np_mp_ns/H5N1_NP_MP_NS"

CONCAT_TREE = f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat_beast.treefile"
CONCAT_ALIGNMENT = f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat_beast.mafft"

PRE_BEAST_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

DATA_BAYESIAN = "data/phylogeny/bayesian"
RESULTS_BAYESIAN = "results/phylogeny/bayesian"

FINAL_PANEL_DIR = "data/phylogeny/final_panel"
FINAL_ALIGNMENT = "data/phylogeny/final_alignment.fasta"
FINAL_PARTITIONS = "data/phylogeny/final_alignment.partitions"

RTT_DIR = "results/phylogeny/rtt"
RTT_DATA_DIR = "data/phylogeny/rtt"

BEAST_PANEL_TAXA = f"{RTT_DATA_DIR}/panel_main_taxa.tsv"
BEAST_PANEL_FILTERED_TAXA = f"{RTT_DATA_DIR}/panel_main_taxa.filtered.tsv"
BEAST_PANEL_RTT_FILTERED_TAXA = f"{DATA_BAYESIAN}/panel_main_taxa.final.tsv"

BEAST_FILTERED_SUBSET_ALIGNMENT = f"{RTT_DATA_DIR}/panel_main_concat.filtered.fasta"
BEAST_FILTERED_SUBSET_TREE = f"{RTT_DATA_DIR}/panel_main_concat.filtered.nwk"
BEAST_FILTERED_SUBSET_AUDIT = f"{RTT_DATA_DIR}/panel_main_concat.filtered.audit.tsv"

BAYESIAN_ALIGNMENT = f"{DATA_BAYESIAN}/bayesian_alignment.fasta"
BAYESIAN_ALIGNMENT_AUDIT = f"{DATA_BAYESIAN}/bayesian_alignment.audit.tsv"
BAYESIAN_ALIGNMENT_PARTITIONS = f"{DATA_BAYESIAN}/bayesian_alignment.partitions"

rule split_h5n1_final_by_segment:
    input:
        final_fasta=f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final.fasta"
    output:
        expand(f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS),
        f"{DATA_PHYLOGENY}/by_segment_summary.csv"
    shell:
        r"""
        python code/01_ml_trees/split_final_fasta_by_segment.py \
            --input-fasta {input.final_fasta} \
            --output-dir {DATA_PHYLOGENY}/by_segment \
            --summary-csv {DATA_PHYLOGENY}/by_segment_summary.csv
        """

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


rule build_single_segment_codon_partition:
    input:
        alignment=f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta"
    output:
        partition=f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.partitions"
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

rule concat_aligned_segments_with_partitions:
    input:
        alignments=expand(f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS)
    output:
        aligned=FULL_CONCAT_ALIGNMENT,
        partitions=FULL_CONCAT_PARTITIONS
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

rule concat_except_ha_with_partitions:
    input:
        alignments=expand(f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta", segment=["PB2", "PB1", "PA", "NP", "NA", "MP", "NS"])
    output:
        aligned=CONCAT_EXCEPT_HA_ALIGNMENT,
        partitions=CONCAT_EXCEPT_HA_PARTITIONS
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

def get_iqtree_replicate_seed(wildcards):
    return RANDOM_SEED + int(wildcards.rep)

rule run_iqtree_codon_segment_replicate:
    input:
        alignment=f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta",
        partitions=f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.partitions"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed
    wildcard_constraints:
        segment="|".join(CODON_SEGMENTS)
    threads: 4
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb 1000 -bnni
        mv {params.prefix}.treefile {output.treefile}
        """

rule run_iqtree_simple_segment_replicate:
    input:
        alignment=f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed
    wildcard_constraints:
        segment="|".join(SIMPLE_SEGMENTS)
    threads: 4
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb 1000 -bnni
        mv {params.prefix}.treefile {output.treefile}
        """

rule concat_except_HA_replicate:
    input:
        alignment=CONCAT_EXCEPT_HA_ALIGNMENT,
        partitions=CONCAT_EXCEPT_HA_PARTITIONS
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/concat_except_HA/rep{{rep}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/concat_except_HA/rep{{rep}}/rep{{rep}}",
        seed=get_iqtree_replicate_seed
    wildcard_constraints:
        segment="concat_except_HA"
    threads: 4
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb 1000 -bnni
        mv {params.prefix}.treefile {output.treefile}
        """

rule iqtree_fast_full_concat:
    input:
        alignment=FULL_CONCAT_ALIGNMENT,
        partitions=FULL_CONCAT_PARTITIONS
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat_beast.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/fast/fast",
        seed=RANDOM_SEED
    threads: MAX_THREADS
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        rm -f {params.prefix}.ckp.gz
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} --pathogen -seed {params.seed} -nt {threads} 
        mv {params.prefix}.treefile {output.treefile}
        """

rule prune_tree_by_distance_to_core:
    input:
        alignment=CONCAT_ALIGNMENT,
        tree=CONCAT_TREE,
    output:
        alignment=BEAST_FILTERED_SUBSET_ALIGNMENT,
        tree=BEAST_FILTERED_SUBSET_TREE,
        audit=BEAST_FILTERED_SUBSET_AUDIT,
        taxa=BEAST_PANEL_FILTERED_TAXA,
        discarded_csv=f"{RESULTS_QC_METRICS}/alignments/pruning_discarded_sequences.csv",
    conda:
        "../envs/01_ml_trees.yml"
    params:
        n_closest=config.get("prunning_parameters", {}).get("n_closest", 200),
        max_distance=config.get("prunning_parameters", {}).get("max_distance", 0.08),
        protect_anchors_per_month=config.get("prunning_parameters", {}).get("protect_anchors_per_month", 2),
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
            --protect-anchors-per-month {params.protect_anchors_per_month}
        """

# ---------------------------------------------------------------------------
# Alias rules — input-only rules with no shell/output.
# Purpose: provide clean named CLI targets (e.g. `snakemake align_all_segments`)
# and cleaner rulegraph nodes instead of embedding full rule names in DAG edges.
# ---------------------------------------------------------------------------

rule align_all_segments:
    input:
        expand(f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS)

rule iqtree_tree_full_concat:
    input:
        f"{FULL_CONCAT_PREFIX}.treefile"


rule generate_itol_colors_full_concat:
    input:
        tree=f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat_beast.treefile"
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


def get_itol_subset_tree(wildcards):
    if wildcards.subdir == "subset_fast":
        return f"{RESULTS_PHYLOGENY}/iq-tree/subset_fast/panel_main_concat.final.treefile"
    return f"{RESULTS_PHYLOGENY}/iq-tree/{wildcards.subdir}/{wildcards.subdir}_final.treefile"

rule itol_color_subsets:
    input:
        tree=get_itol_subset_tree
    output:
        tree_colors=f"{RESULTS_PHYLOGENY}/itol_styles/{{subdir}}/tree_colors.txt",
        colorstrip=f"{RESULTS_PHYLOGENY}/itol_styles/{{subdir}}/dataset_colorstrip.txt"
    wildcard_constraints:
        subdir="|".join(PRE_BEAST_SEGMENTS + ["subset_fast"])
    shell:
        r"""
        python code/01_ml_trees/generate_itol_colors.py \
            -i {input.tree} \
            --tree-colors {output.tree_colors} \
            --colorstrip {output.colorstrip}
        """

rule run_nextclade_alignment:
    input:
        fasta="data/gisaid_h5n1_raw.fasta",
        db="resources/nextclade_h5n1_ha"
    output:
        alignment="results/aligned_ha_cp1.fasta",
        tsv_report="results/nextclade_ha_report.tsv"
    conda:
        "../envs/01_ml_trees.yml"
    threads: 4
    shell:
        """
        nextclade run \
          --input-dataset {input.db} \
          --output-fasta {output.alignment} \
          --output-tsv {output.tsv_report} \
          --jobs {threads} \
          {input.fasta}
        """

rule generate_final_segment_tree:
    input:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep1.treefile"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}_final.treefile"
    shell:
        "cp {input.treefile} {output.treefile}"

rule calculate_concordance_factors:
    input:
        species_tree = f"{RESULTS_PHYLOGENY}/iq-tree/subset_fast/panel_main_concat.final.treefile",
        gene_trees = expand(f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}_final.treefile", segment=PHYLO_SEGMENTS),
        alignments = expand(f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS)
    output:
        gcf = "results/phylogeny/concordance/full_concat.gcs",
        scf = "results/phylogeny/concordance/full_concat.scf"
    params:
        aln_dir = FINAL_PANEL_DIR
    shell:
        r"""
        # Creamos un archivo que liste todos los árboles de los segmentos
        mkdir -p results/phylogeny/concordance
        ls {input.gene_trees} > results/phylogeny/concordance/gene_trees.list
        
        # Ejecutamos IQ-TREE 3
        iqtree -t {input.species_tree} \
               --gcf results/phylogeny/concordance/gene_trees.list \
               -p {params.aln_dir} \
               --scf 100 \
               -pre results/phylogeny/concordance/full_concat
        """

rule calculate_rf_distances:
    input:
        trees=expand(f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile", segment="{segment}", rep=range(1, 6))
    output:
        matrix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rf_dist_matrix.rf"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rf_dist_matrix",
        list_file=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/5_replicas.trees"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        # Concatenamos las 5 réplicas en un único archivo
        cat {input.trees} > {params.list_file}
        
        # Ejecutamos IQ-TREE para calcular las distancias Robinson-Foulds
        iqtree -rf_all {params.list_file} -pre {params.prefix}
        
        # Limpieza
        rm -f {params.list_file}
        """


rule apply_panel_to_segment_alignments:
    # Applies the distance-pruned taxa list to each per-segment MAFFT alignment,
    # and also builds the concatenated alignment and partitions for the subset.
    input:
        alignments=expand("data/processed_alignments/qc_filtered/H5N1_{segment}.mafft", segment=PRE_BEAST_SEGMENTS),
        taxa=BEAST_PANEL_FILTERED_TAXA,
    output:
        segment_alignments=expand(f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta", segment=PRE_BEAST_SEGMENTS),
        aligned=FINAL_ALIGNMENT,
        partitions=FINAL_PARTITIONS,
    params:
        segments=PRE_BEAST_SEGMENTS,
        segment_order=",".join(PRE_BEAST_SEGMENTS),
        codon_segments=",".join(CODON_SEGMENTS),
        final_panel_dir=FINAL_PANEL_DIR
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


rule iqtree_fast_subset:
    input:
        alignment=FINAL_ALIGNMENT,
        partitions=FINAL_PARTITIONS
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/subset_fast/panel_main_concat.final.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/subset_fast/fast/fast",
        seed=RANDOM_SEED
    threads: MAX_THREADS
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        rm -f {params.prefix}.ckp.gz
        iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -fast
        mv {params.prefix}.treefile {output.treefile}
        """

