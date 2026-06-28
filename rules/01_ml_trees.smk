RANDOM_SEED = config.get("random_seed", 39809473)
MAX_THREADS = int(config.get("max_threads", 18))

PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
LINEAGE_SEGMENTS = ["HA", "PB2"]
CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
SIMPLE_SEGMENTS = ["NS", "MP"]
MAIN_PANEL_METADATA = "metadata/H5N1_context.csv"
FILTRADO_CSV = config.get("flu_filtrado", "metadata/flu_filtrado.csv")

# =====================================================================
# =====================================================================

# =====================================================================
# Rule: download_nextclade_db
# =====================================================================

def get_nextclade_dataset(segment):
    if segment == "HA":
        return "community/moncla-lab/iav-h5/ha/2.3.4.4"
    elif segment == "NA":
        return "nextstrain/flu/h1n1pdm/na/MW626056"
    else:
        return f"nextstrain/flu/h1n1pdm/{segment.lower()}"

rule download_nextclade_db:
    output:
        ref=directory("resources/nextclade_{segment}")
    params:
        dataset=lambda wildcards: get_nextclade_dataset(wildcards.segment)
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        """
        mkdir -p {output.ref}
        nextclade dataset get \
          --name '{params.dataset}' \
          --output-dir {output.ref}
        """

# =====================================================================
# Rule: run_nextclade_alignment_all
# =====================================================================

PROCESSED_ALIGNMENTS_CODON_AWARE = config.get("processed_alignments_codon_aware", "data/processed_alignments/codon_aware")
RESULTS_QC_METRICS = config.get("results_qc_metrics", "results/qc_metrics")
# Paso 1: trim gap/N tras MAFFT global (todos los segmentos)
TRIM_GAP_N_QC_DIR = f"{RESULTS_QC_METRICS}/trim_gap_n"
# Paso 2: trim gap/N solo panel HA (BEAST)
HA_ONLY_QC_DIR = f"{RESULTS_QC_METRICS}/HA_only"
TRIM_MAX_DIVERGENCE = 0.10

rule run_nextclade_alignment_all:
    input:
        fasta=f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta",
        db="resources/nextclade_{segment}"
    output:
        trimmed=f"{RESULTS_QC_METRICS}/nextclade/{{segment}}/alignment.fasta",
        report=f"{RESULTS_QC_METRICS}/nextclade/{{segment}}/nextclade_report.csv"
    conda:
        "../envs/01_ml_trees.yml"
    threads: MAX_THREADS
    shell:
        r"""
        nextclade run \
            --input-dataset {input.db} \
            --output-fasta {output.trimmed} \
            --output-csv {output.report} \
            --jobs {threads} \
            {input.fasta}
        """


# =====================================================================
# Rule: filter_alignments_by_nextclade_qc
# =====================================================================

PROCESSED_ALIGNMENTS_QC_FILTERED = config.get("processed_alignments_qc_filtered", "data/processed_alignments/qc_filtered")

rule filter_alignments_by_nextclade_qc_ha:
    input:
        alignment=f"{DATA_PHYLOGENY}/by_segment/H5N1_HA.fasta",
        report=f"{RESULTS_QC_METRICS}/nextclade/HA/nextclade_report.csv",
        metadata=FILTRADO_CSV,
    output:
        filtered=f"{DATA_PHYLOGENY}/by_segment_qc_filtered/H5N1_HA.fasta",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/filter_qc_nextclade.py \
            --input-alignment {input.alignment} \
            --report {input.report} \
            --output-alignment {output.filtered} \
            --metadata {input.metadata}
        """


rule report_nextclade_ha_qc_discarded:
    """Write eliminated-sample CSVs only (does not touch filtered FASTA)."""
    input:
        report=f"{RESULTS_QC_METRICS}/nextclade/HA/nextclade_report.csv",
        metadata=FILTRADO_CSV,
        role_metadata=MAIN_PANEL_METADATA,
    output:
        discarded_csv=f"{RESULTS_QC_METRICS}/nextclade/HA_eliminated.csv",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/filter_qc_nextclade.py \
            --report {input.report} \
            --metadata {input.metadata} \
            --role-metadata {input.role_metadata} \
            --discarded-csv {output.discarded_csv} \
            --filter-step nextclade_ha \
            --discarded-csv-only
        """

rule filter_alignments_by_nextclade_qc_others:
    input:
        alignment=f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta"
    output:
        filtered=f"{DATA_PHYLOGENY}/by_segment_qc_filtered/H5N1_{{segment}}.fasta"
    wildcard_constraints:
        segment="PB2|PB1|PA|NP|NA|MP|NS"
    shell:
        "cp {input.alignment} {output.filtered}"

rule mafft_align_segment:
    input:
        fasta=f"{DATA_PHYLOGENY}/by_segment_qc_filtered/H5N1_{{segment}}.fasta"
    output:
        alignment=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft"
    conda:
        "../envs/01_ml_trees.yml"
    threads: MAX_THREADS
    shell:
        r"""
        mafft --auto --thread {threads} {input.fasta} > {output.alignment}
        """

rule trim_and_filter_mafft:
    input:
        alignment=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft",
    output:
        filtered=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/trim_and_filter_mafft.py \
            --input {input.alignment} \
            --output {output.filtered} \
            --max-divergence {TRIM_MAX_DIVERGENCE}
        """


rule report_trim_gap_n_discarded:
    """Paso 1: CSV de eliminados por trim MAFFT (todos los segmentos); no toca alignments filtrados."""
    input:
        alignment=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft",
        role_metadata=MAIN_PANEL_METADATA,
    output:
        discarded_csv=f"{TRIM_GAP_N_QC_DIR}/{{segment}}_eliminated.csv",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/trim_and_filter_mafft.py \
            --input {input.alignment} \
            --max-divergence {TRIM_MAX_DIVERGENCE} \
            --role-metadata {input.role_metadata} \
            --discarded-csv {output.discarded_csv} \
            --filter-step trim_gap_n_{wildcards.segment} \
            --discarded-csv-only
        """






# =====================================================================
# Rule: build_single_segment_codon_partition
# =====================================================================

MAIN_PANEL_DIR = "data/beast/main_panel"
PRUNED_HA_DIR = f"{MAIN_PANEL_DIR}/pre-alignment"
# HA panel: unaligned subset -> MAFFT -> trim/QC (canonical input for IQ-TREE, RTT, BEAST metadata)
MAIN_PANEL_HA_UNALIGNED = f"{PRUNED_HA_DIR}/H5N1_HA_panel_unaligned.fasta"
MAIN_PANEL_HA_MAFFT = f"{PRUNED_HA_DIR}/H5N1_HA_panel_mafft.mafft"
MAIN_PANEL_HA_POSTQC = f"{MAIN_PANEL_DIR}/H5N1_HA_panel_postQC.fasta"
MAIN_PANEL_CONTEXT_TAXA = f"{PRUNED_HA_DIR}/context_taxa.txt"
MAIN_PANEL_COSTA_AUDIT = f"{PRUNED_HA_DIR}/costa_window.audit.tsv"
MAIN_PANEL_FILTERED_AUDIT = f"{PRUNED_HA_DIR}/pruning.audit.tsv"
# HA panel augur + post-QC (hardcoded; edit here if needed)
PANEL_REGIONAL_CONTEXT_MAX = 250
PANEL_AMERICAN_ANCHOR_MAX = 100
PANEL_REGIONAL_COSTA_ADJACENT_MAX = 0
PANEL_COSTA_WINDOW_PADDING_MONTHS = 2
MAIN_PANEL_ALIGNMENT = f"{MAIN_PANEL_DIR}/main_alignment.fasta"
MAIN_PANEL_PARTITIONS = f"{MAIN_PANEL_DIR}/main_alignment.partitions"


def main_panel_segment_fasta(segment):
    if segment == "HA":
        return MAIN_PANEL_HA_POSTQC
    return f"{MAIN_PANEL_DIR}/H5N1_{segment}.fasta"


rule build_single_segment_codon_partition:
    input:
        alignment=lambda wildcards: main_panel_segment_fasta(wildcards.segment),
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
# Rule: run_iqtree_codon_segment
# =====================================================================

RESULTS_PHYLOGENY = config.get("results_phylogeny", "results/phylogeny")

rule run_iqtree_codon_segment:
    input:
        alignment=lambda wildcards: main_panel_segment_fasta(wildcards.segment),
        partitions=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.partitions"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/workdir/{{segment}}",
        seed=RANDOM_SEED,
        bootstrap=1000
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
# Rule: run_iqtree_simple_segment
# =====================================================================

rule run_iqtree_simple_segment:
    input:
        alignment=lambda wildcards: main_panel_segment_fasta(wildcards.segment),
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/{{segment}}.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/workdir/{{segment}}",
        seed=RANDOM_SEED,
        bootstrap=1000
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
# Rule: build_lineage_codon_partition
# =====================================================================

rule build_lineage_codon_partition:
    input:
        alignment=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft"
    output:
        partition=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.partitions"
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
# Rule: iqtree_fast_codon_segment
# =====================================================================

rule iqtree_fast_codon_segment:
    input:
        alignment=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft",
        partitions=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.partitions"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/H5N1_{{segment}}_fast.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/fast_run/fast_{{segment}}",
        seed=RANDOM_SEED,
        bootstrap=1000
    wildcard_constraints:
        segment="|".join(CODON_SEGMENTS)
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
# Rule: iqtree_fast_simple_segment
# =====================================================================

rule iqtree_fast_simple_segment:
    input:
        alignment=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/H5N1_{{segment}}_fast.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/fast_run/fast_{{segment}}",
        seed=RANDOM_SEED,
        bootstrap=1000
    wildcard_constraints:
        segment="|".join(SIMPLE_SEGMENTS)
    threads: MAX_THREADS
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree -s {input.alignment} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb {params.bootstrap} -bnni
        cp {params.prefix}.treefile {output.treefile}
        """

# =====================================================================
# Rule: segment_analysis_ggtree_lineage (HA / PB2 fast trees)
# =====================================================================

LINEAGE_GGTREE_SEGMENTS = ["HA", "PB2"]
LINEAGE_FIGURES_DIR = "figures/HA_PB2_lineage"

rule segment_analysis_ggtree_lineage:
    input:
        tree=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/H5N1_{{segment}}_fast.treefile",
        metadata=MAIN_PANEL_METADATA,
    output:
        png=f"{LINEAGE_FIGURES_DIR}/{{segment}}_lineage_fast_tree.png",
        rds=f"{LINEAGE_FIGURES_DIR}/{{segment}}_lineage_fast_tree.rds",
    params:
        title=lambda wildcards: f"{wildcards.segment} lineage (PB2+HA panel, fast ML)",
        figures_dir=LINEAGE_FIGURES_DIR,
        outgroup_sample=config.get("outgroup_root_sample", "EPI_ISL_18133416")
    wildcard_constraints:
        segment="|".join(LINEAGE_GGTREE_SEGMENTS),
    conda:
        "../envs/r.yml"
    shell:
        r"""
        mkdir -p {params.figures_dir}
        Rscript code/segment_analysis/plot_ggtree.R \
            {input.tree} \
            {input.metadata} \
            {output.png} \
            "{params.title}" \
            "{params.outgroup_sample}"
        test -s {output.png}
        test -s {output.rds}
        """

# =====================================================================
# Rule: segment_analysis_composite_lineage (heatmap + HA/PB2 tanglegram)
# =====================================================================

rule segment_analysis_composite_lineage_preview:
    input:
        ha_tree=f"{RESULTS_PHYLOGENY}/iq-tree/HA/lineage/H5N1_HA_fast.treefile",
        pb2_tree=f"{RESULTS_PHYLOGENY}/iq-tree/PB2/lineage/H5N1_PB2_fast.treefile",
        panel_metadata=MAIN_PANEL_METADATA,
    output:
        png=f"{LINEAGE_FIGURES_DIR}/HA_PB2_lineage_composite_preview.png",
    params:
        figures_dir=LINEAGE_FIGURES_DIR,
        ribbon_segment="HA",
        outgroup_sample=config.get("outgroup_root_sample", ""),
    conda:
        "../envs/r.yml"
    shell:
        r"""
        ulimit -s unlimited || true
        mkdir -p {params.figures_dir}
        Rscript code/segment_analysis/plot_composite_lineage.R \
            {input.ha_tree} \
            {input.pb2_tree} \
            {input.panel_metadata} \
            {output.png} \
            60 \
            {params.ribbon_segment} \
            "{params.outgroup_sample}"
        test -s {output.png}
        """

rule segment_analysis_composite_lineage:
    input:
        ha_tree=f"{RESULTS_PHYLOGENY}/iq-tree/HA/lineage/H5N1_HA_fast.treefile",
        pb2_tree=f"{RESULTS_PHYLOGENY}/iq-tree/PB2/lineage/H5N1_PB2_fast.treefile",
        panel_metadata=MAIN_PANEL_METADATA,
    output:
        png=f"{LINEAGE_FIGURES_DIR}/HA_PB2_lineage_composite.png",
    params:
        figures_dir=LINEAGE_FIGURES_DIR,
        ribbon_segment="HA",
        outgroup_sample=config.get("outgroup_root_sample", ""),
    conda:
        "../envs/r.yml"
    shell:
        r"""
        ulimit -s unlimited || true
        mkdir -p {params.figures_dir}
        Rscript code/segment_analysis/plot_composite_lineage.R \
            {input.ha_tree} \
            {input.pb2_tree} \
            {input.panel_metadata} \
            {output.png} \
            "" \
            {params.ribbon_segment} \
            "{params.outgroup_sample}"
        test -s {output.png}
        """

# =====================================================================
rule augur_filter_context:
    input:
        metadata=MAIN_PANEL_METADATA,
        ha_fasta=f"{DATA_PHYLOGENY}/by_segment/H5N1_HA.fasta",
    output:
        context_taxa=MAIN_PANEL_CONTEXT_TAXA,
        include_list=f"{PRUNED_HA_DIR}/augur_include.txt",
        costa_audit=MAIN_PANEL_COSTA_AUDIT,
    conda:
        "../envs/01_ml_trees.yml"
    params:
        outgroup_root_sample=config.get("outgroup_root_sample", "EPI_ISL_18133416"),
        seed=RANDOM_SEED,
        regional_context_max=PANEL_REGIONAL_CONTEXT_MAX,
        american_anchor_max=PANEL_AMERICAN_ANCHOR_MAX,
        regional_costa_adjacent_max=PANEL_REGIONAL_COSTA_ADJACENT_MAX,
        costa_window_padding_months=PANEL_COSTA_WINDOW_PADDING_MONTHS,
    shell:
        r"""
        python code/01_ml_trees/build_panel_context_taxa.py \
            --metadata {input.metadata} \
            --ha-fasta {input.ha_fasta} \
            --out-context-taxa {output.context_taxa} \
            --out-include-list {output.include_list} \
            --out-costa-audit {output.costa_audit} \
            --outgroup {params.outgroup_root_sample} \
            --seed {params.seed} \
            --regional-context-max {params.regional_context_max} \
            --american-anchor-max {params.american_anchor_max} \
            --regional-costa-adjacent-max {params.regional_costa_adjacent_max} \
            --costa-window-padding-months {params.costa_window_padding_months}
        """


# =====================================================================
# Rule: subset_segments_pre_alignment
# =====================================================================

# =====================================================================
# Rule: subset_ha_pre_alignment
# =====================================================================

rule subset_ha_pre_alignment:
    input:
        unaligned_fasta=f"{DATA_PHYLOGENY}/by_segment/H5N1_HA.fasta",
        context_list=MAIN_PANEL_CONTEXT_TAXA,
        metadata=MAIN_PANEL_METADATA,
    output:
        fasta=MAIN_PANEL_HA_UNALIGNED,
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/subset_segment_fastas.py \
            --unaligned-fasta {input.unaligned_fasta} \
            --context-list {input.context_list} \
            --metadata {input.metadata} \
            --out-fasta {output.fasta}
        """

# =====================================================================
# Rule: run_nextclade_alignment_pruned_ha
# =====================================================================

rule run_mafft_alignment_pruned_ha:
    input:
        fasta=MAIN_PANEL_HA_UNALIGNED
    output:
        trimmed=MAIN_PANEL_HA_MAFFT
    conda:
        "../envs/01_ml_trees.yml"
    threads: MAX_THREADS
    shell:
        r"""
        mafft --auto --thread {threads} {input.fasta} > {output.trimmed}
        """

rule trim_and_filter_pruned_ha:
    input:
        alignment=MAIN_PANEL_HA_MAFFT,
    output:
        filtered=MAIN_PANEL_HA_POSTQC,
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/trim_and_filter_mafft.py \
            --input {input.alignment} \
            --output {output.filtered} \
            --max-divergence {TRIM_MAX_DIVERGENCE}
        """


rule report_panel_ha_trim_discarded:
    """Paso 2: CSV de eliminados trim panel HA (HA_only); no toca postQC."""
    input:
        alignment=MAIN_PANEL_HA_MAFFT,
        role_metadata=MAIN_PANEL_METADATA,
    output:
        discarded_csv=f"{HA_ONLY_QC_DIR}/panel_HA_eliminated.csv",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/trim_and_filter_mafft.py \
            --input {input.alignment} \
            --max-divergence {TRIM_MAX_DIVERGENCE} \
            --role-metadata {input.role_metadata} \
            --discarded-csv {output.discarded_csv} \
            --filter-step trim_gap_n_panel_HA \
            --discarded-csv-only
        """

# =====================================================================
# Rule: iqtree_pruned_ha (Runs ModelFinder Plus and generates the tree)
# =====================================================================

rule iqtree_pruned_ha:
    input:
        alignment=MAIN_PANEL_HA_POSTQC,
        partitions="data/phylogeny/HA_codons.nex"
    output:
        treefile="results/pre_beast/model_selection/HA.treefile",
        iqtree="results/pre_beast/model_selection/HA.iqtree",
        log="results/pre_beast/model_selection/HA.log"
    params:
        prefix="results/pre_beast/model_selection/HA",
        seed=RANDOM_SEED
    threads: 4
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        rm -f {params.prefix}.*
        iqtree -s {input.alignment} -spp {input.partitions} -m MFP -pre {params.prefix} -seed {params.seed} -nt {threads}
        """

PANEL_SAMPLING_DIR = "results/pre_beast/panel_sampling"

rule plot_panel_sampling_calendar:
    input:
        metadata=MAIN_PANEL_METADATA,
        context_taxa=MAIN_PANEL_CONTEXT_TAXA,
        postqc_fasta=MAIN_PANEL_HA_POSTQC,
        costa_audit=MAIN_PANEL_COSTA_AUDIT,
    output:
        tsv=f"{PANEL_SAMPLING_DIR}/panel_sampling_by_month.tsv",
        png=f"{PANEL_SAMPLING_DIR}/flu_costa_calendar.png",
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/plot_panel_sampling_calendar.py \
            --metadata {input.metadata} \
            --context-taxa {input.context_taxa} \
            --postqc-fasta {input.postqc_fasta} \
            --costa-audit {input.costa_audit} \
            --out-tsv {output.tsv} \
            --out-png {output.png}
        """

# =====================================================================
# Rule: plot_8_segments_composite
# =====================================================================

rule plot_8_segments_composite:
    input:
        trees=expand(f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/lineage/H5N1_{{segment}}_fast.treefile", segment=PHYLO_SEGMENTS),
        metadata=MAIN_PANEL_METADATA
    output:
        png="figures/main_panel_8_segments_collapsed.png"
    conda:
        "../envs/r.yml"
    params:
        outgroup=config.get("outgroup_root_sample", "EPI_ISL_18133416")
    shell:
        r"""
        Rscript code/segment_analysis/plot_8_segments.R \
            {input.metadata} \
            "{params.outgroup}" \
            {output.png} \
            {input.trees}
        """

# =====================================================================
# Rule: flumut_update_db
# =====================================================================
rule flumut_update_db:
    output:
        done="resources/flumut_update.done"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        flumut --update
        touch {output.done}
        """

# =====================================================================
# Rule: prepare_flumut_fasta
# =====================================================================
rule prepare_flumut_fasta:
    input:
        fastas=expand(f"{DATA_PHYLOGENY}/by_segment_qc_filtered/H5N1_{{segment}}.fasta", segment=PHYLO_SEGMENTS)
    output:
        fasta="data/phylogeny/flu_mut/flumut_input.fasta"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/build_flumut_fasta.py {output.fasta} {input.fastas}
        """

# =====================================================================
# Rule: run_flumut
# =====================================================================
rule run_flumut:
    input:
        fasta="data/phylogeny/flu_mut/flumut_input.fasta",
        db_done="resources/flumut_update.done"
    output:
        markers="results/phylogeny/flu_mut/flumut_report_markers.tsv",
        mutations="results/phylogeny/flu_mut/flumut_report_mutations.tsv"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        flumut -m {output.markers} -M {output.mutations} --relaxed {input.fasta}
        """
