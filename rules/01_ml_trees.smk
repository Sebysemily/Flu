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
        metadata=FILTRADO_CSV
    output:
        filtered=f"{DATA_PHYLOGENY}/by_segment_qc_filtered/H5N1_HA.fasta",
        discarded_csv=f"{RESULTS_QC_METRICS}/nextclade/HA_eliminated.csv"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/build_inputs/filter_qc_nextclade.py \
            --input-alignment {input.alignment} \
            --report {input.report} \
            --output-alignment {output.filtered} \
            --metadata {input.metadata} \
            --discarded-csv {output.discarded_csv}
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
        alignment=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft"
    output:
        filtered=f"{PROCESSED_ALIGNMENTS_QC_FILTERED}/H5N1_{{segment}}.mafft"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/trim_and_filter_mafft.py \
            --input {input.alignment} \
            --output {output.filtered} \
            --max-divergence 0.1
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
# Rule: run_iqtree_codon_segment
# =====================================================================

RESULTS_PHYLOGENY = config.get("results_phylogeny", "results/phylogeny")

rule run_iqtree_codon_segment:
    input:
        alignment=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta",
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
        alignment=f"{MAIN_PANEL_DIR}/H5N1_{{segment}}.fasta"
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
        heatmap_rds="figures/build_inputs/flu_lineage.rds",
        ha_tree=f"{RESULTS_PHYLOGENY}/iq-tree/HA/lineage/H5N1_HA_fast.treefile",
        pb2_tree=f"{RESULTS_PHYLOGENY}/iq-tree/PB2/lineage/H5N1_PB2_fast.treefile",
        panel_metadata=MAIN_PANEL_METADATA,
    output:
        png=f"{LINEAGE_FIGURES_DIR}/HA_PB2_lineage_composite_preview.png",
    params:
        figures_dir=LINEAGE_FIGURES_DIR,
        max_tips=60,
        ribbon_segment="HA",
        outgroup_sample=config.get("outgroup_root_sample", ""),
    conda:
        "../envs/r.yml"
    shell:
        r"""
        ulimit -s unlimited || true
        mkdir -p {params.figures_dir}
        Rscript code/segment_analysis/plot_composite_lineage.R \
            {input.heatmap_rds} \
            {input.ha_tree} \
            {input.pb2_tree} \
            {input.panel_metadata} \
            {output.png} \
            {params.max_tips} \
            {params.ribbon_segment} \
            "{params.outgroup_sample}"
        test -s {output.png}
        """

rule segment_analysis_composite_lineage:
    input:
        heatmap_rds="figures/build_inputs/flu_lineage.rds",
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
            {input.heatmap_rds} \
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
PRUNED_HA_DIR = "data/phylogeny/main_panel/pre-alignment"
MAIN_PANEL_CONTEXT_TAXA = f"{PRUNED_HA_DIR}/context_taxa.txt"
MAIN_PANEL_FILTERED_AUDIT = f"{PRUNED_HA_DIR}/pruning.audit.tsv"

rule augur_filter_context:
    input:
        metadata=MAIN_PANEL_METADATA,
    output:
        context_taxa=MAIN_PANEL_CONTEXT_TAXA,
        include_list=f"{PRUNED_HA_DIR}/augur_include.txt",
    conda:
        "../envs/01_ml_trees.yml"
    params:
        outgroup_root_sample=config.get("outgroup_root_sample", "EPI_ISL_18133416"),
        seed=RANDOM_SEED
    shell:
        r"""
        python -c "
import pandas as pd
df = pd.read_csv('{input.metadata}', dtype=str)
df['date'] = df['collection_date']

# Select up to 5 american anchors per month/year to protect
anchors = df[df['expected_role'] == 'american_anchor'].copy()
def get_year_month(d):
    if pd.isna(d):
        return 'unknown'
    parts = d.split('-')
    if len(parts) >= 2:
        return parts[0] + '-' + parts[1]
    return parts[0]
anchors['year_month'] = anchors['collection_date'].apply(get_year_month)

protected_anchors = []
for ym, gp in anchors.groupby('year_month'):
    protected_anchors.extend(gp.head(5)['file_name'].tolist())

# Write all protected strains to include_list (Ecuador + monthly anchors + outgroup)
flu_strains = df[df['expected_role'].str.startswith('flu_', na=False)]['file_name'].tolist()
final_include = list(set(['{params.outgroup_root_sample}'] + flu_strains + protected_anchors))

with open('{output.include_list}', 'w') as f:
    for strain in sorted(final_include):
        f.write(strain + '\n')

# Prepare temporary metadata file with is_subsamplable column
df['is_subsamplable'] = df['expected_role'].apply(
    lambda r: 'yes' if r in ['american_anchor', 'regional_context'] else 'no'
)
df.to_csv('{input.metadata}.tmp.csv', index=False)

# Create weights TSV file with integer weights
weights_content = 'expected_role\tweight\nregional_context\t3\namerican_anchor\t1\ndefault\t1\n'
with open('{output.context_taxa}.weights.tsv', 'w') as f:
    f.write(weights_content)
"

        # Execute single weighted augur filter call
        augur filter \
            --metadata {input.metadata}.tmp.csv \
            --metadata-id-columns file_name \
            --exclude-where is_subsamplable=no \
            --include {output.include_list} \
            --group-by expected_role \
            --group-by-weights {output.context_taxa}.weights.tsv \
            --subsample-max-sequences 400 \
            --subsample-seed {params.seed} \
            --output-strains {output.context_taxa}

        # Cleanup temporary files
        rm -f {input.metadata}.tmp.csv {output.context_taxa}.weights.tsv
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
        fasta=f"{PRUNED_HA_DIR}/H5N1_HA_pruned.fasta",
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
        fasta=f"{PRUNED_HA_DIR}/H5N1_HA_pruned.fasta"
    output:
        trimmed=f"{PRUNED_HA_DIR}/H5N1_HA_pruned.mafft"
    conda:
        "../envs/01_ml_trees.yml"
    threads: MAX_THREADS
    shell:
        r"""
        mafft --auto --thread {threads} {input.fasta} > {output.trimmed}
        """

rule trim_and_filter_pruned_ha:
    input:
        alignment=f"{PRUNED_HA_DIR}/H5N1_HA_pruned.mafft"
    output:
        filtered=f"{MAIN_PANEL_DIR}/H5N1_HA.fasta"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/01_ml_trees/trim_and_filter_mafft.py \
            --input {input.alignment} \
            --output {output.filtered} \
            --max-divergence 0.05
        """

# =====================================================================
# Rule: iqtree_pruned_ha
# =====================================================================

rule iqtree_pruned_ha:
    input:
        alignment=f"{MAIN_PANEL_DIR}/H5N1_HA.fasta"
    output:
        treefile=f"{RESULTS_PHYLOGENY}/iq-tree/HA/HA.treefile"
    params:
        prefix=f"{RESULTS_PHYLOGENY}/iq-tree/HA/workdir/HA",
        seed=RANDOM_SEED,
        bootstrap=1000
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
