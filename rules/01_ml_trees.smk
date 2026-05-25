PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
SIMPLE_SEGMENTS = ["NS", "MP"]
NP_MP_NS_SEGMENTS = ["NP", "MP", "NS"]
NP_MP_NS_CODON_SEGMENTS = ["NP"]
RANDOM_SEED = config.get("random_seed", 39809473)
MAX_THREADS = int(config.get("max_threads", 18))

DATA_COMBINED_CONTEXT_EC = config.get("data_combined_context_ecuador", "data/final")
DATA_PHYLOGENY = config.get("data_phylogeny", "data/phylogeny")
RESULTS_PHYLOGENY = config.get("results_phylogeny", "results/phylogeny")
DATA_PRE_BEAST = config.get("data_pre_beast", "data/pre_beast")
PROCESSED_ALIGNMENTS_CODON_AWARE = "data/processed_alignments/codon_aware"
QC_METRICS_DIR = "results/qc_metrics/processed_alignments/codon_aware"

FULL_CONCAT_FASTA = f"{DATA_COMBINED_CONTEXT_EC}/H5N1_final_beast.fasta"
FULL_CONCAT_ALIGNMENT = f"{DATA_PHYLOGENY}/aligned/H5N1_full_concat_beast.mafft"
FULL_CONCAT_PARTITIONS = f"{DATA_PHYLOGENY}/H5N1_full_concat_beast.partitions"
FULL_CONCAT_PREFIX = f"{RESULTS_PHYLOGENY}/iq-tree/full_concat/H5N1_full_concat_beast"
NP_MP_NS_ALIGNMENT = f"{DATA_PHYLOGENY}/aligned/H5N1_NP_MP_NS.mafft"
NP_MP_NS_PARTITIONS = f"{DATA_PHYLOGENY}/H5N1_NP_MP_NS.partitions"
NP_MP_NS_PREFIX = f"{RESULTS_PHYLOGENY}/iq-tree/np_mp_ns/H5N1_NP_MP_NS"



CONCAT_TREE = "results/phylogeny/iq-tree/full_concat/H5N1_full_concat_beast.treefile"
CONCAT_ALIGNMENT = "data/phylogeny/aligned/H5N1_full_concat_beast.mafft"

PRE_BEAST_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

BEAST_PRE_DATA_DIR = "data/pre_beast/panels"
BEAST_EXPORT_DIR = "data/beast"
BEAST_FINAL_DATA_DIR = f"{BEAST_EXPORT_DIR}/final_panel_segment"
BEAST_RAW_FINAL_SEGMENT_DIR = f"{BEAST_PRE_DATA_DIR}/raw_segment_subset"
BEAST_PRE_RESULTS_DIR = "results/beast_pre"
BEAST_PRE_QC_DIR = f"{BEAST_PRE_RESULTS_DIR}/qc_validation"
BEAST_SEGMENT_QC_DIR = f"{BEAST_PRE_QC_DIR}/segment_alignment"

BEAST_PANEL_TAXA = f"{BEAST_PRE_DATA_DIR}/panel_main_taxa.tsv"
BEAST_PANEL_FILTERED_TAXA = f"{BEAST_PRE_DATA_DIR}/panel_main_taxa.filtered.tsv"
BEAST_PANEL_RTT_FILTERED_TAXA = f"{BEAST_EXPORT_DIR}/panel_main_taxa.final.tsv"

BEAST_FILTERED_SUBSET_ALIGNMENT = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.filtered.fasta"
BEAST_FILTERED_SUBSET_TREE = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.filtered.nwk"
BEAST_FILTERED_SUBSET_AUDIT = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.filtered.audit.tsv"
BEAST_FINAL_SUBSET_ALIGNMENT = f"{BEAST_EXPORT_DIR}/panel_main_concat.final.fasta"
BEAST_FINAL_SUBSET_AUDIT = f"{BEAST_EXPORT_DIR}/panel_main_concat.final.audit.tsv"



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


rule mafft_align_per_segment:
	input:
		fasta=f"{DATA_PHYLOGENY}/by_segment/H5N1_{{segment}}.fasta"
	output:
		aligned=temp(f"{DATA_PHYLOGENY}/aligned/H5N1_{{segment}}.mafft")
	threads: MAX_THREADS
	conda:
		"../envs/01_ml_trees.yml"
	shell:
		r"""
		mkdir -p {DATA_PHYLOGENY}/aligned
		mafft --auto --thread {threads} {input.fasta} > {output.aligned}
		"""


rule codon_qc_trim_per_segment:
	input:
		aligned=f"{DATA_PHYLOGENY}/aligned/H5N1_{{segment}}.mafft"
	output:
		trimmed=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft",
		qc_csv=f"{QC_METRICS_DIR}/H5N1_{{segment}}.qc_metrics.csv"
	shell:
		r"""
		python code/build_inputs/codon_qc_trim.py \
			--alignment {input.aligned} \
			--segment {wildcards.segment} \
			--output {output.trimmed} \
			--qc-summary {output.qc_csv}
		"""


rule build_segment_codon_partitions:
	input:
		alignment=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft"
	output:
		partitions=f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.codon.partitions"
	wildcard_constraints:
		segment="PB2|PB1|PA|HA|NP|NA"
	shell:
		r"""
		python code/01_ml_trees/build_single_segment_codon_partition.py \
			--alignment {input.alignment} \
			--segment {wildcards.segment} \
			--output {output.partitions}
		"""


rule concat_aligned_segments_with_partitions:
	input:
		alignments=expand(f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS)
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


rule concat_np_mp_ns_with_partitions:
	input:
		alignments=expand(f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft", segment=NP_MP_NS_SEGMENTS)
	output:
		aligned=NP_MP_NS_ALIGNMENT,
		partitions=NP_MP_NS_PARTITIONS
	params:
		segment_order=",".join(NP_MP_NS_SEGMENTS),
		codon_segments=",".join(NP_MP_NS_CODON_SEGMENTS)
	shell:
		r"""
		python code/01_ml_trees/build_concat_codon_partitions.py \
			--segment-order {params.segment_order} \
			--codon-segments {params.codon_segments} \
			--output-alignment {output.aligned} \
			--output-partitions {output.partitions} \
			{input.alignments}
		"""



rule iqtree_np_mp_ns:
	input:
		alignment=NP_MP_NS_ALIGNMENT,
		partitions=NP_MP_NS_PARTITIONS
	output:
		treefile=f"{NP_MP_NS_PREFIX}.treefile"
	params:
		prefix=NP_MP_NS_PREFIX,
		seed=RANDOM_SEED
	threads: MAX_THREADS
	conda:
		"../envs/01_ml_trees.yml"
	shell:
		r"""
		mkdir -p {RESULTS_PHYLOGENY}/iq-tree/np_mp_ns
		rm -f {params.prefix}.ckp.gz
		iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -m GTR+G -fast
		"""


rule iqtree_fast_full_concat:
	input:
		alignment=FULL_CONCAT_ALIGNMENT,
		partitions=FULL_CONCAT_PARTITIONS
	output:
		treefile=f"{FULL_CONCAT_PREFIX}.treefile"
	params:
		prefix=FULL_CONCAT_PREFIX,
		seed=RANDOM_SEED
	threads: MAX_THREADS
	conda:
		"../envs/01_ml_trees.yml"
	shell:
		r"""
		mkdir -p {RESULTS_PHYLOGENY}/iq-tree/full_concat
		rm -f {params.prefix}.ckp.gz
		iqtree -s {input.alignment} -spp {input.partitions} -pre {params.prefix} -seed {params.seed} -nt {threads} -m GTR+G -fast
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
	conda:
		"../envs/01_ml_trees.yml"
	params:
		n_closest=config.get("beast_panel", {}).get("n_closest", 200),
		max_distance=config.get("beast_panel", {}).get("max_distance", 0.08),
		protect_anchors_count=config.get("beast_panel", {}).get("protect_anchors_count", 10),
		temp_audit=f"{BEAST_PRE_DATA_DIR}/panel_distance_audit.tsv",
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
			--n-closest {params.n_closest} \
			--max-distance {params.max_distance} \
			--protect-anchors-count {params.protect_anchors_count}
		"""


# ---------------------------------------------------------------------------
# Alias rules — input-only rules with no shell/output.
# Purpose: provide clean named CLI targets (e.g. `snakemake align_all_segments`)
# and cleaner rulegraph nodes instead of embedding full rule names in DAG edges.
# ---------------------------------------------------------------------------

rule align_all_segments:
	input:
		expand(f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft", segment=PHYLO_SEGMENTS)





rule iqtree_tree_full_concat:
	input:
		f"{FULL_CONCAT_PREFIX}.treefile"


rule iqtree_tree_np_mp_ns:
	input:
		f"{NP_MP_NS_PREFIX}.treefile"


