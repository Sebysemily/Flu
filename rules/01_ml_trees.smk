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

BEAST_PRE_DATA_DIR = f"{DATA_PRE_BEAST}/panels"
BEAST_EXPORT_DIR = DATA_BEAST
BEAST_FINAL_DATA_DIR = f"{BEAST_EXPORT_DIR}/final_panel_segment"
BEAST_RAW_FINAL_SEGMENT_DIR = f"{BEAST_PRE_DATA_DIR}/raw_segment_subset"
BEAST_PRE_RESULTS_DIR = RESULTS_BEAST_PRE
BEAST_PRE_QC_DIR = f"{RESULTS_QC_METRICS}/beast"
BEAST_SEGMENT_QC_DIR = f"{BEAST_PRE_QC_DIR}/final_panel_segment"

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


rule concat_except_ha_with_partitions:
	input:
		alignments=expand(f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{{segment}}.mafft", segment=["PB2", "PB1", "PA", "NP", "NA", "MP", "NS"])
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


def get_iqtree_replicate_inputs(wildcards):
	if wildcards.segment == "concat_except_HA":
		return {
			"alignment": CONCAT_EXCEPT_HA_ALIGNMENT,
			"partitions": CONCAT_EXCEPT_HA_PARTITIONS
		}
	inputs = {
		"alignment": f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{wildcards.segment}.mafft"
	}
	if wildcards.segment in CODON_SEGMENTS:
		inputs["partitions"] = f"{PROCESSED_ALIGNMENTS_CODON_AWARE}/H5N1_{wildcards.segment}.codon.partitions"
	return inputs

def get_iqtree_replicate_seed(wildcards):
	return RANDOM_SEED + int(wildcards.rep)

def get_iqtree_replicate_part_arg(wildcards, input):
	if "partitions" in input.keys():
		return f"-spp {input.partitions}"
	return ""

rule run_iqtree_segment_replicate:
	input:
		unpack(get_iqtree_replicate_inputs)
	output:
		treefile=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}.treefile"
	params:
		prefix=f"{RESULTS_PHYLOGENY}/iq-tree/{{segment}}/rep{{rep}}/rep{{rep}}",
		seed=get_iqtree_replicate_seed,
		part_arg=get_iqtree_replicate_part_arg
	threads: 4
	conda:
		"../envs/01_ml_trees.yml"
	shell:
		r"""
		mkdir -p $(dirname {params.prefix})
		iqtree -s {input.alignment} {params.part_arg} -pre {params.prefix} -seed {params.seed} -nt {threads} -bb 1000 -bnni
		mv {params.prefix}.treefile {output.treefile}
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


