PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
CODON_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA"]
SIMPLE_SEGMENTS = ["NS", "MP"]
RANDOM_SEED = config.get("random_seed", 39809473)
MAX_THREADS = int(config.get("max_threads", 20))
MAFFT_THREADS = int(config.get("mafft_threads", MAX_THREADS))
SEGMENT_RAXML_THREADS = int(config.get("raxml_segment_threads", MAX_THREADS))
FULL_RAXML_THREADS = int(config.get("raxml_full_concat_threads", MAX_THREADS))
RF_RAXML_THREADS = int(config.get("raxml_rf_threads", MAX_THREADS))

FULL_CONCAT_FASTA = "data/final/H5N1_final_beast.fasta"
FULL_CONCAT_ALIGNMENT = "data/phylogeny/aligned/H5N1_full_concat_beast.mafft"
FULL_CONCAT_PARTITIONS = "data/phylogeny/H5N1_full_concat_beast.partitions"
FULL_CONCAT_PREFIX = "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast"


rule split_h5n1_final_by_segment:
	input:
		final_fasta="data/final/H5N1_final.fasta"
	output:
		# produce one fasta file per segment plus a summary CSV
		expand("data/phylogeny/by_segment/H5N1_{segment}.fasta", segment=PHYLO_SEGMENTS),
		"data/phylogeny/by_segment_summary.csv"
	shell:
		r"""
		python code/01_ml_trees/split_final_fasta_by_segment.py \
			--input-fasta {input.final_fasta} \
			--output-dir data/phylogeny/by_segment \
			--summary-csv data/phylogeny/by_segment_summary.csv
		"""


rule mafft_align_per_segment:
	input:
		fasta="data/phylogeny/by_segment/H5N1_{segment}.fasta"
	output:
		aligned="data/phylogeny/aligned/H5N1_{segment}.mafft"
	threads: MAFFT_THREADS
	conda:
		"../envs/ml_per_segment.yml"
	shell:
		r"""
		mkdir -p data/phylogeny/aligned
		mafft --auto --thread {threads} {input.fasta} > {output.aligned}
		"""


rule align_all_segments:
	input:
		expand("data/phylogeny/aligned/H5N1_{segment}.mafft", segment=PHYLO_SEGMENTS)


rule build_segment_codon_partitions:
	input:
		alignment="data/phylogeny/aligned/H5N1_{segment}.mafft"
	output:
		partitions="data/phylogeny/codon_partitions/H5N1_{segment}.codon.partitions"
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
		alignments=expand("data/phylogeny/aligned/H5N1_{segment}.mafft", segment=PHYLO_SEGMENTS),
		segment_trees=expand(
			"results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
			segment=PHYLO_SEGMENTS,
		)
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


rule raxml_ng_tree_per_segment_codon:
	input:
		alignment="data/phylogeny/aligned/H5N1_{segment}.mafft",
		partitions="data/phylogeny/codon_partitions/H5N1_{segment}.codon.partitions"
	output:
		best_tree="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.bestTreeCollapsed",
		support_tbe="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE"
	params:
		prefix=lambda wildcards: f"results/phylogeny/raxml/{wildcards.segment}/H5N1_{wildcards.segment}",
		extra=lambda wildcards: f"--seed {RANDOM_SEED} --bs-metric tbe --force perf_threads --tree 'pars{{20}},rand{{20}}'"
	threads: SEGMENT_RAXML_THREADS
	wildcard_constraints:
		segment="PB2|PB1|PA|HA|NP|NA"
	conda:
		"../envs/ml_per_segment.yml"
	shell:
		r"""
		mkdir -p results/phylogeny/raxml/{wildcards.segment}
		raxml-ng \
			--all \
			--msa {input.alignment} \
			--model {input.partitions} \
			--prefix {params.prefix} \
			--threads {threads} \
			{params.extra}
		"""


rule raxml_ng_tree_per_segment_simple:
	input:
		alignment="data/phylogeny/aligned/H5N1_{segment}.mafft"
	output:
		best_tree="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.bestTreeCollapsed",
		support_tbe="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE"
	params:
		prefix=lambda wildcards: f"results/phylogeny/raxml/{wildcards.segment}/H5N1_{wildcards.segment}",
		model="GTR+G",
		extra=lambda wildcards: f"--seed {RANDOM_SEED} --bs-metric tbe --force perf_threads --tree 'pars{{20}},rand{{20}}'"
	threads: SEGMENT_RAXML_THREADS
	wildcard_constraints:
		segment="NS|MP"
	conda:
		"../envs/ml_per_segment.yml"
	shell:
		r"""
		mkdir -p results/phylogeny/raxml/{wildcards.segment}
		raxml-ng \
			--all \
			--msa {input.alignment} \
			--model {params.model} \
			--prefix {params.prefix} \
			--threads {threads} \
			{params.extra}
		"""


rule raxml_trees_all_segments:
	input:
		expand(
			"results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE",
			segment=PHYLO_SEGMENTS
		)


rule raxml_ng_tree_full_concat_beast:
	input:
		alignment=FULL_CONCAT_ALIGNMENT,
		partitions=FULL_CONCAT_PARTITIONS
	output:
		best_tree=f"{FULL_CONCAT_PREFIX}.raxml.bestTreeCollapsed",
		support_tbe=f"{FULL_CONCAT_PREFIX}.raxml.supportTBE"
	params:
		prefix=FULL_CONCAT_PREFIX,
		extra=lambda wildcards: f"--seed {RANDOM_SEED} --bs-metric tbe --force perf_threads --tree 'pars{{30}},rand{{30}}'"
	threads: FULL_RAXML_THREADS
	conda:
		"../envs/ml_per_segment.yml"
	shell:
		r"""
		mkdir -p results/phylogeny/raxml/full_concat
		raxml-ng \
			--all \
			--msa {input.alignment} \
			--model {input.partitions} \
			--prefix {params.prefix} \
			--threads {threads} \
			{params.extra}
		"""


rule raxml_tree_full_concat:
	input:
		f"{FULL_CONCAT_PREFIX}.raxml.supportTBE"
