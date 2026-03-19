PHYLO_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


rule split_h5n1_final_by_segment:
	input:
		final_fasta="data/final/H5N1_final.fasta"
	output:
		directory("data/phylogeny/by_segment"),
		"data/phylogeny/by_segment_summary.csv"
	shell:
		r"""
		python code/split_final_fasta_by_segment.py \
			--input-fasta {input.final_fasta} \
			--output-dir data/phylogeny/by_segment \
			--summary-csv data/phylogeny/by_segment_summary.csv
		"""


rule mafft_align_per_segment:
	input:
		fasta="data/phylogeny/by_segment/H5N1_{segment}.fasta"
	output:
		aligned="data/phylogeny/aligned/H5N1_{segment}.mafft"
	threads: 4
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


rule raxml_ng_tree_per_segment:
	input:
		alignment="data/phylogeny/aligned/H5N1_{segment}.mafft"
	output:
		best_tree="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.bestTreeCollapsed",
		support_fbp="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportFBP",
		support_tbe="results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportTBE"
	params:
		prefix=lambda wildcards: f"results/phylogeny/raxml/{wildcards.segment}/H5N1_{wildcards.segment}",
		model="GTR+G",
		bs_trees=300,
		extra="--bs-metric fbp,tbe"
	threads: 12
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
			--bs-trees {params.bs_trees} \
			{params.extra}
		"""


rule raxml_trees_all_segments:
	input:
		expand(
			"results/phylogeny/raxml/{segment}/H5N1_{segment}.raxml.supportFBP",
			segment=PHYLO_SEGMENTS
		)
