CONCAT_TREE = "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE"
CONCAT_ALIGNMENT = "data/phylogeny/aligned/H5N1_full_concat_beast.mafft"
BEAST_SUMMARY = "data/final/H5N1_final_beast_summary.csv"
CONTEXT_META = config.get("context_metadata_tsv", "config/final_metadata_50_per_country_isolates.tsv")
ECUADOR_DATE_SOURCE = config.get("ecuador_date_source", "reception")
RANDOM_SEED = config.get("random_seed", 39809473)


rule build_beast_panels:
    input:
        tree=CONCAT_TREE,
        beast_summary=BEAST_SUMMARY,
        context_metadata=CONTEXT_META,
    output:
        panel_main="data/beast_pre/panels/panel_main_taxa.tsv",
        audit="data/beast_pre/panels/panel_selection_audit.tsv",
    conda:
        "../envs/02_beast.yml"
    params:
        max_cluster_dist=config.get("beast_pre", {}).get("max_cluster_dist", 0.08),
        n_per_country=config.get("beast_pre", {}).get("n_per_country", 4),
        n_total=config.get("beast_pre", {}).get("n_total", 60),
    shell:
        r"""
        mkdir -p data/beast_pre/panels
        python code/build_beast_panels.py \
            --tree {input.tree} \
            --beast-summary {input.beast_summary} \
            --context-metadata {input.context_metadata} \
            --panel-main-out {output.panel_main} \
            --audit-out {output.audit} \
            --max-cluster-dist {params.max_cluster_dist} \
            --n-per-country {params.n_per_country} \
            --n-total {params.n_total}
        """


rule subset_alignment_and_prune_tree:
    input:
        alignment=CONCAT_ALIGNMENT,
        tree=CONCAT_TREE,
        taxa="data/beast_pre/panels/panel_main_taxa.tsv",
    output:
        alignment="data/beast_pre/panels/panel_main.subset.mafft",
        tree="data/beast_pre/panels/panel_main.subset.nwk",
        audit="data/beast_pre/panels/panel_main.subset.audit.tsv",
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        python code/subset_alignment_and_prune_tree.py \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment} \
            --out-tree {output.tree} \
            --audit {output.audit}
        """


rule run_root_to_tip:
    input:
        tree="data/beast_pre/panels/panel_main.subset.nwk",
        alignment="data/beast_pre/panels/panel_main.subset.mafft",
        flu_filtrado="config/flu_filtrado.csv",
        context_metadata=CONTEXT_META,
    output:
        log="data/beast_pre/root_to_tip/treetime_clock.log",
        done="data/beast_pre/root_to_tip/treetime_clock.done",
    params:
        ecuador_date_source=ECUADOR_DATE_SOURCE,
        random_seed=RANDOM_SEED,
    conda:
        "../envs/02_beast.yml"
    shell:
        r"""
        python code/run_treetime_clock.py \
            --tree {input.tree} \
            --aln {input.alignment} \
            --flu-filtrado-csv {input.flu_filtrado} \
            --ecuador-date-source {params.ecuador_date_source} \
            --context-metadata {input.context_metadata} \
            --seed {params.random_seed} \
            --outdir data/beast_pre/root_to_tip \
            --log {output.log} \
            --done {output.done}
        """
