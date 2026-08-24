rule plot_critical_mutations:
    input:
        mutations_tsv="results/phylogeny/flu_mut/flumut_report_mutations.tsv",
        metadata="metadata/H5N1_context.csv"
    output:
        plot="figures/main_panel_critical_mutations.png",
        marker_list="results/phylogeny/flu_mut/critical_markers.txt"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        "python code/segment_analysis/plot_critical_mutations.py "
        "--mutations_tsv {input.mutations_tsv} "
        "--metadata {input.metadata} "
        "--out_png {output.plot} "
        "--out_list {output.marker_list}"

rule plot_private_mutations:
    input:
        metadata="metadata/H5N1_context.csv",
        nextclade_ha="results/qc_metrics/nextclade/HA/nextclade_report.csv"
    output:
        plot="figures/main_panel_private_mutations.png"
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        "python code/segment_analysis/plot_private_mutations.py "
        "--nextclade_ha {input.nextclade_ha} "
        "--metadata {input.metadata} "
        "--out_png {output.plot}"
