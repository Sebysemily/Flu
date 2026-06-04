RTT_DIR = "results/pre_beast/rtt_HA"
PRE_BEAST_DATES = f"{RTT_DIR}/dates_from_headers.tsv"
PRE_BEAST_ROOT_TO_TIP_OUTLIERS = f"{RTT_DIR}/outliers.tsv"
PRE_BEAST_ROOT_TO_TIP_LOG = f"{RTT_DIR}/treetime_clock.log"
PRE_BEAST_ROOT_TO_TIP_DONE = f"{RTT_DIR}/treetime_clock.done"

rule build_treetime_dates_subsets:
    input:
        metadata="metadata/H5N1_context.csv",
    output:
        dates=PRE_BEAST_DATES,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.dates})
        python code/02_Beast/build_treetime_dates.py \
            --metadata {input.metadata} \
            --out {output.dates}
        """

rule run_root_to_tip_subsets_HA:
    input:
        tree="results/pre_beast/model_selection/HA.treefile",
        alignment=MAIN_PANEL_HA_POSTQC,
        dates=PRE_BEAST_DATES,
    output:
        outliers=PRE_BEAST_ROOT_TO_TIP_OUTLIERS,
        log=PRE_BEAST_ROOT_TO_TIP_LOG,
        done=PRE_BEAST_ROOT_TO_TIP_DONE,
    conda:
        "../envs/02_pre_beast.yml"
    params:
        clock_filter=config.get("treetime_parameters", {}).get("clock_filter", 3.0),
    shell:
        r"""
        python code/02_Beast/run_treetime_clock.py \
            --tree {input.tree} \
            --aln {input.alignment} \
            --dates {input.dates} \
            --outdir $(dirname {output.log}) \
            --clock-filter {params.clock_filter} \
            --log {output.log} \
            --done {output.done}
        """

rule build_beast_metadata:
    input:
        metadata="metadata/H5N1_context.csv",
        alignment=MAIN_PANEL_HA_POSTQC
    output:
        tsv="metadata/beast/metadata_beast.tsv"
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/build_beast_metadata.py \
            --metadata {input.metadata} \
            --aln {input.alignment} \
            --out {output.tsv}
        """
