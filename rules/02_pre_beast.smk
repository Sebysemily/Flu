
PRE_BEAST_FINAL_SEGMENT_FASTAS = expand(
    f"{FINAL_PANEL_DIR}/H5N1_{{segment}}.fasta",
    segment=PRE_BEAST_SEGMENTS,
)

RTT_QC_DIR = "results/qc_metrics/phylogeny/rtt"

PRE_BEAST_RTT_EXCLUSIONS = f"{RTT_QC_DIR}/rtt_outlier_exclusions.tsv"
PRE_BEAST_RTT_EXCLUSIONS_SUMMARY = f"{RTT_QC_DIR}/rtt_outlier_exclusions.summary.tsv"

PRE_BEAST_DATES = f"{RTT_DIR}/dates_from_headers.tsv"
PRE_BEAST_ROOT_TO_TIP_OUTLIERS = f"{RTT_DIR}/outliers.tsv"
PRE_BEAST_ROOT_TO_TIP_LOG = f"{RTT_DIR}/treetime_clock.log"
PRE_BEAST_ROOT_TO_TIP_DONE = f"{RTT_DIR}/treetime_clock.done"
BEAST_FINAL_DATES = f"{DATA_BAYESIAN}/panel_main_dates.final.tsv"

BEAST_CONFIG = config.get("beast", {})
BEAST_LOG_EVERY = int(BEAST_CONFIG.get("log_every", 2000))
BEAST_TREE_EVERY = int(BEAST_CONFIG.get("tree_every", BEAST_LOG_EVERY))
BEAST_ECHO_EVERY = int(BEAST_CONFIG.get("echo_every", BEAST_LOG_EVERY))

BEAST_XML_DIR = f"{RESULTS_BAYESIAN}/xml"
BEAST_TEMPLATE_XMLS = {
    "strict_constant": "template_beast/concat_strict_constant.base.xml",
    "strict_constant_lugar": "template_beast/concat_strict_const_lugar.xml",
    "ucln_constant": "template_beast/concat_ucln_constant.base.xml",
    "strict_exp": "template_beast/concat_strict_exp.base.xml",
    "ucln_exp": "template_beast/concat_ucln_exp.base.xml",
}
BEAST_RUN_SCENARIOS = list(BEAST_TEMPLATE_XMLS.keys())
BEAST_CHAIN_LENGTH_CONFIG = BEAST_CONFIG.get("chain_length", 20000000)
if isinstance(BEAST_CHAIN_LENGTH_CONFIG, dict):
    BEAST_CHAIN_LENGTHS = {
        scenario: int(BEAST_CHAIN_LENGTH_CONFIG.get(scenario, 20000000))
        for scenario in BEAST_RUN_SCENARIOS
    }
else:
    BEAST_CHAIN_LENGTHS = {
        scenario: int(BEAST_CHAIN_LENGTH_CONFIG)
        for scenario in BEAST_RUN_SCENARIOS
    }
BEAST_PREPARED_XMLS = {
    scenario: f"{BEAST_XML_DIR}/{scenario}.xml"
    for scenario in BEAST_TEMPLATE_XMLS
}
PREPARED_BEAST_XMLS = list(BEAST_PREPARED_XMLS.values())













rule build_treetime_dates_subsets:
    input:
        alignment=FINAL_ALIGNMENT,
    output:
        dates=PRE_BEAST_DATES,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.dates})
        python code/02_Beast/build_treetime_dates.py \
            --aln {input.alignment} \
            --out {output.dates}
        """


rule run_root_to_tip_subsets:
    input:
        tree=f"{RESULTS_PHYLOGENY}/iq-tree/subset_fast/panel_main_concat.final.treefile",
        alignment=FINAL_ALIGNMENT,
        dates=PRE_BEAST_DATES,
    output:
        outliers=PRE_BEAST_ROOT_TO_TIP_OUTLIERS,
        log=PRE_BEAST_ROOT_TO_TIP_LOG,
        done=PRE_BEAST_ROOT_TO_TIP_DONE,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/run_treetime_clock.py \
            --tree {input.tree} \
            --aln {input.alignment} \
            --dates {input.dates} \
            --outdir $(dirname {output.log}) \
            --log {output.log} \
            --done {output.done}
        """


rule filter_rtt_outliers_subsets:
    input:
        panel_tsv=BEAST_PANEL_FILTERED_TAXA,
        rtt_outliers=PRE_BEAST_ROOT_TO_TIP_OUTLIERS,
        dates_in=PRE_BEAST_DATES,
    output:
        filtered_panel=BEAST_PANEL_RTT_FILTERED_TAXA,
        dates_out=BEAST_FINAL_DATES,
        exclusions=PRE_BEAST_RTT_EXCLUSIONS,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/filter_beast_panel_by_rtt_outliers.py \
            --panel-taxa {input.panel_tsv} \
            --rtt-outliers {input.rtt_outliers} \
            --dates-in {input.dates_in} \
            --dates-out {output.dates_out} \
            --filtered-panel-out {output.filtered_panel} \
            --exclusions-out {output.exclusions}
        """


rule generate_bayesian_alignment:
    input:
        alignment=FINAL_ALIGNMENT,
        taxa=BEAST_PANEL_RTT_FILTERED_TAXA
    output:
        alignment=BAYESIAN_ALIGNMENT
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/02_Beast/subset_alignment_by_taxa.py \
            --alignment {input.alignment} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment}
        """


rule prepare_beast_run_xml:
    input:
        template_xml=lambda wildcards: BEAST_TEMPLATE_XMLS[wildcards.scenario],
        panel_tsv=BEAST_PANEL_RTT_FILTERED_TAXA,
        dates=BEAST_FINAL_DATES,
        subset_alignment=BAYESIAN_ALIGNMENT,
    output:
        xml=f"{BEAST_XML_DIR}/{{scenario}}.xml",
    wildcard_constraints:
        scenario="|".join(BEAST_TEMPLATE_XMLS.keys()),
    params:
        output_prefix=lambda wildcards: f"results/phylogeny/bayesian/runs/{wildcards.scenario}/{wildcards.scenario}",
        chain_length=lambda wildcards: BEAST_CHAIN_LENGTHS[wildcards.scenario],
        log_every=BEAST_LOG_EVERY,
        tree_every=BEAST_TREE_EVERY,
        echo_every=BEAST_ECHO_EVERY,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/prepare_beast_run_xml.py \
            --template-xml {input.template_xml} \
            --scenario-name {wildcards.scenario} \
            --output-xml {output.xml} \
            --output-prefix {params.output_prefix} \
            --panel-taxa {input.panel_tsv} \
            --chain-length {params.chain_length} \
            --log-every {params.log_every} \
            --tree-every {params.tree_every} \
            --echo-every {params.echo_every}
        """

