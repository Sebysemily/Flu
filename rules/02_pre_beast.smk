
PRE_BEAST_RAW_FINAL_SEGMENT_FASTAS = expand(
    f"{BEAST_RAW_FINAL_SEGMENT_DIR}/H5N1_{{segment}}.fasta",
    segment=PRE_BEAST_SEGMENTS,
)
PRE_BEAST_FINAL_SEGMENT_FASTAS = expand(
    f"{BEAST_FINAL_DATA_DIR}/H5N1_{{segment}}.fasta",
    segment=PRE_BEAST_SEGMENTS,
)

PRE_BEAST_SOURCE_QC_METRICS = f"{BEAST_PRE_QC_DIR}/source_panel_qc.metrics.tsv"
PRE_BEAST_SOURCE_QC_OUTLIERS = f"{BEAST_PRE_QC_DIR}/source_panel_qc.outliers.tsv"
PRE_BEAST_SOURCE_QC_SUMMARY = f"{BEAST_PRE_QC_DIR}/source_panel_qc.summary.tsv"
PRE_BEAST_SOURCE_QC_REPORT = f"{BEAST_PRE_QC_DIR}/source_panel_qc.report.md"

PRE_BEAST_SEGMENT_QC_METRICS = expand(
    f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.metrics.tsv",
    segment=PRE_BEAST_SEGMENTS,
)
PRE_BEAST_SEGMENT_QC_OUTLIERS = expand(
    f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.outliers.tsv",
    segment=PRE_BEAST_SEGMENTS,
)
PRE_BEAST_SEGMENT_QC_SUMMARIES = expand(
    f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.summary.tsv",
    segment=PRE_BEAST_SEGMENTS,
)
PRE_BEAST_SEGMENT_QC_REPORTS = expand(
    f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.report.md",
    segment=PRE_BEAST_SEGMENTS,
)
PRE_BEAST_FINAL_SEGMENT_QC_METRICS = f"{BEAST_PRE_QC_DIR}/final_segment_panel_qc.metrics.tsv"
PRE_BEAST_FINAL_SEGMENT_QC_OUTLIERS = f"{BEAST_PRE_QC_DIR}/final_segment_panel_qc.outliers.tsv"
PRE_BEAST_FINAL_SEGMENT_QC_SUMMARY = f"{BEAST_PRE_QC_DIR}/final_segment_panel_qc.summary.tsv"
PRE_BEAST_FINAL_SEGMENT_QC_REPORT = f"{BEAST_PRE_QC_DIR}/final_segment_panel_qc.report.md"

PRE_BEAST_PANEL_EXCLUSIONS = f"{BEAST_PRE_QC_DIR}/panel_exclusions.tsv"
PRE_BEAST_PANEL_EXCLUSIONS_SUMMARY = f"{BEAST_PRE_QC_DIR}/panel_exclusions.summary.tsv"

PRE_BEAST_RTT_EXCLUSIONS = f"{BEAST_PRE_QC_DIR}/rtt_outlier_exclusions.tsv"
PRE_BEAST_RTT_EXCLUSIONS_SUMMARY = f"{BEAST_PRE_QC_DIR}/rtt_outlier_exclusions.summary.tsv"

PRE_BEAST_RTT_DIR = f"{BEAST_PRE_RESULTS_DIR}/rtt"
PRE_BEAST_DATES = f"{PRE_BEAST_RTT_DIR}/dates_from_headers.tsv"
PRE_BEAST_ROOT_TO_TIP_OUTLIERS = f"{PRE_BEAST_RTT_DIR}/outliers.tsv"
PRE_BEAST_ROOT_TO_TIP_LOG = f"{PRE_BEAST_RTT_DIR}/treetime_clock.log"
PRE_BEAST_ROOT_TO_TIP_DONE = f"{PRE_BEAST_RTT_DIR}/treetime_clock.done"
BEAST_FINAL_DATES = f"{BEAST_EXPORT_DIR}/panel_main_dates.final.tsv"

NP_MP_NS_SENSITIVITY_CONFIG = config.get("np_mp_ns_sensitivity", {})
NP_MP_NS_FORCED_AMERICAN_ANCHOR = NP_MP_NS_SENSITIVITY_CONFIG.get(
    "forced_american_anchor_accession",
    "OQ968009",
)
NP_MP_NS_AMERICAN_ANCHOR_TOTAL = int(
    NP_MP_NS_SENSITIVITY_CONFIG.get("american_anchor_total", 6)
)
NP_MP_NS_SENSITIVITY_DATA_DIR = f"{DATA_PRE_BEAST}/sensitivity_np_mp_ns"
NP_MP_NS_SENSITIVITY_RESULTS_DIR = f"{RESULTS_QC_METRICS}/pre_beast/sensitivity_np_mp_ns"
NP_MP_NS_SENSITIVITY_RTT_DIR = f"{NP_MP_NS_SENSITIVITY_RESULTS_DIR}/rtt"
NP_MP_NS_SENSITIVITY_PANEL_TAXA = f"{NP_MP_NS_SENSITIVITY_DATA_DIR}/panel_taxa.tsv"
NP_MP_NS_SENSITIVITY_PANEL_AUDIT = f"{NP_MP_NS_SENSITIVITY_RESULTS_DIR}/panel_audit.tsv"
NP_MP_NS_SENSITIVITY_SUBSET_ALIGNMENT = (
    f"{NP_MP_NS_SENSITIVITY_DATA_DIR}/panel_np_mp_ns.fasta"
)
NP_MP_NS_SENSITIVITY_SUBSET_TREE = (
    f"{NP_MP_NS_SENSITIVITY_DATA_DIR}/panel_np_mp_ns.nwk"
)
NP_MP_NS_SENSITIVITY_SUBSET_AUDIT = (
    f"{NP_MP_NS_SENSITIVITY_RESULTS_DIR}/panel_np_mp_ns.audit.tsv"
)
NP_MP_NS_SENSITIVITY_DATES = f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/dates_from_headers.tsv"
NP_MP_NS_SENSITIVITY_RTT_OUTLIERS = f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/outliers.tsv"
NP_MP_NS_SENSITIVITY_RTT_LOG = f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/treetime_clock.log"
NP_MP_NS_SENSITIVITY_RTT_DONE = f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/treetime_clock.done"
NP_MP_NS_SENSITIVITY_RTT_CSV = f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/rtt.csv"
NP_MP_NS_SENSITIVITY_RTT_PDF = (
    f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/root_to_tip_regression.pdf"
)
NP_MP_NS_SENSITIVITY_RTT_FILTERED_PANEL = (
    f"{NP_MP_NS_SENSITIVITY_DATA_DIR}/panel_taxa.rtt_filtered.tsv"
)
NP_MP_NS_SENSITIVITY_RTT_FILTERED_DATES = (
    f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/dates.rtt_filtered.tsv"
)
NP_MP_NS_SENSITIVITY_RTT_EXCLUSIONS = (
    f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/rtt_outlier_exclusions.tsv"
)
NP_MP_NS_SENSITIVITY_RTT_SUMMARY = f"{NP_MP_NS_SENSITIVITY_RTT_DIR}/summary.tsv"

NP_MP_NS_RTT_SENSITIVITY_TARGETS = [
    NP_MP_NS_SENSITIVITY_PANEL_TAXA,
    NP_MP_NS_SENSITIVITY_PANEL_AUDIT,
    NP_MP_NS_SENSITIVITY_SUBSET_ALIGNMENT,
    NP_MP_NS_SENSITIVITY_SUBSET_TREE,
    NP_MP_NS_SENSITIVITY_SUBSET_AUDIT,
    NP_MP_NS_SENSITIVITY_DATES,
    NP_MP_NS_SENSITIVITY_RTT_OUTLIERS,
    NP_MP_NS_SENSITIVITY_RTT_CSV,
    NP_MP_NS_SENSITIVITY_RTT_PDF,
    NP_MP_NS_SENSITIVITY_RTT_LOG,
    NP_MP_NS_SENSITIVITY_RTT_DONE,
    NP_MP_NS_SENSITIVITY_RTT_FILTERED_PANEL,
    NP_MP_NS_SENSITIVITY_RTT_FILTERED_DATES,
    NP_MP_NS_SENSITIVITY_RTT_EXCLUSIONS,
    NP_MP_NS_SENSITIVITY_RTT_SUMMARY,
]

BEAST_CONFIG = config.get("beast", {})
BEAST_LOG_EVERY = int(BEAST_CONFIG.get("log_every", 2000))
BEAST_TREE_EVERY = int(BEAST_CONFIG.get("tree_every", BEAST_LOG_EVERY))
BEAST_ECHO_EVERY = int(BEAST_CONFIG.get("echo_every", BEAST_LOG_EVERY))

RESULTS_BEAST = config.get("results_beast", "results/beast")
BEAST_XML_DIR = f"{RESULTS_BEAST}/xml"
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


rule build_np_mp_ns_sensitivity_panel:
    input:
        tree=CONCAT_TREE,
    output:
        panel=NP_MP_NS_SENSITIVITY_PANEL_TAXA,
        audit=NP_MP_NS_SENSITIVITY_PANEL_AUDIT,
    params:
        forced_american_anchor=NP_MP_NS_FORCED_AMERICAN_ANCHOR,
        american_anchor_total=NP_MP_NS_AMERICAN_ANCHOR_TOTAL,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/build_np_mp_ns_sensitivity_panel.py \
            --tree {input.tree} \
            --panel-out {output.panel} \
            --audit-out {output.audit} \
            --forced-american-anchor-accession {params.forced_american_anchor} \
            --american-anchor-total {params.american_anchor_total}
        """


rule subset_np_mp_ns_sensitivity_alignment_and_prune_tree:
    input:
        alignment=CONCAT_ALIGNMENT,
        tree=CONCAT_TREE,
        taxa=NP_MP_NS_SENSITIVITY_PANEL_TAXA,
    output:
        alignment=NP_MP_NS_SENSITIVITY_SUBSET_ALIGNMENT,
        tree=NP_MP_NS_SENSITIVITY_SUBSET_TREE,
        audit=NP_MP_NS_SENSITIVITY_SUBSET_AUDIT,
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/02_Beast/subset_alignment_and_prune_tree.py \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment} \
            --out-tree {output.tree} \
            --audit {output.audit}
        """


rule build_np_mp_ns_sensitivity_treetime_dates:
    input:
        alignment=NP_MP_NS_SENSITIVITY_SUBSET_ALIGNMENT,
    output:
        dates=NP_MP_NS_SENSITIVITY_DATES,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.dates})
        python code/02_Beast/build_treetime_dates.py \
            --aln {input.alignment} \
            --out {output.dates}
        """


rule run_np_mp_ns_sensitivity_root_to_tip:
    input:
        tree=NP_MP_NS_SENSITIVITY_SUBSET_TREE,
        alignment=NP_MP_NS_SENSITIVITY_SUBSET_ALIGNMENT,
        dates=NP_MP_NS_SENSITIVITY_DATES,
    output:
        outliers=NP_MP_NS_SENSITIVITY_RTT_OUTLIERS,
        rtt_csv=NP_MP_NS_SENSITIVITY_RTT_CSV,
        rtt_pdf=NP_MP_NS_SENSITIVITY_RTT_PDF,
        log=NP_MP_NS_SENSITIVITY_RTT_LOG,
        done=NP_MP_NS_SENSITIVITY_RTT_DONE,
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


rule filter_np_mp_ns_sensitivity_panel_by_rtt_outliers:
    input:
        panel_tsv=NP_MP_NS_SENSITIVITY_PANEL_TAXA,
        rtt_outliers=NP_MP_NS_SENSITIVITY_RTT_OUTLIERS,
        dates_in=NP_MP_NS_SENSITIVITY_DATES,
    output:
        filtered_panel=NP_MP_NS_SENSITIVITY_RTT_FILTERED_PANEL,
        dates_out=NP_MP_NS_SENSITIVITY_RTT_FILTERED_DATES,
        exclusions=NP_MP_NS_SENSITIVITY_RTT_EXCLUSIONS,
        summary=NP_MP_NS_SENSITIVITY_RTT_SUMMARY,
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
            --exclusions-out {output.exclusions} \
            --summary-out {output.summary}
        """


rule np_mp_ns_rtt_sensitivity:
    input:
        NP_MP_NS_RTT_SENSITIVITY_TARGETS





rule subset_filtered_raw_segment_fasta:
    input:
        alignment="data/phylogeny/by_segment/H5N1_{segment}.fasta",
        taxa=BEAST_PANEL_RTT_FILTERED_TAXA,
    output:
        alignment=f"{BEAST_RAW_FINAL_SEGMENT_DIR}/H5N1_{{segment}}.fasta",
    wildcard_constraints:
        segment="|".join(PRE_BEAST_SEGMENTS),
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/02_Beast/subset_alignment_by_taxa.py \
            --alignment {input.alignment} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment}
        """





rule build_treetime_dates:
    input:
        alignment=BEAST_FILTERED_SUBSET_ALIGNMENT,
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


rule run_root_to_tip:
    input:
        tree=BEAST_FILTERED_SUBSET_TREE,
        alignment=BEAST_FILTERED_SUBSET_ALIGNMENT,
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


rule filter_beast_panel_by_rtt_outliers:
    input:
        panel_tsv=BEAST_PANEL_FILTERED_TAXA,
        rtt_outliers=PRE_BEAST_ROOT_TO_TIP_OUTLIERS,
        dates_in=PRE_BEAST_DATES,
    output:
        filtered_panel=BEAST_PANEL_RTT_FILTERED_TAXA,
        dates_out=BEAST_FINAL_DATES,
        exclusions=PRE_BEAST_RTT_EXCLUSIONS,
        summary=PRE_BEAST_RTT_EXCLUSIONS_SUMMARY,
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
            --exclusions-out {output.exclusions} \
            --summary-out {output.summary}
        """


rule publish_final_panel_concat_alignment:
    input:
        alignment=BEAST_FILTERED_SUBSET_ALIGNMENT,
        taxa=BEAST_PANEL_RTT_FILTERED_TAXA,
    output:
        alignment=BEAST_FINAL_SUBSET_ALIGNMENT,
        audit=BEAST_FINAL_SUBSET_AUDIT,
    conda:
        "../envs/01_ml_trees.yml"
    shell:
        r"""
        python code/02_Beast/subset_alignment_by_taxa.py \
            --alignment {input.alignment} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment} \
            --audit {output.audit}
        """


rule subset_filtered_segment_alignment:
    input:
        alignment="data/phylogeny/aligned/H5N1_{segment}.mafft",
        taxa=BEAST_PANEL_RTT_FILTERED_TAXA,
    output:
        alignment=f"{BEAST_FINAL_DATA_DIR}/H5N1_{{segment}}.fasta",
    wildcard_constraints:
        segment="|".join(PRE_BEAST_SEGMENTS),
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
        subset_alignment=BEAST_FINAL_SUBSET_ALIGNMENT,
    output:
        xml=f"{BEAST_XML_DIR}/{{scenario}}.xml",
    wildcard_constraints:
        scenario="|".join(BEAST_TEMPLATE_XMLS.keys()),
    params:
        output_prefix=lambda wildcards: f"results/beast/runs/{wildcards.scenario}/{wildcards.scenario}",
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


BEAST_PANEL_STATS_CSV = f"{RESULTS_QC_METRICS}/phylogeny/main_analysis_panel.csv"
BEAST_PANEL_STATS_SUMMARY = f"{RESULTS_QC_METRICS}/phylogeny/main_analysis_panel_summary.txt"


rule generate_main_analysis_panel_statistics:
    input:
        panel_tsv=BEAST_PANEL_RTT_FILTERED_TAXA,
        metadata_csv=config["flu_filtrado"],
    output:
        csv=BEAST_PANEL_STATS_CSV,
        summary=BEAST_PANEL_STATS_SUMMARY,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/generate_panel_statistics.py {input.panel_tsv} {output.csv} {input.metadata_csv}
        """
