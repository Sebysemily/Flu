CONCAT_TREE = "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE"
CONCAT_ALIGNMENT = "data/phylogeny/aligned/H5N1_full_concat_beast.mafft"
CONTEXT_META = config.get("context_metadata_tsv", "config/final_metadata_50_per_country_isolates.tsv")

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

BEAST_CONFIG = config.get("beast", {})
BEAST_LOG_EVERY = int(BEAST_CONFIG.get("log_every", 2000))
BEAST_TREE_EVERY = int(BEAST_CONFIG.get("tree_every", BEAST_LOG_EVERY))
BEAST_ECHO_EVERY = int(BEAST_CONFIG.get("echo_every", BEAST_LOG_EVERY))

BEAST_XML_DIR = "results/beast/xml"
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


rule build_beast_panels:
    input:
        tree=CONCAT_TREE,
        context_metadata=CONTEXT_META,
    output:
        panel_main=BEAST_PANEL_TAXA,
        audit=f"{BEAST_PRE_DATA_DIR}/panel_selection_audit.tsv",
        country_month_audit=f"{BEAST_PRE_DATA_DIR}/panel_country_month_coverage.tsv",
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        mkdir -p {BEAST_PRE_DATA_DIR}
        python code/02_Beast/build_beast_panels.py \
            --tree {input.tree} \
            --context-metadata {input.context_metadata} \
            --panel-main-out {output.panel_main} \
            --audit-out {output.audit} \
            --country-month-audit-out {output.country_month_audit}
        """


rule observe_beast_subset_source_qc:
    input:
        final_fasta="data/final/H5N1_final.fasta",
        panel_tsv=BEAST_PANEL_TAXA,
        ecuador_summary="data/input/H5N1_EC_summary.csv",
        context_summary="data/input/H5N1_context_summary.csv",
        ecuador_audit="data/assembled/ecuador_intermediate_audit.csv",
    output:
        metrics=PRE_BEAST_SOURCE_QC_METRICS,
        outliers=PRE_BEAST_SOURCE_QC_OUTLIERS,
        summary=PRE_BEAST_SOURCE_QC_SUMMARY,
        report=PRE_BEAST_SOURCE_QC_REPORT,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/observe_beast_subset_source_qc.py \
            --final-fasta {input.final_fasta} \
            --panel-taxa {input.panel_tsv} \
            --ecuador-summary {input.ecuador_summary} \
            --context-summary {input.context_summary} \
            --ecuador-audit {input.ecuador_audit} \
            --out-metrics {output.metrics} \
            --out-outliers {output.outliers} \
            --out-summary {output.summary} \
            --out-report {output.report}
        """


rule subset_filtered_raw_segment_fasta:
    input:
        alignment="data/phylogeny/by_segment/H5N1_{segment}.fasta",
        taxa=BEAST_PANEL_RTT_FILTERED_TAXA,
    output:
        alignment=f"{BEAST_RAW_FINAL_SEGMENT_DIR}/H5N1_{{segment}}.fasta",
    wildcard_constraints:
        segment="|".join(PRE_BEAST_SEGMENTS),
    conda:
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        python code/02_Beast/subset_alignment_by_taxa.py \
            --alignment {input.alignment} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment}
        """


rule filter_beast_panel_by_qc:
    input:
        panel_tsv=BEAST_PANEL_TAXA,
        source_qc_metrics=PRE_BEAST_SOURCE_QC_METRICS,
    output:
        filtered_panel=BEAST_PANEL_FILTERED_TAXA,
        exclusions=PRE_BEAST_PANEL_EXCLUSIONS,
        summary=PRE_BEAST_PANEL_EXCLUSIONS_SUMMARY,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/filter_beast_panel_by_qc.py \
            --panel-taxa {input.panel_tsv} \
            --source-qc-metrics {input.source_qc_metrics} \
            --filtered-panel-out {output.filtered_panel} \
            --exclusions-out {output.exclusions} \
            --summary-out {output.summary}
        """


rule subset_filtered_panel_concat_alignment_and_prune_tree:
    input:
        alignment=CONCAT_ALIGNMENT,
        tree=CONCAT_TREE,
        taxa=BEAST_PANEL_FILTERED_TAXA,
    output:
        alignment=BEAST_FILTERED_SUBSET_ALIGNMENT,
        tree=BEAST_FILTERED_SUBSET_TREE,
        audit=BEAST_FILTERED_SUBSET_AUDIT,
    conda:
        "../envs/ml_per_segment.yml"
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
        "../envs/ml_per_segment.yml"
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
        "../envs/ml_per_segment.yml"
    shell:
        r"""
        python code/02_Beast/subset_alignment_by_taxa.py \
            --alignment {input.alignment} \
            --taxa {input.taxa} \
            --out-alignment {output.alignment}
        """


rule observe_segment_alignment_qc:
    input:
        before_alignment=f"{BEAST_RAW_FINAL_SEGMENT_DIR}/H5N1_{{segment}}.fasta",
        after_alignment=f"{BEAST_FINAL_DATA_DIR}/H5N1_{{segment}}.fasta",
        taxa_tsv=BEAST_PANEL_RTT_FILTERED_TAXA,
    output:
        metrics=f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.metrics.tsv",
        outliers=f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.outliers.tsv",
        summary=f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.summary.tsv",
        report=f"{BEAST_SEGMENT_QC_DIR}/H5N1_{{segment}}.report.md",
    wildcard_constraints:
        segment="|".join(PRE_BEAST_SEGMENTS),
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/observe_subset_alignment_qc.py \
            --before-alignment {input.before_alignment} \
            --after-alignment {input.after_alignment} \
            --taxa-tsv {input.taxa_tsv} \
            --out-metrics {output.metrics} \
            --out-outliers {output.outliers} \
            --out-summary {output.summary} \
            --out-report {output.report}
        """


rule summarize_final_segment_qc:
    input:
        segment_metrics=PRE_BEAST_SEGMENT_QC_METRICS,
    output:
        metrics=PRE_BEAST_FINAL_SEGMENT_QC_METRICS,
        outliers=PRE_BEAST_FINAL_SEGMENT_QC_OUTLIERS,
        summary=PRE_BEAST_FINAL_SEGMENT_QC_SUMMARY,
        report=PRE_BEAST_FINAL_SEGMENT_QC_REPORT,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/summarize_final_segment_qc.py \
            --segment-metrics {input.segment_metrics} \
            --out-metrics {output.metrics} \
            --out-outliers {output.outliers} \
            --out-summary {output.summary} \
            --out-report {output.report}
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
