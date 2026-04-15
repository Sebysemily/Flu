import json
import shlex
from pathlib import Path

CONCAT_TREE = "results/phylogeny/raxml/full_concat/H5N1_full_concat_beast.raxml.supportTBE"
CONCAT_ALIGNMENT = "data/phylogeny/aligned/H5N1_full_concat_beast.mafft"
CONCAT_PARTITIONS = "data/phylogeny/H5N1_full_concat_beast.partitions"
CONTEXT_META = config.get("context_metadata_tsv", "config/final_metadata_50_per_country_isolates.tsv")
IQTREE_MODELTEST_THREADS = int(config.get("iqtree_modeltest_threads", 8))
BETS_CONFIG = config.get("beast_pre", {}).get("bets", {})
BETS_THREADS = int(BETS_CONFIG.get("threads", 4))

BEAST_PRE_DATA_DIR = "data/beast_pre/panels"
BEAST_PRE_RESULTS_DIR = "results/beast_pre"
BEAST_PANEL_TAXA = "data/beast_pre/panels/panel_main_taxa.tsv"
BEAST_SUBSET_ALIGNMENT = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.subset.fasta"
BEAST_SUBSET_TREE = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.subset.nwk"
BEAST_SUBSET_AUDIT = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.subset.audit.tsv"
IQTREE_CONCAT_PARTITIONS = f"{BEAST_PRE_DATA_DIR}/panel_main_concat.subset.iqtree.nex"
BEAST_MODEL_TEST_PREFIX = f"{BEAST_PRE_RESULTS_DIR}/model_test/panel_main_concat"
PRE_BEAST_MODEL_TEST_SUMMARY = f"{BEAST_MODEL_TEST_PREFIX}.iqtree"
PRE_BEAST_MODEL_TEST_BEST_MODEL = f"{BEAST_MODEL_TEST_PREFIX}.best_model.txt"
PRE_BEAST_DATES = f"{BEAST_PRE_RESULTS_DIR}/root_to_tip/dates_from_headers.tsv"
PRE_BEAST_ROOT_TO_TIP_LOG = f"{BEAST_PRE_RESULTS_DIR}/root_to_tip/treetime_clock.log"
PRE_BEAST_ROOT_TO_TIP_DONE = f"{BEAST_PRE_RESULTS_DIR}/root_to_tip/treetime_clock.done"
PRE_BEAST_BETS_DIR = f"{BEAST_PRE_RESULTS_DIR}/bets"
PRE_BEAST_BETS_PACKAGE_DIR = BETS_CONFIG.get("package_dir", f"{PRE_BEAST_BETS_DIR}/packages")
PRE_BEAST_BETS_PACKAGE_READY = f"{PRE_BEAST_BETS_DIR}/model_selection_package.ready"
PRE_BEAST_BETS_RUN_SUMMARY = f"{PRE_BEAST_BETS_DIR}/bets_runs.tsv"
PRE_BEAST_BETS_BAYES_FACTORS = f"{PRE_BEAST_BETS_DIR}/bets_bayes_factors.tsv"

BETS_SCENARIO_XMLS = {
    "strict_constant_heterochronous": BETS_CONFIG.get(
        "strict_constant_heterochronous_xml",
        "config/beast_pre/bets/strict_constant_heterochronous.xml",
    ),
    "strict_constant_isochronous": BETS_CONFIG.get(
        "strict_constant_isochronous_xml",
        "config/beast_pre/bets/strict_constant_isochronous.xml",
    ),
    "ucln_constant_heterochronous": BETS_CONFIG.get(
        "ucln_constant_heterochronous_xml",
        "config/beast_pre/bets/ucln_constant_heterochronous.xml",
    ),
    "ucln_constant_isochronous": BETS_CONFIG.get(
        "ucln_constant_isochronous_xml",
        "config/beast_pre/bets/ucln_constant_isochronous.xml",
    ),
}

BETS_SCENARIOS = list(BETS_SCENARIO_XMLS.keys())
BETS_ENABLED_CONFIG = BETS_CONFIG.get("enabled", None)
BETS_XMLS_PRESENT = all(Path(xml_path).exists() for xml_path in BETS_SCENARIO_XMLS.values())
BETS_ENABLED = BETS_XMLS_PRESENT if BETS_ENABLED_CONFIG is None else bool(BETS_ENABLED_CONFIG)
BETS_SCENARIO_STDOUT = {
    scenario: f"{PRE_BEAST_BETS_DIR}/{scenario}/beast.stdout.log"
    for scenario in BETS_SCENARIOS
}
BETS_SCENARIO_STDERR = {
    scenario: f"{PRE_BEAST_BETS_DIR}/{scenario}/beast.stderr.log"
    for scenario in BETS_SCENARIOS
}
BETS_SCENARIO_DONE = {
    scenario: f"{PRE_BEAST_BETS_DIR}/{scenario}/run.done"
    for scenario in BETS_SCENARIOS
}


rule build_beast_panels:
    input:
        tree=CONCAT_TREE,
        context_metadata=CONTEXT_META,
    output:
        panel_main=BEAST_PANEL_TAXA,
        audit="data/beast_pre/panels/panel_selection_audit.tsv",
        country_month_audit="data/beast_pre/panels/panel_country_month_coverage.tsv",
    conda:
        "../envs/02_pre_beast.yml"
    params:
        max_cluster_dist=config.get("beast_pre", {}).get("max_cluster_dist", 0.08),
        n_per_country=config.get("beast_pre", {}).get("n_per_country", 4),
        n_total=config.get("beast_pre", {}).get("n_total", 60),
        min_mrca_support=config.get("beast_pre", {}).get("min_mrca_support", 70.0),
        max_per_country_month=config.get("beast_pre", {}).get("max_per_country_month", 2),
        usa_distal_quota=config.get("beast_pre", {}).get("usa_distal_quota", 0),
        additional_american_anchor_quota=config.get("beast_pre", {}).get(
            "additional_american_anchor_quota", 1
        ),
        forced_american_anchor_accession=config.get("beast_pre", {}).get(
            "forced_american_anchor_accession", "OQ968009"
        ),
        relaxed_min_mrca_support=config.get("beast_pre", {}).get(
            "relaxed_min_mrca_support", 50.0
        ),
        relaxed_fill=config.get("beast_pre", {}).get("relaxed_fill", 0),
        ecuador_core_json=shlex.quote(json.dumps(config.get("beast_pre", {}).get(
            "ecuador_core_ids",
            [
                "Flu-0316", "Flu-0317", "Flu-0580", "Flu-0582", "Flu-0583", "Flu-0584",
                "Flu-0586", "Flu-0589", "Flu-0592", "Flu-0593", "Flu-0596", "Flu-0599",
                "Flu-0600", "Flu-0604", "Flu-0608", "Flu-0610", "Flu-0611", "Flu-0613",
                "Flu-0614", "Flu-0619", "Flu-0621", "Flu-0622", "Flu-0623", "Flu-0630",
                "Flu-0641", "Flu-0652", "Flu-0653", "Flu-0654",
            ],
        ))),
        regional_blacklist_tokens_json=shlex.quote(json.dumps(config.get("beast_pre", {}).get(
            "regional_blacklist_tokens",
            ["__eurasian_anchor", "__american_anchor", "__usa_"],
        ))),
    shell:
        r"""
        mkdir -p data/beast_pre/panels
        python code/02_Beast/build_beast_panels.py \
            --tree {input.tree} \
            --context-metadata {input.context_metadata} \
            --panel-main-out {output.panel_main} \
            --audit-out {output.audit} \
            --country-month-audit-out {output.country_month_audit} \
            --max-cluster-dist {params.max_cluster_dist} \
            --n-per-country {params.n_per_country} \
            --n-total {params.n_total} \
            --min-mrca-support {params.min_mrca_support} \
            --max-per-country-month {params.max_per_country_month} \
            --usa-distal-quota {params.usa_distal_quota} \
            --additional-american-anchor-quota {params.additional_american_anchor_quota} \
            --forced-american-anchor-accession {params.forced_american_anchor_accession} \
            --relaxed-min-mrca-support {params.relaxed_min_mrca_support} \
            --relaxed-fill {params.relaxed_fill} \
            --ecuador-core-json {params.ecuador_core_json} \
            --regional-blacklist-tokens-json {params.regional_blacklist_tokens_json}
        """


rule subset_panel_concat_alignment_and_prune_tree:
    input:
        alignment=CONCAT_ALIGNMENT,
        tree=CONCAT_TREE,
        taxa=BEAST_PANEL_TAXA,
    output:
        alignment=BEAST_SUBSET_ALIGNMENT,
        tree=BEAST_SUBSET_TREE,
        audit=BEAST_SUBSET_AUDIT,
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


rule build_concat_iqtree_partition_nexus:
    input:
        partitions=CONCAT_PARTITIONS,
    output:
        partitions=IQTREE_CONCAT_PARTITIONS,
    shell:
        r"""
        python - <<'PY'
from pathlib import Path

src = Path("{input.partitions}")
dst = Path("{output.partitions}")
dst.parent.mkdir(parents=True, exist_ok=True)

lines = []
with src.open("r", encoding="utf-8") as handle:
    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        _, rhs = line.split(",", 1)
        name, ranges = rhs.split("=", 1)
        name = name.strip()
        ranges = " ".join(part.strip() for part in ranges.split(","))
        lines.append("    charset " + name + " = DNA, " + ranges + ";")

with dst.open("w", encoding="utf-8") as handle:
    handle.write("#nexus\n")
    handle.write("begin sets;\n")
    for line in lines:
        handle.write(line + "\n")
    handle.write("end;\n")
PY
        """


rule build_treetime_dates:
    input:
        alignment=BEAST_SUBSET_ALIGNMENT,
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
        tree=BEAST_SUBSET_TREE,
        alignment=BEAST_SUBSET_ALIGNMENT,
        dates=PRE_BEAST_DATES,
    output:
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


rule ensure_beast_model_selection_package:
    output:
        ready=PRE_BEAST_BETS_PACKAGE_READY,
    params:
        package_dir=PRE_BEAST_BETS_PACKAGE_DIR,
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        mkdir -p {params.package_dir}
        export JAVA_TOOL_OPTIONS="${{JAVA_TOOL_OPTIONS:-}} -Dbeast.user.package.dir=$(realpath {params.package_dir})"
        packagemanager -add MODEL_SELECTION
        touch {output.ready}
        """


rule run_bets_scenario:
    wildcard_constraints:
        scenario="|".join(BETS_SCENARIOS),
    input:
        xml=lambda wildcards: BETS_SCENARIO_XMLS[wildcards.scenario],
        package_ready=PRE_BEAST_BETS_PACKAGE_READY,
        root_to_tip_done=PRE_BEAST_ROOT_TO_TIP_DONE,
    output:
        stdout=f"{PRE_BEAST_BETS_DIR}/{{scenario}}/beast.stdout.log",
        stderr=f"{PRE_BEAST_BETS_DIR}/{{scenario}}/beast.stderr.log",
        done=f"{PRE_BEAST_BETS_DIR}/{{scenario}}/run.done",
    params:
        package_dir=PRE_BEAST_BETS_PACKAGE_DIR,
    threads: BETS_THREADS
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/run_beast_bets.py \
            --xml {input.xml} \
            --package-dir {params.package_dir} \
            --threads {threads} \
            --stdout {output.stdout} \
            --stderr {output.stderr} \
            --done {output.done}
        """


rule summarize_bets:
    input:
        stdout=list(BETS_SCENARIO_STDOUT.values()),
        stderr=list(BETS_SCENARIO_STDERR.values()),
        done=list(BETS_SCENARIO_DONE.values()),
    output:
        summary=PRE_BEAST_BETS_RUN_SUMMARY,
        bayes_factors=PRE_BEAST_BETS_BAYES_FACTORS,
    params:
        strict_hetero_stdout=BETS_SCENARIO_STDOUT["strict_constant_heterochronous"],
        strict_iso_stdout=BETS_SCENARIO_STDOUT["strict_constant_isochronous"],
        ucln_hetero_stdout=BETS_SCENARIO_STDOUT["ucln_constant_heterochronous"],
        ucln_iso_stdout=BETS_SCENARIO_STDOUT["ucln_constant_isochronous"],
        strict_hetero_stderr=BETS_SCENARIO_STDERR["strict_constant_heterochronous"],
        strict_iso_stderr=BETS_SCENARIO_STDERR["strict_constant_isochronous"],
        ucln_hetero_stderr=BETS_SCENARIO_STDERR["ucln_constant_heterochronous"],
        ucln_iso_stderr=BETS_SCENARIO_STDERR["ucln_constant_isochronous"],
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        python code/02_Beast/summarize_bets.py \
            --out-summary {output.summary} \
            --out-bayes-factors {output.bayes_factors} \
            --scenario-log strict_constant_heterochronous={params.strict_hetero_stdout} \
            --scenario-log strict_constant_isochronous={params.strict_iso_stdout} \
            --scenario-log ucln_constant_heterochronous={params.ucln_hetero_stdout} \
            --scenario-log ucln_constant_isochronous={params.ucln_iso_stdout} \
            --scenario-err strict_constant_heterochronous={params.strict_hetero_stderr} \
            --scenario-err strict_constant_isochronous={params.strict_iso_stderr} \
            --scenario-err ucln_constant_heterochronous={params.ucln_hetero_stderr} \
            --scenario-err ucln_constant_isochronous={params.ucln_iso_stderr}
        """


rule iqtree_model_test_concat_subset:
    input:
        alignment=BEAST_SUBSET_ALIGNMENT,
        partitions=IQTREE_CONCAT_PARTITIONS,
        root_to_tip_done=PRE_BEAST_ROOT_TO_TIP_DONE,
        bets_summary=(PRE_BEAST_BETS_RUN_SUMMARY if BETS_ENABLED else []),
    output:
        summary=PRE_BEAST_MODEL_TEST_SUMMARY,
        best_model=PRE_BEAST_MODEL_TEST_BEST_MODEL,
    params:
        prefix=BEAST_MODEL_TEST_PREFIX,
    threads: IQTREE_MODELTEST_THREADS
    conda:
        "../envs/02_pre_beast.yml"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        iqtree2 \
            -s {input.alignment} \
            -p {input.partitions} \
            -m MFP \
            -mset HKY,GTR \
            -redo \
            -pre {params.prefix} \
            -nt {threads}
        if [ -f "{params.prefix}.iqtree" ]; then
            grep -m1 'Best-fit model' "{params.prefix}.iqtree" | sed -E 's/^[^:]+:[[:space:]]*//;s/[[:space:]]*$//' > "{params.prefix}.best_model.txt" || true
        else
            echo "" > "{params.prefix}.best_model.txt"
        fi
        """
