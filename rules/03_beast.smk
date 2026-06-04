BEAST_ENABLED = bool(config.get("beast", {}).get("enabled", True))
# BEAST 1.x uses ~1 CPU per chain; Snakemake threads reserve a core per job.
BEAST_THREADS_PER_CHAIN = int(config.get("beast", {}).get("threads_per_chain", 1))
BEAST_MAX_PARALLEL_CHAINS = int(config.get("beast", {}).get("max_parallel_chains", 2))
BEAST_THREADS = int(config.get("max_threads", 16))
BEAST_MAX_HOURS = int(config.get("beast", {}).get("max_hours", 12))
BEAST_BINARY = config.get("beast", {}).get("binary") or ""
BEAST_SEED_OFFSET = 100000
BEAST_REPLICATES = ["r1", "r2"]
BEAST_PREPARED_XMLS = {
    "strict_const": "template_beast/strict_const.xml",
    "strict_exp": "template_beast/strict_exp.xml",
    "strict_bay": "template_beast/strict_bay.xml",
    "ucln_const": "template_beast/ucln_const.xml",
    "ucln_exp": "template_beast/ucln_exp.xml",
    "ucln_bay": "template_beast/ucln_bay.xml"
}
BEAST_RUN_SCENARIOS = list(BEAST_PREPARED_XMLS.keys())
BEAST_SEEDS_CONFIG = config.get("beast", {}).get("seeds", {})
BEAST_BEAGLE_CONFIG = config.get("beast", {}).get("beagle", {})
BEAST_CHAIN_LENGTHS = config.get("beast", {}).get("chain_length", {})
BEAST_DEFAULT_SEEDS = {
    "strict_constant": 1001,
    "strict_constant_lugar": 1005,
    "ucln_constant": 1002,
    "strict_exp": 1003,
    "ucln_exp": 1004,
}
RESULTS_BEAST = "results/phylogeny/bayesian"
STRICT_CONSTANT_FINAL_DIR = f"{RESULTS_BEAST}/final/time"
STRICT_CONSTANT_FINAL_PREFIX = f"{STRICT_CONSTANT_FINAL_DIR}/strict_constant"
STRICT_CONSTANT_FINAL_COMBINED_LOG = f"{STRICT_CONSTANT_FINAL_PREFIX}.combined.log"
STRICT_CONSTANT_FINAL_COMBINED_TREES = f"{STRICT_CONSTANT_FINAL_PREFIX}.combined.trees"
STRICT_CONSTANT_FINAL_TREE = f"{STRICT_CONSTANT_FINAL_PREFIX}.mcc.mean.tree"
STRICT_CONSTANT_FINAL_DONE = f"{STRICT_CONSTANT_FINAL_DIR}/run.done"
STRICT_CONSTANT_CHAIN_LENGTH = int(BEAST_CHAIN_LENGTHS.get("strict_constant", 100000000))
STRICT_CONSTANT_BURNIN_STATES = int(STRICT_CONSTANT_CHAIN_LENGTH * 0.10)
STRICT_CONSTANT_LUGAR_FINAL_DIR = f"{RESULTS_BEAST}/final/geography"
STRICT_CONSTANT_LUGAR_FINAL_PREFIX = f"{STRICT_CONSTANT_LUGAR_FINAL_DIR}/strict_constant_lugar"
STRICT_CONSTANT_LUGAR_FINAL_COMBINED_LOG = f"{STRICT_CONSTANT_LUGAR_FINAL_PREFIX}.combined.log"
STRICT_CONSTANT_LUGAR_FINAL_COMBINED_TREES = f"{STRICT_CONSTANT_LUGAR_FINAL_PREFIX}.combined.trees"
STRICT_CONSTANT_LUGAR_FINAL_TREE = f"{STRICT_CONSTANT_LUGAR_FINAL_PREFIX}.mcc.mean.tree"
STRICT_CONSTANT_LUGAR_FINAL_DONE = f"{STRICT_CONSTANT_LUGAR_FINAL_DIR}/run.done"
STRICT_CONSTANT_LUGAR_CHAIN_LENGTH = int(BEAST_CHAIN_LENGTHS.get("strict_constant_lugar", 100000000))
STRICT_CONSTANT_LUGAR_BURNIN_STATES = int(STRICT_CONSTANT_LUGAR_CHAIN_LENGTH * 0.10)


def beast_seed_for(scenario, replicate):
    explicit_key = f"{scenario}_rep{1 if replicate == 'r1' else 2}"
    if explicit_key in BEAST_SEEDS_CONFIG:
        return int(BEAST_SEEDS_CONFIG[explicit_key])

    base_default = BEAST_DEFAULT_SEEDS.get(scenario, 1001)
    base_seed = int(BEAST_SEEDS_CONFIG.get(scenario, base_default))
    return base_seed + (BEAST_SEED_OFFSET if replicate == "r2" else 0)


BEAST_RUN_TARGETS = expand(f"{RESULTS_BEAST}/runs/{{scenario}}/run.done", scenario=BEAST_RUN_SCENARIOS)
BEAST_FINAL_TARGETS = [STRICT_CONSTANT_FINAL_DONE, STRICT_CONSTANT_LUGAR_FINAL_DONE]

# Only target scenarios for GSS test runs that actually have an XML file created
BEAST_GSS_SCENARIOS = [s for s, xml in BEAST_PREPARED_XMLS.items() if os.path.exists(xml)]
BEAST_GSS_TARGETS = expand("results/beast/GSS/{scenario}/run.done", scenario=BEAST_GSS_SCENARIOS)

BEAST_PUBLIC_TARGETS = (
    [*BEAST_RUN_TARGETS, *BEAST_FINAL_TARGETS, *BEAST_GSS_TARGETS]
    if BEAST_ENABLED
    else []
)


# =====================================================================
# Rule: update_beast_xml_params
# Updates the XML in-place if config.yml changes.
# =====================================================================
rule update_beast_xml_params:
    input:
        config="config/config.yml"
    output:
        stamp="results/beast/GSS/{scenario}/xml_params.json",
    params:
        xml=lambda wildcards: BEAST_PREPARED_XMLS[wildcards.scenario]
    shell:
        r"""
        python code/03_beast/update_beast_xmls.py {wildcards.scenario} \
            --config {input.config} \
            --stamp {output.stamp}
        """

# =====================================================================
# Rule: run_beast_single
# Single chain runs with -resume support, reading directly from template_beast.
# =====================================================================
rule run_beast_single:
    input:
        xml=lambda wildcards: BEAST_PREPARED_XMLS[wildcards.scenario],
        stamp="results/beast/GSS/{scenario}/xml_params.json",
    output:
        done="results/beast/GSS/{scenario}/run.done"
    params:
        outdir="results/beast/GSS/{scenario}",
        seed=config.get("random_seed", 1001)
    threads: BEAST_THREADS_PER_CHAIN
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        # Calculate relative path to xml from the working directory
        REL_XML="../../../../{input.xml}"
        
        # Find checkpoint file if it exists
        CHKPT=$(ls *.state *.chkpt 2>/dev/null | head -n 1 || true)
        
        if [ -n "$CHKPT" ]; then
            echo "Resuming from checkpoint $CHKPT..."
            # BEAST 1.x does not append natively without truncating. Protect existing data by renaming it.
            TIMESTAMP=$(date +%s)
            mkdir -p logs
            for f in *.log *.trees *.ops; do
                if [ -f "$f" ]; then
                    mv "$f" "logs/$f.part_$TIMESTAMP"
                fi
            done
            # Also move any existing part files to keep it tidy
            for f in *.part_*; do
                if [ -f "$f" ]; then
                    mv "$f" "logs/"
                fi
            done

            # Calculate remaining states needed to reach the target chainLength
            TARGET_LENGTH=$(grep -oP 'chainLength="\K[0-9]+' $REL_XML | head -n 1)
            LAST_STATE=$(cat logs/*.log.part_* 2>/dev/null | awk '/^[0-9]/ {{print $1}}' | sort -n | tail -n 1 || echo 0)
            
            if [ -z "$LAST_STATE" ]; then LAST_STATE=0; fi

            REMAINING=$(( TARGET_LENGTH - LAST_STATE ))
            echo "Target length: $TARGET_LENGTH, Last state: $LAST_STATE, Remaining: $REMAINING"

            if [ "$REMAINING" -le 0 ]; then
                echo "Target chainLength already reached or exceeded ($LAST_STATE >= $TARGET_LENGTH). Finishing gracefully."
            else
                # Modify a local copy of the XML to only run the remaining states
                cp $REL_XML local.xml
                sed -i "s/chainLength=\"[0-9]*\"/chainLength=\"$REMAINING\"/" local.xml
                beast -load_state "$CHKPT" -seed {params.seed} local.xml
            fi
        else
            if ls *.log *.trees >/dev/null 2>&1; then
                echo "ERROR: Logs exist but no checkpoint was found! Aborting to protect your data."
                exit 1
            fi
            echo "No checkpoint found. Starting fresh BEAST chain..."
            beast -seed {params.seed} $REL_XML
        fi
        
        touch run.done
        """


rule run_beast_replicate:
    input:
        xml=lambda wildcards: BEAST_PREPARED_XMLS[wildcards.scenario],
        previous_done=lambda wildcards: []
        if wildcards.replicate == "r1"
        else f"{RESULTS_BEAST}/runs/{wildcards.scenario}/r1/run.done",
    output:
        status=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/status.log",
        stdout=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/stdout.log",
        stderr=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/stderr.log",
        beast_log=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/{{scenario}}.log",
        beast_trees=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/{{scenario}}.trees",
        beast_chkpt=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/{{scenario}}.chkpt",
        beast_ops=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/{{scenario}}.ops",
        done=f"{RESULTS_BEAST}/runs/{{scenario}}/{{replicate}}/run.done",
    params:
        beast_binary=BEAST_BINARY,
        output_prefix=lambda wildcards: f"{RESULTS_BEAST}/runs/{wildcards.scenario}/{wildcards.scenario}",
        replicate_dir=lambda wildcards: f"{RESULTS_BEAST}/runs/{wildcards.scenario}/{wildcards.replicate}",
        seed=lambda wildcards: beast_seed_for(wildcards.scenario, wildcards.replicate),
        beagle_mode=lambda wildcards: str(
            BEAST_BEAGLE_CONFIG.get(
                "mode",
                "auto" if bool(BEAST_BEAGLE_CONFIG.get("enabled", False)) else "off",
            )
        ),
        beagle_resource=lambda wildcards: str(BEAST_BEAGLE_CONFIG.get("resource", "auto")),
        beagle_vendor=lambda wildcards: str(BEAST_BEAGLE_CONFIG.get("vendor", "any")),
        beagle_platform=lambda wildcards: str(BEAST_BEAGLE_CONFIG.get("platform", "auto")),
        beagle_precision=lambda wildcards: str(BEAST_BEAGLE_CONFIG.get("precision", "auto")),
        beagle_scaling=lambda wildcards: str(BEAST_BEAGLE_CONFIG.get("scaling", "default")),
        beagle_info=lambda wildcards: str(bool(BEAST_BEAGLE_CONFIG.get("info", False))).lower(),
        beagle_threads=lambda wildcards: str(BEAST_BEAGLE_CONFIG.get("threads", "")),
        beagle_fallback_to_cpu=lambda wildcards: str(
            bool(BEAST_BEAGLE_CONFIG.get("fallback_to_cpu", False))
        ).lower(),
    wildcard_constraints:
        scenario="|".join(BEAST_RUN_SCENARIOS),
        replicate="|".join(BEAST_REPLICATES),
    threads: BEAST_THREADS
    resources:
        runtime=BEAST_MAX_HOURS * 60
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        python code/03_beast/run_beast_replicate.py \
            --xml {input.xml} \
            --beast-binary '{params.beast_binary}' \
            --scenario {wildcards.scenario} \
            --replicate {wildcards.replicate} \
            --output-prefix {params.output_prefix} \
            --replicate-dir {params.replicate_dir} \
            --seed {params.seed} \
            --threads {threads} \
            --beagle-mode {params.beagle_mode} \
            --beagle-resource {params.beagle_resource} \
            --beagle-vendor {params.beagle_vendor} \
            --beagle-platform {params.beagle_platform} \
            --beagle-precision {params.beagle_precision} \
            --beagle-scaling {params.beagle_scaling} \
            --beagle-info {params.beagle_info} \
            --beagle-threads '{params.beagle_threads}' \
            --beagle-fallback-to-cpu {params.beagle_fallback_to_cpu} \
            --status {output.status} \
            --stdout {output.stdout} \
            --stderr {output.stderr} \
            --done {output.done}
        """


rule summarize_beast_run:
    input:
        run_r1=f"{RESULTS_BEAST}/runs/{{scenario}}/r1/run.done",
        run_r2=f"{RESULTS_BEAST}/runs/{{scenario}}/r2/run.done",
    output:
        done=f"{RESULTS_BEAST}/runs/{{scenario}}/run.done",
    params:
        seed_r1=lambda wildcards: beast_seed_for(wildcards.scenario, "r1"),
        seed_r2=lambda wildcards: beast_seed_for(wildcards.scenario, "r2"),
    wildcard_constraints:
        scenario="|".join(BEAST_RUN_SCENARIOS),
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        python code/02_Beast/summarize_beast_run.py \
            --scenario {wildcards.scenario} \
            --run-done {input.run_r1} \
            --run-done {input.run_r2} \
            --seed {params.seed_r1} \
            --seed {params.seed_r2} \
            --out {output.done}
        """


rule combine_strict_constant_logs:
    input:
        scenario_done=f"{RESULTS_BEAST}/runs/strict_constant/run.done",
        log_r1=f"{RESULTS_BEAST}/runs/strict_constant/r1/strict_constant.log",
        log_r2=f"{RESULTS_BEAST}/runs/strict_constant/r2/strict_constant.log",
    output:
        combined=STRICT_CONSTANT_FINAL_COMBINED_LOG,
    params:
        burnin=STRICT_CONSTANT_BURNIN_STATES,
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_FINAL_DIR}
        rm -f {output.combined}
        logcombiner \
            -burnin {params.burnin} \
            {input.log_r1} \
            {input.log_r2} \
            {output.combined}
        """


rule combine_strict_constant_trees:
    input:
        scenario_done=f"{RESULTS_BEAST}/runs/strict_constant/run.done",
        trees_r1=f"{RESULTS_BEAST}/runs/strict_constant/r1/strict_constant.trees",
        trees_r2=f"{RESULTS_BEAST}/runs/strict_constant/r2/strict_constant.trees",
    output:
        combined=STRICT_CONSTANT_FINAL_COMBINED_TREES,
    params:
        burnin=STRICT_CONSTANT_BURNIN_STATES,
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_FINAL_DIR}
        rm -f {output.combined}
        logcombiner \
            -trees \
            -burnin {params.burnin} \
            {input.trees_r1} \
            {input.trees_r2} \
            {output.combined}
        """


rule annotate_strict_constant_final_tree:
    input:
        trees=STRICT_CONSTANT_FINAL_COMBINED_TREES,
    output:
        tree=STRICT_CONSTANT_FINAL_TREE,
    threads: BEAST_THREADS
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_FINAL_DIR}
        rm -f {output.tree}
        treeannotator \
            -type mcc \
            -heights mean \
            -threads {threads} \
            {input.trees} \
            {output.tree}
        """


rule summarize_strict_constant_final_time:
    input:
        combined_log=STRICT_CONSTANT_FINAL_COMBINED_LOG,
        combined_trees=STRICT_CONSTANT_FINAL_COMBINED_TREES,
        tree=STRICT_CONSTANT_FINAL_TREE,
    output:
        done=STRICT_CONSTANT_FINAL_DONE,
    params:
        burnin=STRICT_CONSTANT_BURNIN_STATES,
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_FINAL_DIR}
        cat > {output.done} <<'EOF'
combined_log	{input.combined_log}
combined_trees	{input.combined_trees}
final_tree	{input.tree}
burnin_states_per_replicate	{params.burnin}
tree_summary	MCC_mean_heights
EOF
        """


rule combine_strict_constant_lugar_logs:
    input:
        scenario_done=f"{RESULTS_BEAST}/runs/strict_constant_lugar/run.done",
        log_r1=f"{RESULTS_BEAST}/runs/strict_constant_lugar/r1/strict_constant_lugar.log",
        log_r2=f"{RESULTS_BEAST}/runs/strict_constant_lugar/r2/strict_constant_lugar.log",
    output:
        combined=STRICT_CONSTANT_LUGAR_FINAL_COMBINED_LOG,
    params:
        burnin=STRICT_CONSTANT_LUGAR_BURNIN_STATES,
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_LUGAR_FINAL_DIR}
        rm -f {output.combined}
        logcombiner \
            -burnin {params.burnin} \
            {input.log_r1} \
            {input.log_r2} \
            {output.combined}
        """


rule combine_strict_constant_lugar_trees:
    input:
        scenario_done=f"{RESULTS_BEAST}/runs/strict_constant_lugar/run.done",
        trees_r1=f"{RESULTS_BEAST}/runs/strict_constant_lugar/r1/strict_constant_lugar.trees",
        trees_r2=f"{RESULTS_BEAST}/runs/strict_constant_lugar/r2/strict_constant_lugar.trees",
    output:
        combined=STRICT_CONSTANT_LUGAR_FINAL_COMBINED_TREES,
    params:
        burnin=STRICT_CONSTANT_LUGAR_BURNIN_STATES,
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_LUGAR_FINAL_DIR}
        rm -f {output.combined}
        logcombiner \
            -trees \
            -burnin {params.burnin} \
            {input.trees_r1} \
            {input.trees_r2} \
            {output.combined}
        """


rule annotate_strict_constant_lugar_final_tree:
    input:
        trees=STRICT_CONSTANT_LUGAR_FINAL_COMBINED_TREES,
    output:
        tree=STRICT_CONSTANT_LUGAR_FINAL_TREE,
    threads: BEAST_THREADS
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_LUGAR_FINAL_DIR}
        rm -f {output.tree}
        treeannotator \
            -type mcc \
            -heights mean \
            -threads {threads} \
            {input.trees} \
            {output.tree}
        """


rule summarize_strict_constant_lugar_final_geography:
    input:
        combined_log=STRICT_CONSTANT_LUGAR_FINAL_COMBINED_LOG,
        combined_trees=STRICT_CONSTANT_LUGAR_FINAL_COMBINED_TREES,
        tree=STRICT_CONSTANT_LUGAR_FINAL_TREE,
    output:
        done=STRICT_CONSTANT_LUGAR_FINAL_DONE,
    params:
        burnin=STRICT_CONSTANT_LUGAR_BURNIN_STATES,
    shell:
        r"""
        mkdir -p {STRICT_CONSTANT_LUGAR_FINAL_DIR}
        cat > {output.done} <<'EOF'
combined_log	{input.combined_log}
combined_trees	{input.combined_trees}
final_tree	{input.tree}
burnin_states_per_replicate	{params.burnin}
tree_summary	MCC_mean_heights
EOF
        """
