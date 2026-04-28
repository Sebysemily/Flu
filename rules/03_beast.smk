BEAST_ENABLED = bool(config.get("beast", {}).get("enabled", True))
BEAST_THREADS = int(config.get("beast", {}).get("threads", 4))
BEAST_MAX_HOURS = int(config.get("beast", {}).get("max_hours", 12))
BEAST_BINARY = config.get("beast", {}).get("binary") or ""
BEAST_SEED_OFFSET = 100000
BEAST_REPLICATES = ["r1", "r2"]
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
STRICT_CONSTANT_FINAL_DIR = "results/beast/final/time"
STRICT_CONSTANT_FINAL_PREFIX = f"{STRICT_CONSTANT_FINAL_DIR}/strict_constant"
STRICT_CONSTANT_FINAL_COMBINED_LOG = f"{STRICT_CONSTANT_FINAL_PREFIX}.combined.log"
STRICT_CONSTANT_FINAL_COMBINED_TREES = f"{STRICT_CONSTANT_FINAL_PREFIX}.combined.trees"
STRICT_CONSTANT_FINAL_TREE = f"{STRICT_CONSTANT_FINAL_PREFIX}.mcc.mean.tree"
STRICT_CONSTANT_FINAL_DONE = f"{STRICT_CONSTANT_FINAL_DIR}/run.done"
STRICT_CONSTANT_CHAIN_LENGTH = int(BEAST_CHAIN_LENGTHS.get("strict_constant", 100000000))
STRICT_CONSTANT_BURNIN_STATES = int(STRICT_CONSTANT_CHAIN_LENGTH * 0.10)
STRICT_CONSTANT_LUGAR_FINAL_DIR = "results/beast/final/geography"
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


BEAST_RUN_TARGETS = expand("results/beast/runs/{scenario}/run.done", scenario=BEAST_RUN_SCENARIOS)
BEAST_FINAL_TARGETS = [STRICT_CONSTANT_FINAL_DONE, STRICT_CONSTANT_LUGAR_FINAL_DONE]
BEAST_PUBLIC_TARGETS = (
    [*BEAST_RUN_TARGETS, *BEAST_FINAL_TARGETS]
    if BEAST_ENABLED
    else []
)


rule validate_beast_xml:
    input:
        xml=lambda wildcards: BEAST_PREPARED_XMLS[wildcards.scenario],
    output:
        validated="results/beast/xml/{scenario}.validated",
    wildcard_constraints:
        scenario="|".join(BEAST_RUN_SCENARIOS),
    conda:
        "../envs/03_beast.yml"
    shell:
        r"""
        python code/02_Beast/validate_beast_xml.py \
            --xml {input.xml} \
            --out {output.validated}
        """


rule run_beast_replicate:
    input:
        xml=lambda wildcards: BEAST_PREPARED_XMLS[wildcards.scenario],
        validated="results/beast/xml/{scenario}.validated",
        previous_done=lambda wildcards: []
        if wildcards.replicate == "r1"
        else f"results/beast/runs/{wildcards.scenario}/r1/run.done",
    output:
        status="results/beast/runs/{scenario}/{replicate}/status.log",
        stdout="results/beast/runs/{scenario}/{replicate}/stdout.log",
        stderr="results/beast/runs/{scenario}/{replicate}/stderr.log",
        beast_log="results/beast/runs/{scenario}/{replicate}/{scenario}.log",
        beast_trees="results/beast/runs/{scenario}/{replicate}/{scenario}.trees",
        beast_chkpt="results/beast/runs/{scenario}/{replicate}/{scenario}.chkpt",
        beast_ops="results/beast/runs/{scenario}/{replicate}/{scenario}.ops",
        done="results/beast/runs/{scenario}/{replicate}/run.done",
    params:
        beast_binary=BEAST_BINARY,
        output_prefix=lambda wildcards: f"results/beast/runs/{wildcards.scenario}/{wildcards.scenario}",
        replicate_dir=lambda wildcards: f"results/beast/runs/{wildcards.scenario}/{wildcards.replicate}",
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
        run_r1="results/beast/runs/{scenario}/r1/run.done",
        run_r2="results/beast/runs/{scenario}/r2/run.done",
    output:
        done="results/beast/runs/{scenario}/run.done",
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
        scenario_done="results/beast/runs/strict_constant/run.done",
        log_r1="results/beast/runs/strict_constant/r1/strict_constant.log",
        log_r2="results/beast/runs/strict_constant/r2/strict_constant.log",
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
        scenario_done="results/beast/runs/strict_constant/run.done",
        trees_r1="results/beast/runs/strict_constant/r1/strict_constant.trees",
        trees_r2="results/beast/runs/strict_constant/r2/strict_constant.trees",
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
        scenario_done="results/beast/runs/strict_constant_lugar/run.done",
        log_r1="results/beast/runs/strict_constant_lugar/r1/strict_constant_lugar.log",
        log_r2="results/beast/runs/strict_constant_lugar/r2/strict_constant_lugar.log",
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
        scenario_done="results/beast/runs/strict_constant_lugar/run.done",
        trees_r1="results/beast/runs/strict_constant_lugar/r1/strict_constant_lugar.trees",
        trees_r2="results/beast/runs/strict_constant_lugar/r2/strict_constant_lugar.trees",
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
