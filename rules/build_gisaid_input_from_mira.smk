import os
import glob

MIRA_BASE = config["mira_base_dir"]
ECUADOR_INPUT_SOURCE = config.get("ecuador_input_source", "assembled")
VALID_ECUADOR_INPUT_SOURCES = {"assembled", "mira_raw"}

if ECUADOR_INPUT_SOURCE not in VALID_ECUADOR_INPUT_SOURCES:
    raise ValueError(
        "config.ecuador_input_source must be one of: "
        + ", ".join(sorted(VALID_ECUADOR_INPUT_SOURCES))
    )

MIRA_CFG = config.get("mira_cli", {})
MIRA_RAW_INPUT = ECUADOR_INPUT_SOURCE == "mira_raw" or bool(MIRA_CFG.get("enabled", False))
MIRA_DATA_DIR = MIRA_CFG.get("data_dir", MIRA_BASE)
MIRA_RUN_GLOB = MIRA_CFG.get("run_glob", "run*")
MIRA_REQUIRE_SAMPLESHEET = MIRA_CFG.get("require_samplesheet", True)


def get_existing_mira_fastas():
    mira_fastas = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
    mira_fastas += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
    return sorted(set(mira_fastas))


def get_mira_cli_runs():
    runs = sorted(glob.glob(os.path.join(MIRA_DATA_DIR, MIRA_RUN_GLOB)))
    selected_runs = []
    for run in runs:
        if not os.path.isdir(run):
            continue
        if os.path.basename(run) == "run_agro":
            continue
        if MIRA_REQUIRE_SAMPLESHEET and not os.path.exists(os.path.join(run, "samplesheet.csv")):
            continue
        selected_runs.append(run)
    return selected_runs


def get_mira_cli_consensus():
    return [os.path.join(run, "amended_consensus.fasta") for run in get_mira_cli_runs()]

MIRA_FASTAS = get_mira_cli_consensus() if MIRA_RAW_INPUT else get_existing_mira_fastas()

FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
ECUADOR_DATE_SOURCE = config.get("ecuador_date_source", "reception")


rule run_mira_cli_for_run:
    input:
        samplesheet=os.path.join(MIRA_DATA_DIR, "{run_id}", "samplesheet.csv")
    output:
        consensus=os.path.join(MIRA_DATA_DIR, "{run_id}", "amended_consensus.fasta"),
        stdout=os.path.join(MIRA_DATA_DIR, "{run_id}", "mira_batch_resume.out"),
        stderr=os.path.join(MIRA_DATA_DIR, "{run_id}", "mira_batch_resume.err")
    params:
        container=MIRA_CFG.get("container_name", "mira"),
        workdir=MIRA_CFG.get("container_workdir", "/data"),
        experiment=MIRA_CFG.get("experiment", "Flu-ONT")
    shell:
        r"""
        set -euo pipefail

        docker ps --format '{{{{.Names}}}}' | grep -qx "{params.container}" || (
            echo "ERROR: MIRA container '{params.container}' is not running." >&2
            echo "Run: cd ~/MIRA_NGS && docker compose up -d" >&2
            exit 1
        )

        docker exec -w {params.workdir} {params.container} MIRA.sh \
            -s {wildcards.run_id}/samplesheet.csv \
            -r {wildcards.run_id} \
            -e {params.experiment} \
            > {output.stdout} \
            2> {output.stderr}

        test -s {output.consensus}
        """


rule build_ecuador_intermediate_input:
    input:
        MIRA_FASTAS,
        FILTRADO_CSV
    output:
        directory("data/all_amended_fasta"),
        "data/assembled/ecuador_intermediate_raw_segments.csv",
        "data/assembled/ecuador_intermediate_expected_lengths.csv",
        "data/assembled/ecuador_intermediate_summary.csv",
        "data/assembled/ecuador_intermediate_audit.csv",
        "data/assembled/ecuador_intermediate_issues.csv",
        "data/assembled/ecuador_intermediate_sequences.fasta",
        directory("data/assembled/ecuador_intermediate_per_sample")
    shell:
        r"""
        python code/build_gisaid_input_from_mira/build_gisaid_input_from_mira.py
        """


rule build_h5n1_ec_fasta:
    input:
        per_sample_fastas="data/assembled/ecuador_intermediate_per_sample",
        audit_csv="data/assembled/ecuador_intermediate_audit.csv",
        metadata_csv=FILTRADO_CSV
    output:
        fasta="data/input/H5N1_EC.fasta",
        summary="data/input/H5N1_EC_summary.csv"
    params:
        ecuador_date_source=ECUADOR_DATE_SOURCE
    shell:
        r"""
        export PYTHONPATH=code:${{PYTHONPATH:-}}
        python code/build_gisaid_input_from_mira/build_denv2_style_fasta_from_assembled.py \
            --per-sample-dir {input.per_sample_fastas} \
            --audit-csv {input.audit_csv} \
            --metadata-csv {input.metadata_csv} \
            --ecuador-date-source {params.ecuador_date_source} \
            --output-fasta {output.fasta} \
            --summary-csv {output.summary}
        """
