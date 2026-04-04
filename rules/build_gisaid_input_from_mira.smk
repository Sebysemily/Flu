import os
import glob

MIRA_BASE = config["mira_base_dir"]

MIRA_FASTAS = sorted(glob.glob(os.path.join(MIRA_BASE, "run*", "amended_consensus.fasta")))
MIRA_FASTAS += sorted(glob.glob(os.path.join(MIRA_BASE, "run_agro", "amended_consensus.fasta")))
MIRA_FASTAS = sorted(set(MIRA_FASTAS))

FILTRADO_CSV = config.get("flu_filtrado", "config/flu_filtrado.csv")
ECUADOR_DATE_SOURCE = config.get("ecuador_date_source", "reception")


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