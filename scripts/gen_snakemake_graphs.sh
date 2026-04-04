#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scripts/gen_snakemake_graphs.sh [TARGET] [OUTDIR]
# Example (full pipeline, including ML trees): ./scripts/gen_snakemake_graphs.sh all analysis
# Example (only until final FASTA): ./scripts/gen_snakemake_graphs.sh data/final/H5N1_final.fasta analysis

TARGET=${1:-all}
OUTDIR=${2:-analysis}
mkdir -p "$OUTDIR"

echo "Generating Snakemake DOT graphs for target: $TARGET"

resolve_rule_inputs() {
  local rule_name="$1"
  local tmp_rule
  local input_line
  tmp_rule=$(mktemp)

  snakemake -n -p "$rule_name" > "$tmp_rule"

  input_line=$(awk -v rule="rule ${rule_name}:" '
    $0 == rule {in_rule=1; next}
    in_rule && /^rule /{in_rule=0}
    in_rule && /^[[:space:]]*input:/{sub(/^[[:space:]]*input:[[:space:]]*/, ""); print; exit}
  ' "$tmp_rule")
  rm -f "$tmp_rule"

  if [[ -z "${input_line:-}" ]]; then
    return 1
  fi

  printf '%s\n' "$input_line" | tr ',' '\n' | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//' | sed '/^$/d'
}

resolve_targets_for_selection() {
  local selected_target="$1"

  if [[ "$selected_target" == "all" ]]; then
    if ! mapfile -t GRAPH_TARGETS < <(resolve_rule_inputs all); then
      echo "Error: could not resolve 'rule all' inputs dynamically from Snakemake output" >&2
      exit 1
    fi
  else
    GRAPH_TARGETS=("$selected_target")
  fi

  if [[ ${#GRAPH_TARGETS[@]} -eq 0 ]]; then
    echo "Error: resolved empty target list" >&2
    exit 1
  fi
}

resolve_targets_for_selection "$TARGET"

echo "Graph targets: ${GRAPH_TARGETS[*]}"

# Produce DOT files (no execution). Some Snakemake versions print status lines
# before the digraph; keep only the DOT block so Graphviz parses it.
generate_dot() {
  local mode="$1"
  local out="$2"
  local tmp
  tmp=$(mktemp)

  snakemake -n "${GRAPH_TARGETS[@]}" "$mode" > "$tmp"

  awk 'BEGIN{keep=0} /^digraph /{keep=1} keep{print}' "$tmp" > "$out"
  rm -f "$tmp"

  if [[ ! -s "$out" ]]; then
    echo "Error: could not extract DOT content for $mode" >&2
    exit 1
  fi
}

generate_dot --dag "$OUTDIR/pipeline_dag.dot"
generate_dot --rulegraph "$OUTDIR/pipeline_rulegraph.dot"
generate_dot --filegraph "$OUTDIR/pipeline_filegraph.dot"

# Build an additional "all pipeline options" graph from all top-level rules
# defined in snakefile (e.g., all, build_gisaid, rf_validation).
mapfile -t SNAKEFILE_RULES < <(awk '/^rule[[:space:]]+[A-Za-z0-9_]+:/{gsub(":", "", $2); print $2}' snakefile)

ALL_OPTIONS_TARGETS=()
for rule_name in "${SNAKEFILE_RULES[@]}"; do
  if mapfile -t _rule_targets < <(resolve_rule_inputs "$rule_name"); then
    ALL_OPTIONS_TARGETS+=("${_rule_targets[@]}")
  fi
done

if [[ ${#ALL_OPTIONS_TARGETS[@]} -gt 0 ]]; then
  mapfile -t ALL_OPTIONS_TARGETS < <(printf '%s\n' "${ALL_OPTIONS_TARGETS[@]}" | awk '!seen[$0]++')

  generate_all_options_dot() {
    local mode="$1"
    local out="$2"
    local tmp
    tmp=$(mktemp)

    snakemake -n "${ALL_OPTIONS_TARGETS[@]}" "$mode" > "$tmp"
    awk 'BEGIN{keep=0} /^digraph /{keep=1} keep{print}' "$tmp" > "$out"
    rm -f "$tmp"

    if [[ ! -s "$out" ]]; then
      echo "Error: could not extract DOT content for all pipeline options ($mode)" >&2
      exit 1
    fi
  }

  generate_all_options_dot --dag "$OUTDIR/pipeline_all_pipeline_options_dag.dot"
  generate_all_options_dot --rulegraph "$OUTDIR/pipeline_all_pipeline_options_rulegraph.dot"
  generate_all_options_dot --filegraph "$OUTDIR/pipeline_all_pipeline_options_filegraph.dot"
fi

# Convert to SVG (requires Graphviz dot)
if command -v dot >/dev/null 2>&1; then
  echo "Converting DOT -> SVG using dot"
  dot -Tsvg "$OUTDIR/pipeline_dag.dot" -o "$OUTDIR/pipeline_dag.svg"
  dot -Tsvg "$OUTDIR/pipeline_rulegraph.dot" -o "$OUTDIR/pipeline_rulegraph.svg"
  dot -Tsvg "$OUTDIR/pipeline_filegraph.dot" -o "$OUTDIR/pipeline_filegraph.svg"
  if [[ -s "$OUTDIR/pipeline_all_pipeline_options_dag.dot" ]]; then
    dot -Tsvg "$OUTDIR/pipeline_all_pipeline_options_dag.dot" -o "$OUTDIR/pipeline_all_pipeline_options_dag.svg"
    dot -Tsvg "$OUTDIR/pipeline_all_pipeline_options_rulegraph.dot" -o "$OUTDIR/pipeline_all_pipeline_options_rulegraph.svg"
    dot -Tsvg "$OUTDIR/pipeline_all_pipeline_options_filegraph.dot" -o "$OUTDIR/pipeline_all_pipeline_options_filegraph.svg"
  fi
else
  echo "Warning: 'dot' not found in PATH. DOT files created but SVG conversion skipped. Install graphviz to convert."
fi

# Summaries and HTML report
echo "Generating summaries and HTML report"
snakemake -n "${GRAPH_TARGETS[@]}" --summary > "$OUTDIR/pipeline_summary.tsv"
snakemake -n "${GRAPH_TARGETS[@]}" --detailed-summary > "$OUTDIR/pipeline_detailed_summary.tsv"

# Snakemake report (may run lightweight checks)
snakemake --report "$OUTDIR/snakemake_report.html" "${GRAPH_TARGETS[@]}" || echo "--report failed (requires full snakemake environment), DOT + summaries still available."

echo "Done. Outputs in: $OUTDIR"
ls -lh "$OUTDIR"
