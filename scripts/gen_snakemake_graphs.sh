#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scripts/gen_snakemake_graphs.sh [TARGET] [OUTDIR]
# Example (full pipeline, including ML trees): ./scripts/gen_snakemake_graphs.sh all analysis
# Example (only until final FASTA): ./scripts/gen_snakemake_graphs.sh data/final/H5N1_final.fasta analysis


TARGET=${1:-all}
OUTDIR=${2:-analysis}
mkdir -p "$OUTDIR"

if ! command -v snakemake >/dev/null 2>&1; then
  echo "Error: 'snakemake' not found in PATH. Activate the conda environment that contains snakemake or install it." >&2
  exit 1
fi

echo "Generating Snakemake DOT graphs for target: $TARGET"

GRAPH_TARGETS=("$TARGET")
STRIP_ROOT_ALL=0
if [[ "$TARGET" == "all" ]]; then
  STRIP_ROOT_ALL=1
fi

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

strip_rule_all_node() {
  local file="$1"
  python3 - "$file" <<'PY'
import re
import sys
fname = sys.argv[1]
with open(fname, 'r') as fh:
    lines = fh.readlines()
all_ids = set()
for line in lines:
    m = re.match(r'^\s*"?([0-9A-Za-z_]+)"?\s*\[.*label\s*=\s*"all"', line)
    if m:
        all_ids.add(m.group(1))
if not all_ids:
    # Fallback: remove bare label lines with all
    with open(fname, 'w') as fh:
        for line in lines:
            if 'label="all"' in line or 'label = "all"' in line:
                continue
            fh.write(line)
    sys.exit(0)
new_lines = []
for line in lines:
    if any(re.match(r'^\s*"?%s"?\s*\[' % re.escape(all_id), line) for all_id in all_ids):
        continue
    skip = False
    for all_id in all_ids:
        if re.search(r'(^|\s)"?%s"?\s*->' % re.escape(all_id), line):
            skip = True
            break
        if re.search(r'->\s*"?%s"?(\s*$|\s)' % re.escape(all_id), line):
            skip = True
            break
    if skip:
        continue
    new_lines.append(line)
with open(fname, 'w') as fh:
    fh.writelines(new_lines)
PY
}

generate_dot --dag "$OUTDIR/pipeline_dag.dot"
generate_dot --rulegraph "$OUTDIR/pipeline_rulegraph.dot"
generate_dot --filegraph "$OUTDIR/pipeline_filegraph.dot"

if [[ $STRIP_ROOT_ALL -eq 1 ]]; then
  strip_rule_all_node "$OUTDIR/pipeline_dag.dot"
  strip_rule_all_node "$OUTDIR/pipeline_rulegraph.dot"
  strip_rule_all_node "$OUTDIR/pipeline_filegraph.dot"
fi

# Build an additional "all pipeline options" graph from all top-level rules
# defined in snakefile (exclude the controller rule 'all' to avoid a cluttered final edge).
mapfile -t SNAKEFILE_RULES < <(awk '/^rule[[:space:]]+[A-Za-z0-9_]+:/{gsub(":", "", $2); print $2}' snakefile)
ALL_OPTIONS_TARGETS=()
for rule_name in "${SNAKEFILE_RULES[@]}"; do
  if [[ "$rule_name" == "all" ]]; then
    continue
  fi
  ALL_OPTIONS_TARGETS+=("$rule_name")
 done

if [[ ${#ALL_OPTIONS_TARGETS[@]} -gt 0 ]]; then
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
