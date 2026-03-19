#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scripts/gen_snakemake_graphs.sh [TARGET] [OUTDIR]
# Example: ./scripts/gen_snakemake_graphs.sh data/final/H5N1_final.fasta analysis

TARGET=${1:-data/final/H5N1_final.fasta}
OUTDIR=${2:-analysis}
mkdir -p "$OUTDIR"

echo "Generating Snakemake DOT graphs for target: $TARGET"

# Produce DOT files (no execution). Some Snakemake versions print status lines
# before the digraph; keep only the DOT block so Graphviz parses it.
generate_dot() {
  local mode="$1"
  local out="$2"
  local tmp
  tmp=$(mktemp)

  snakemake -n "$TARGET" "$mode" > "$tmp"

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

# Convert to SVG (requires Graphviz dot)
if command -v dot >/dev/null 2>&1; then
  echo "Converting DOT -> SVG using dot"
  dot -Tsvg "$OUTDIR/pipeline_dag.dot" -o "$OUTDIR/pipeline_dag.svg"
  dot -Tsvg "$OUTDIR/pipeline_rulegraph.dot" -o "$OUTDIR/pipeline_rulegraph.svg"
  dot -Tsvg "$OUTDIR/pipeline_filegraph.dot" -o "$OUTDIR/pipeline_filegraph.svg"
else
  echo "Warning: 'dot' not found in PATH. DOT files created but SVG conversion skipped. Install graphviz to convert."
fi

# Summaries and HTML report
echo "Generating summaries and HTML report"
snakemake -n "$TARGET" --summary > "$OUTDIR/pipeline_summary.tsv"
snakemake -n "$TARGET" --detailed-summary > "$OUTDIR/pipeline_detailed_summary.tsv"

# Snakemake report (may run lightweight checks)
snakemake --report "$OUTDIR/snakemake_report.html" || echo "--report failed (requires full snakemake environment), DOT + summaries still available."

echo "Done. Outputs in: $OUTDIR"
ls -lh "$OUTDIR"
