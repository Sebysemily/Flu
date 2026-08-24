#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Generate Snakemake SVG Graphs and HTML Report
# ==============================================================================
# Outputs:
#   1. pipeline_rulegraph.svg                      (Main pipeline without 'all' node)
#   2. pipeline_all_pipeline_options_rulegraph.svg (All options including optional calls)
#   3. pipeline_filegraph.svg                      (File dependencies graph)
#   4. snakemake_report.html                       (Snakemake interactive HTML report)
#
# Usage:
#   ./scripts/gen_snakemake_graphs.sh [TARGET] [OUTDIR]
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_OUTDIR="$(cd "$SCRIPT_DIR/.." && pwd)"
WORKSPACE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
SNAKEFILE="$WORKSPACE_DIR/snakefile"

TARGET="${1:-all}"
OUTDIR="${2:-$DEFAULT_OUTDIR}"
mkdir -p "$OUTDIR"

# Ensure conda env with snakemake / dot is in PATH
if ! command -v snakemake >/dev/null 2>&1; then
  for candidate in \
    "$HOME/miniforge3/envs/snakemake/bin" \
    "$HOME/miniconda3/envs/snakemake/bin" \
    "$HOME/anaconda3/envs/snakemake/bin" \
    "/opt/conda/envs/snakemake/bin"; do
    if [[ -x "$candidate/snakemake" ]]; then
      export PATH="$candidate:$PATH"
      break
    fi
  done
fi

if ! command -v snakemake >/dev/null 2>&1; then
  echo "Error: 'snakemake' not found in PATH. Activate your conda environment." >&2
  exit 1
fi

if ! command -v dot >/dev/null 2>&1; then
  echo "Error: 'dot' (Graphviz) not found in PATH." >&2
  exit 1
fi

if [[ ! -f "$SNAKEFILE" ]]; then
  echo "Error: Snakefile not found at $SNAKEFILE" >&2
  exit 1
fi

echo "============================================================"
echo "Snakemake Workspace: $WORKSPACE_DIR"
echo "Snakefile:           $SNAKEFILE"
echo "Target:              $TARGET"
echo "Output Directory:    $OUTDIR"
echo "============================================================"

# ------------------------------------------------------------------------------
# Helper: Generate DOT, filter 'all' node if requested, and convert to SVG
# ------------------------------------------------------------------------------
generate_graph_svg() {
  local mode="$1"
  local svg_out="$2"
  local strip_all="$3"
  shift 3
  local targets=("$@")

  local tmp_raw
  tmp_raw=$(mktemp)
  local tmp_dot
  tmp_dot=$(mktemp)

  snakemake --snakefile "$SNAKEFILE" --directory "$WORKSPACE_DIR" -n "${targets[@]}" "$mode" > "$tmp_raw" 2>&1 || {
    echo "Error: Snakemake failed generating $mode for targets: ${targets[*]}" >&2
    cat "$tmp_raw" >&2
    rm -f "$tmp_raw" "$tmp_dot"
    exit 1
  }

  python3 - "$tmp_raw" "$tmp_dot" "$strip_all" <<'PY'
import sys
import re

raw_file = sys.argv[1]
out_file = sys.argv[2]
strip_all = (sys.argv[3] == "1")

with open(raw_file, "r", encoding="utf-8", errors="replace") as fh:
    lines = fh.readlines()

digraph_lines = []
keep = False
for line in lines:
    if line.lstrip().startswith("digraph"):
        keep = True
    if keep:
        digraph_lines.append(line)

if not digraph_lines:
    sys.exit("Error: No digraph block found in Snakemake output.")

if strip_all:
    # Find node IDs with label="all"
    all_ids = set()
    for line in digraph_lines:
        m = re.search(r'^\s*"?([0-9A-Za-z_]+)"?\s*\[.*label\s*=\s*"all"', line)
        if m:
            all_ids.add(m.group(1))

    if all_ids:
        filtered = []
        for line in digraph_lines:
            # Drop node definition of 'all'
            if any(re.match(r'^\s*"?%s"?\s*\[' % re.escape(aid), line) for aid in all_ids):
                continue
            # Drop edges pointing to or from 'all'
            skip = False
            for aid in all_ids:
                if re.search(r'(?:^|\s)"?%s"?\s*->' % re.escape(aid), line) or re.search(r'->\s*"?%s"?(?:\s*$|\s|;)' % re.escape(aid), line):
                    skip = True
                    break
            if skip:
                continue
            filtered.append(line)
        digraph_lines = filtered

with open(out_file, "w", encoding="utf-8") as fh:
    fh.writelines(digraph_lines)
PY

  dot -Tsvg "$tmp_dot" -o "$svg_out"
  echo "  [OK] Created: $(basename "$svg_out")"

  rm -f "$tmp_raw" "$tmp_dot"
}

# ------------------------------------------------------------------------------
# 1. Main Pipeline Rulegraph (without 'all' node)
# ------------------------------------------------------------------------------
echo ""
echo "[1/4] Generating Main Pipeline Rulegraph..."
generate_graph_svg --rulegraph "$OUTDIR/pipeline_rulegraph.svg" 1 "$TARGET"

# ------------------------------------------------------------------------------
# 2. All Options Rulegraph (including optional branches, without 'all' node)
# ------------------------------------------------------------------------------
echo ""
echo "[2/4] Generating All Pipeline Options Rulegraph..."

mapfile -t SNAKEFILE_RULES < <(awk '/^rule[[:space:]]+[A-Za-z0-9_]+:/{gsub(":", "", $2); print $2}' "$SNAKEFILE")

ALL_OPTIONS_TARGETS=("$TARGET")
for rule_name in "${SNAKEFILE_RULES[@]}"; do
  if [[ "$rule_name" != "all" ]]; then
    ALL_OPTIONS_TARGETS+=("$rule_name")
  fi
done

OPTIONAL_TARGETS=(
  "plot_flu_lineage"
  "segment_analysis_composite_lineage_preview"
  "plot_panel_sampling_calendar"
  "report_panel_ha_trim_discarded"
  "results/qc_metrics/trim_gap_n/HA_eliminated.csv"
  "results/beast/GSS/model_selection.csv"
)

for opt in "${OPTIONAL_TARGETS[@]}"; do
  ALL_OPTIONS_TARGETS+=("$opt")
done

generate_graph_svg --rulegraph "$OUTDIR/pipeline_all_pipeline_options_rulegraph.svg" 1 "${ALL_OPTIONS_TARGETS[@]}"

# ------------------------------------------------------------------------------
# 3. Pipeline Filegraph
# ------------------------------------------------------------------------------
echo ""
echo "[3/4] Generating Pipeline Filegraph..."
generate_graph_svg --filegraph "$OUTDIR/pipeline_filegraph.svg" 0 "$TARGET"

# ------------------------------------------------------------------------------
# 4. Snakemake HTML Report
# ------------------------------------------------------------------------------
echo ""
echo "[4/4] Generating Snakemake HTML Report..."
snakemake --snakefile "$SNAKEFILE" --directory "$WORKSPACE_DIR" --report "$OUTDIR/snakemake_report.html" "$TARGET" 2>/dev/null || {
  echo "  [INFO] Snakemake HTML report skipped or requires executed outputs."
}
if [[ -f "$OUTDIR/snakemake_report.html" ]]; then
  echo "  [OK] Created: snakemake_report.html"
fi

echo ""
echo "============================================================"
echo "Done! Final files in: $OUTDIR"
echo "============================================================"
ls -lh "$OUTDIR"/*.svg "$OUTDIR"/*.html 2>/dev/null || true
