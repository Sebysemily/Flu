#!/usr/bin/env Rscript
# plot_composite_lineage.R
# Composite: flu_lineage heatmap (unchanged) + HA/PB2 ribbon tanglegram.
#
# Usage:
#   Rscript plot_composite_lineage.R <heatmap_rds> <ha_tree> <pb2_tree> \
#       <panel_metadata_csv> <output_png> [max_tips] [ribbon_segment]
#
# Optional max_tips: style-preview subset (keeps all Ecuador costa/sierra tips,
# fills remaining slots with context taxa). Omit for the full panel.
# Optional ribbon_segment: HA or PB2 — lineage column for context tip node colors.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(png)
  library(ape)
  library(dendextend)
  library(dplyr)
  library(scales)
  library(stringr)
})

source("code/segment_analysis/tree_aesthetics.R")
source("code/segment_analysis/lineage_palette.R")
source("code/segment_analysis/tanglegram_legend.R")

select_preview_tips <- function(all_tips, roles, cap) {
  if (is.na(cap) || length(all_tips) <= cap) {
    return(all_tips)
  }

  ecuador_tips <- all_tips[
    all_tips %in% names(roles) & vapply(roles[all_tips], is_ecuador_tip, logical(1))
  ]
  context_tips <- setdiff(all_tips, ecuador_tips)

  if (length(ecuador_tips) >= cap) {
    return(ecuador_tips[seq_len(cap)])
  }

  n_context <- cap - length(ecuador_tips)
  context_pick <- if (length(context_tips) <= n_context) {
    context_tips
  } else {
    sample(context_tips, n_context)
  }

  c(ecuador_tips, context_pick)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript plot_composite_lineage.R ",
    "<heatmap_rds> <ha_tree> <pb2_tree> <panel_metadata_csv> <output_png> ",
    "[max_tips] [ribbon_segment]"
  )
}

heatmap_rds <- args[1]
ha_tree_path <- args[2]
pb2_tree_path <- args[3]
panel_metadata_path <- args[4]
output_png <- args[5]
max_tips <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else NA_integer_
ribbon_segment <- if (length(args) >= 7 && nzchar(args[7])) toupper(args[7]) else "HA"
if (!is.na(max_tips) && max_tips < 2) {
  stop("max_tips must be at least 2")
}
if (!ribbon_segment %in% c("HA", "PB2")) {
  stop("ribbon_segment must be HA or PB2")
}
lineage_col <- paste0(ribbon_segment, "_lineage")

for (path in c(heatmap_rds, ha_tree_path, pb2_tree_path, panel_metadata_path)) {
  if (!file.exists(path)) stop("File not found: ", path)
}

heatmap_obj <- readRDS(heatmap_rds)
p_heatmap <- heatmap_obj$plot
lineage_palette <- heatmap_obj$palette

panel_meta <- read.csv(panel_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"file_name" %in% colnames(panel_meta)) {
  stop("panel metadata must contain file_name column")
}
if (!"expected_role" %in% colnames(panel_meta)) {
  stop("panel metadata must contain expected_role column")
}
if (!"HA_lineage" %in% colnames(panel_meta)) {
  stop("panel metadata missing column: HA_lineage")
}
if (!"PB2_lineage" %in% colnames(panel_meta)) {
  stop("panel metadata missing column: PB2_lineage")
}

tip_roles <- setNames(panel_meta$expected_role, panel_meta$file_name)
tip_lineages_ha <- setNames(
  str_trim(panel_meta$HA_lineage),
  panel_meta$file_name
)
tip_lineages_ha[is.na(tip_lineages_ha) | tip_lineages_ha == ""] <- "unknown"

tip_lineages_pb2 <- setNames(
  str_trim(panel_meta$PB2_lineage),
  panel_meta$file_name
)
tip_lineages_pb2[is.na(tip_lineages_pb2) | tip_lineages_pb2 == ""] <- "unknown"

tree_ha <- read.tree(ha_tree_path)
tree_pb2 <- read.tree(pb2_tree_path)

common_tips <- intersect(tree_ha$tip.label, tree_pb2$tip.label)
if (length(common_tips) < 2) {
  stop("Need at least 2 shared tips between HA and PB2 trees")
}

n_all_tips <- length(common_tips)
common_tips <- select_preview_tips(common_tips, tip_roles, max_tips)
if (!is.na(max_tips)) {
  cat(
    "Style preview:", length(common_tips), "tips",
    "(subset of", n_all_tips, "; all Ecuador costa/sierra kept)\n"
  )
}

tree_ha <- keep.tip(tree_ha, common_tips)
tree_pb2 <- keep.tip(tree_pb2, common_tips)

phylo_to_dend <- function(tree) {
  dist_mat <- cophenetic(tree)
  as.dendrogram(hclust(as.dist(dist_mat), method = "average"))
}

decorate_tangle_leaves <- function(dend, roles, lineages, palette) {
  labs <- labels(dend)
  n <- length(labs)
  pch <- rep(CONTEXT_TIP_SHAPE, n)
  col <- rep("grey80", n)
  cex <- rep(CONTEXT_TIP_CEX, n)
  label_cex <- rep(0, n)

  for (i in seq_along(labs)) {
    tip <- labs[i]
    role <- roles[tip]
    if (is_ecuador_tip(role)) {
      pch[i] <- flu_tip_shapes[role]
      col[i] <- unname(flu_tip_colors[role])
      cex[i] <- FLU_TIP_CEX
    } else {
      pch[i] <- CONTEXT_TIP_SHAPE
      col[i] <- lineage_color_for(lineages[tip], palette)
      cex[i] <- CONTEXT_TIP_CEX
    }
  }

  dend %>%
    set("leaves_pch", pch) %>%
    set("leaves_col", col) %>%
    set("leaves_cex", cex) %>%
    set("labels_cex", label_cex)
}

context_lineages_in_panel <- function(tips, roles, lineages_ha, lineages_pb2) {
  context_tips <- tips[!vapply(roles[tips], is_ecuador_tip, logical(1))]
  vals <- c(lineages_ha[context_tips], lineages_pb2[context_tips])
  unique(vals[!is.na(vals) & vals != ""])
}

dend_ha <- decorate_tangle_leaves(
  phylo_to_dend(tree_ha),
  tip_roles,
  tip_lineages_ha,
  lineage_palette
)
dend_pb2 <- decorate_tangle_leaves(
  phylo_to_dend(tree_pb2),
  tip_roles,
  tip_lineages_pb2,
  lineage_palette
)

dl <- dendlist(dend_ha, dend_pb2)
if (length(common_tips) <= 500) {
  dl <- dl %>% untangle(method = "step1side")
}

tip_order <- labels(dl[[1]])
all_present_lineages <- c(tip_lineages_ha[tip_order], tip_lineages_pb2[tip_order])
ribbon_palette <- extend_lineage_palette(
  lineage_palette,
  all_present_lineages
)

ribbon_colors <- vapply(
  tip_order,
  function(tip) ribbon_color_for_type(tip_roles[tip]),
  character(1)
)

# Do not strip dendrogram labels before tanglegram plotting.
# tanglegram() uses the labels to correctly match leaves between HA and PB2
# to draw the connecting ribbons. The labels themselves are made invisible
# via lab.cex = 0.
dl_plot <- dl

n_tips <- length(common_tips)
tangle_height_in <- if (is.na(max_tips)) {
  42
} else {
  max(8, min(24, n_tips * 0.18))
}
composite_height_in <- if (is.na(max_tips)) 14 else max(10, 3.2 + tangle_height_in * 0.28)

context_lineages <- context_lineages_in_panel(tip_order, tip_roles, tip_lineages_ha, tip_lineages_pb2)
cat(
  "Ribbon colors by expected_role;",
  "context tip nodes by HA & PB2 lineage\n"
)
cat(
  "  types:", paste(names(table(tip_roles[tip_order])), collapse = ", "), "\n"
)
cat(
  "  context lineages:", length(context_lineages), "\n"
)

tangle_tmp <- tempfile(fileext = ".png")
on.exit(unlink(tangle_tmp), add = TRUE)

cat("Rendering tanglegram (", length(common_tips), " tips )...\n")
grDevices::png(
  filename = tangle_tmp,
  width = 16,
  height = tangle_height_in,
  units = "in",
  res = 300,
  bg = "white"
)
tanglegram(
  dl_plot,
  color_lines = ribbon_colors,
  lwd = 6,
  edge.lwd = 1.5,
  columns_width = c(4, 1.2, 4),
  highlight_distinct_edges = FALSE,
  common_subtrees_color_lines = FALSE,
  main_left = "HA",
  main_right = "PB2",
  margin_inner = 0.15,
  margin_outer = 0.1,
  dLeaf = 0,
  lab.cex = 0,
  axes = FALSE,
  left_dendo_mar = c(0.4, 0.05, 0.4, 0.05),
  right_dendo_mar = c(0.4, 0.05, 0.4, 0.05)
)
grDevices::dev.off()

if (!file.exists(tangle_tmp) || file.info(tangle_tmp)$size < 1000) {
  stop("Tanglegram PNG was not created or is empty")
}

tangle_img <- png::readPNG(tangle_tmp)
p_tanglegram <- cowplot::ggdraw() +
  cowplot::draw_grob(
    grid::rasterGrob(
      tangle_img,
      interpolate = TRUE,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc")
    ),
    0, 0, 1, 1
  )

p_tangle_legend <- build_tanglegram_legend(
  ribbon_roles = tip_roles[tip_order],
  lineage_col_label = "lineage",
  context_lineages = context_lineages,
  lineage_palette = ribbon_palette
)

composite_plot <- p_heatmap /
  (p_tanglegram + p_tangle_legend + patchwork::plot_layout(widths = c(4, 1))) +
  patchwork::plot_layout(heights = c(1, 3)) +
  patchwork::plot_annotation(tag_levels = "A")

dir.create(dirname(output_png), showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = output_png,
  plot = composite_plot,
  width = 16,
  height = composite_height_in,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

cat(
  "Composite saved to", output_png,
  "(", length(common_tips), "shared tips; ribbons = type, nodes = HA & PB2 lineages)\n"
)
