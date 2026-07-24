#!/usr/bin/env Rscript
# plot_composite_lineage.R
# Composite: flu_lineage heatmap + HA/PB2 cophylogenetic tanglegram (phytools).
#
# Usage:
#   Rscript plot_composite_lineage.R <ha_tree> <pb2_tree> \
#       <panel_metadata_csv> <output_png> [max_tips] [ribbon_segment]
#
# max_tips: preview subset (keeps all Ecuador tips, fills with context).
# ribbon_segment: HA or PB2 — ignored for ribbons (colored by type); kept for CLI compat.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(png)
  library(ape)
  library(phytools)
  library(scales)
  library(stringr)
})

source("code/segment_analysis/tree_aesthetics.R")
source("code/segment_analysis/lineage_palette.R")
source("code/segment_analysis/tanglegram_legend.R")
source("code/segment_analysis/plot_cophylo_tanglegram.R")

select_preview_tips <- function(all_tips, roles, cap, outgroup_sample = NA_character_) {
  if (is.na(cap) || length(all_tips) <= cap) {
    return(all_tips)
  }

  ecuador_tips <- all_tips[
    all_tips %in% names(roles) & vapply(roles[all_tips], is_ecuador_tip, logical(1))
  ]
  
  # Ensure outgroup sample is kept if specified and present
  outgroup_to_keep <- NULL
  if (!is.na(outgroup_sample) && outgroup_sample %in% all_tips) {
    outgroup_to_keep <- outgroup_sample
  }

  context_tips <- setdiff(all_tips, c(ecuador_tips, outgroup_to_keep))

  needed_tips <- c(ecuador_tips, outgroup_to_keep)
  if (length(needed_tips) >= cap) {
    return(needed_tips[seq_len(cap)])
  }

  n_context <- cap - length(needed_tips)
  context_pick <- if (length(context_tips) <= n_context) {
    context_tips
  } else {
    sample(context_tips, n_context)
  }

  c(needed_tips, context_pick)
}

context_lineages_in_panel <- function(tips, roles, lineages_ha, lineages_pb2) {
  context_tips <- tips[!vapply(roles[tips], is_ecuador_tip, logical(1))]
  vals <- c(lineages_ha[context_tips], lineages_pb2[context_tips])
  unique(vals[!is.na(vals) & vals != "" & vals != "unknown"])
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_composite_lineage.R <ha_tree> <pb2_tree> <panel_metadata_csv> <output_png> [max_tips] [ribbon_segment]")
}

ha_tree_path <- args[1]
pb2_tree_path <- args[2]
panel_metadata_path <- args[3]
output_png <- args[4]
max_tips <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_
ribbon_segment <- if (length(args) >= 6) args[6] else "HA"
outgroup_sample <- if (length(args) >= 7 && nzchar(args[7])) args[7] else NA_character_

# We no longer load the heatmap or its palette.
lineage_palette <- character()

for (path in c(ha_tree_path, pb2_tree_path, panel_metadata_path)) {
  if (!file.exists(path)) stop("File not found: ", path)
}

panel_meta <- read.csv(panel_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

tip_roles <- normalize_role_vector(setNames(panel_meta$expected_role, panel_meta$file_name))
tip_genotype <- setNames(str_trim(panel_meta$genotype), panel_meta$file_name)
tip_genotype[is.na(tip_genotype) | tip_genotype == ""] <- "unknown"

tree_ha <- read.tree(ha_tree_path)
tree_pb2 <- read.tree(pb2_tree_path)
common_tips <- intersect(tree_ha$tip.label, tree_pb2$tip.label)

n_all_tips <- length(common_tips)
common_tips <- select_preview_tips(common_tips, tip_roles, max_tips, outgroup_sample)
if (!is.na(max_tips)) {
  cat(
    "Style preview:", length(common_tips), "tips",
    "(subset of", n_all_tips, "; all Ecuador costa/sierra kept)\n"
  )
}

tree_ha <- keep.tip(tree_ha, common_tips)
tree_pb2 <- keep.tip(tree_pb2, common_tips)

# Root trees with outgroup sample if available
if (!is.na(outgroup_sample) && outgroup_sample %in% tree_ha$tip.label) {
  tree_ha <- ape::root(tree_ha, outgroup = outgroup_sample, resolve.root = TRUE)
} else {
  tree_ha <- phytools::midpoint.root(tree_ha)
}

if (!is.na(outgroup_sample) && outgroup_sample %in% tree_pb2$tip.label) {
  tree_pb2 <- ape::root(tree_pb2, outgroup = outgroup_sample, resolve.root = TRUE)
} else {
  tree_pb2 <- phytools::midpoint.root(tree_pb2)
}

tree_ha <- ape::ladderize(tree_ha, right = FALSE)   # largest clades at bottom
tree_pb2 <- ape::ladderize(tree_pb2, right = FALSE)  # largest clades at bottom

genotype_palette <- extend_lineage_palette(
   lineage_palette,
   tip_genotype[common_tips]
 )
genotype_palette["B3.2"] <- "#0000FF"
genotype_palette["b3.2"] <- "#0000FF"
genotype_palette["B4.1"] <- "#00BFFF"
genotype_palette["b4.1"] <- "#00BFFF"

n_tips <- length(common_tips)
tangle_height_in <- if (is.na(max_tips)) 42 else max(8, min(24, n_tips * 0.18))
composite_height_in <- if (is.na(max_tips)) 14 else max(10, 3.2 + tangle_height_in * 0.28)

tangle_tmp <- tempfile(fileext = ".png")
on.exit(unlink(tangle_tmp), add = TRUE)

cat("Rendering cophylo tanglegram (", n_tips, " tips )...\n")
grDevices::png(
  filename = tangle_tmp,
  width = 16,
  height = tangle_height_in,
  units = "in",
  res = 300,
  bg = "white"
)

cophylo_plot <- render_cophylo_tanglegram(
   tree_ha = tree_ha,
   tree_pb2 = tree_pb2,
   tip_roles = tip_roles,
   tip_lineages_ha = tip_genotype,
   tip_lineages_pb2 = tip_genotype,
   lineage_palette = genotype_palette
 )
tip_order <- cophylo_plot$tip_order

grDevices::dev.off()

context_lineages <- unique(tip_genotype[tip_order])

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
  ) +
  cowplot::draw_label(
    label = "HA",
    x = 0.06,
    y = 0.985,
    hjust = 0.5,
    vjust = 1,
    fontface = "bold",
    size = 18
  ) +
  cowplot::draw_label(
    label = "PB2",
    x = 0.90,
    y = 0.985,
    hjust = 0.5,
    vjust = 1,
    fontface = "bold",
    size = 18
  )

p_tangle_legend <- build_tanglegram_legend(
   ribbon_roles = tip_roles[tip_order],
   lineage_col_label = "genotype",
   context_lineages = context_lineages,
   lineage_palette = genotype_palette
 ) +
   ggplot2::theme(plot.tag = ggplot2::element_blank())

composite_plot <- (p_tanglegram + p_tangle_legend + patchwork::plot_layout(widths = c(4, 1)))

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
  "(", n_tips, "tips; cophylo; ribbons = type, nodes = HA/PB2 lineage)\n"
)
