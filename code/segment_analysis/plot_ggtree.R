#!/usr/bin/env Rscript
# plot_ggtree.R
# Phylogenetic tree with triangle markers on flu_costa / flu_sierra tips only.
#
# Usage: Rscript plot_ggtree.R <tree_file> <metadata_csv> <output_png> [title]

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(phytools)
  library(ggtree)
  library(treeio)
  library(ggplot2)
})

source("code/segment_analysis/tree_aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_ggtree.R <tree_file> <metadata_csv> <output_file> [title]")
}

tree_file     <- args[1]
metadata_file <- args[2]
output_file   <- args[3]
title         <- if (length(args) >= 4) args[4] else "Phylogenetic Tree"

if (!file.exists(tree_file))     stop(paste("Tree file not found:", tree_file))
if (!file.exists(metadata_file)) stop(paste("Metadata file not found:", metadata_file))

tree <- read.newick(tree_file)
tree <- midpoint.root(tree)

metadata <- read_csv(metadata_file, show_col_types = FALSE) |>
  select(file_name, expected_role) |>
  distinct(file_name, .keep_all = TRUE)

metadata <- metadata[metadata$file_name %in% tree$tip.label, ]

tip_df <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) |>
  left_join(metadata, by = c("label" = "file_name")) |>
  mutate(
    type_c = if_else(
      expected_role %in% FLU_TIP_ROLES,
      expected_role,
      NA_character_
    )
  )

flu_tips <- tip_df |> filter(!is.na(type_c))

p <- ggtree(tree, layout = "roundrect", color = "grey40", linewidth = 0.8) %<+% tip_df

if (nrow(flu_tips) > 0) {
  present <- intersect(FLU_TIP_ROLES, unique(flu_tips$type_c))
  p <- p +
    geom_tippoint(
      aes(subset = !is.na(type_c), color = type_c, shape = type_c),
      size = 2.5
    ) +
    scale_color_manual(
      values = flu_tip_colors[present],
      labels = flu_tip_labels[present],
      name = "Ecuador",
      na.value = NA
    ) +
    scale_shape_manual(
      values = flu_tip_shapes[present],
      labels = flu_tip_labels[present],
      name = "Ecuador",
      na.value = NA
    )
}

p <- p +
  theme_tree2() +
  ggtitle(title) +
  theme(
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position  = "right",
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11, face = "bold"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
n_tips <- Ntip(tree)
fig_w <- 14
fig_h <- min(80, max(8, n_tips * 0.012))
ggsave(output_file, plot = p, width = fig_w, height = fig_h, dpi = 300, limitsize = FALSE)

rds_file <- sub("\\.[^.]+$", ".rds", output_file)
saveRDS(p, rds_file)

cat("Saved tree plot →", output_file, "(", n_tips, "tips,", fig_w, "x", fig_h, "in)\n")
cat("Saved ggtree object →", rds_file, "\n")
