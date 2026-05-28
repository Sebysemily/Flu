#!/usr/bin/env Rscript
# plot_ggtree.R
# Plots a phylogenetic tree with tip colours from unified metadata/H5N1_context.csv
#
# CLI args (supplied by Snakemake):
#   1  tree_file     – path to Newick tree
#   2  metadata_file – path to metadata/H5N1_context.csv
#   3  output_file   – path for the PNG output
#   4  title         – (optional) figure title

suppressPackageStartupMessages({
  library(phytools)
  library(tidyverse)
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

metadata <- read_csv(metadata_file, show_col_types = FALSE) |>
  select(file_name, collection_date, country, expected_role) |>
  distinct(file_name, .keep_all = TRUE)

extra_in_metadata <- setdiff(metadata$file_name, tree$tip.label)
if (length(extra_in_metadata) > 0) {
  message("Dropping metadata rows not in tree: ", paste(extra_in_metadata, collapse = ", "))
}
metadata <- metadata[metadata$file_name %in% tree$tip.label, ]

tip_df <- tibble(label = tree$tip.label) |>
  left_join(metadata, by = c("label" = "file_name")) |>
  mutate(
    type_c = if_else(is.na(expected_role), "unknown", expected_role)
  )

present_types <- unique(tip_df$type_c)
aes <- filter_aesthetics_for_types(present_types)
type_colors <- aes$colors
type_shapes <- aes$shapes
type_sizes  <- aes$sizes
type_alphas <- aes$alphas
type_labels <- aes$labels

flu_tips <- tip_df |>
  filter(str_starts(type_c, "flu_")) |>
  pull(label)

green_nodes <- integer(0)
blue_nodes  <- integer(0)

if (length(flu_tips) > 1) {
  flu_mrca  <- getMRCA(tree, flu_tips)
  green_nodes <- flu_mrca

  n_tips    <- Ntip(tree)
  root_node <- n_tips + 1L
  cur <- flu_mrca
  for (lvl in 1:3) {
    par <- tree$edge[tree$edge[, 2] == cur, 1]
    if (length(par) == 0 || par == root_node) break
    children  <- tree$edge[tree$edge[, 1] == par, 2]
    blue_nodes <- c(blue_nodes, children[children != cur])
    cur <- par
  }

  red_flu <- tip_df |>
    filter(type_c %in% c("flu_costa", "flu_epi_isl")) |>
    pull(label)
  if (length(red_flu) > 1) {
    red_mrca <- getMRCA(tree, red_flu)
    if (!red_mrca %in% green_nodes) {
      green_nodes <- c(green_nodes, red_mrca)
    }
  }
}

if (length(green_nodes) + length(blue_nodes) > 0) {
  marker_df <- tibble(
    node   = c(green_nodes, blue_nodes),
    marker = c(rep("Flu clade", length(green_nodes)),
               rep("Sister clade", length(blue_nodes)))
  )
} else {
  marker_df <- tibble(node = integer(0), marker = character(0))
}

p <- ggtree(tree, color = "grey40", linewidth = 0.3) %<+% tip_df %<+% marker_df

if (nrow(marker_df) > 0) {
  p <- p +
    geom_point2(
      aes(subset = !is.na(marker), fill = marker),
      shape = 22, size = 6, color = "grey20", stroke = 1.2
    ) +
    scale_fill_manual(values = mrk_fill, name = "Clade markers")
}

for (tc in names(type_colors)) {
  tc_data <- tip_df |> filter(type_c == tc)
  if (nrow(tc_data) == 0) next
  p <- p +
    geom_tippoint(
      data       = . %>% filter(label %in% tc_data$label),
      aes(color  = type_c, shape = type_c),
      size  = type_sizes[tc],
      alpha = type_alphas[tc]
    )
}

p <- p +
  scale_color_manual(values = type_colors, labels = type_labels, name = "Samples") +
  scale_shape_manual(values = type_shapes, labels = type_labels, name = "Samples") +
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
fig_h <- max(6, Ntip(tree) * 0.10)
ggsave(output_file, plot = p, width = 14, height = fig_h, dpi = 300, limitsize = FALSE)
cat(paste("Saved tree plot →", output_file, "\n"))
