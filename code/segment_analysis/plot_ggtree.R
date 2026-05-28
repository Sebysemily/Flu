#!/usr/bin/env Rscript
# plot_ggtree.R
# -------------
# Plots a rooted phylogenetic tree with tip colours derived from the
# main-panel metadata CSV (type_c column), replacing the old iTOL-file
# parsing workflow.
#
# CLI args (supplied by Snakemake):
#   1  tree_file     – path to Newick tree
#   2  metadata_file – path to metadata/main_panel.csv
#   3  output_file   – path for the PNG output
#   4  title         – (optional) figure title

suppressPackageStartupMessages({
  library(phytools)   # read.newick
  library(tidyverse)  # dplyr, tibble, etc.
  library(ggtree)     # ggtree, geom_tippoint, …
  library(ape)        # getMRCA, Ntip
})

# ── CLI args ───────────────────────────────────────────────────────────────
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

# ── Load tree ──────────────────────────────────────────────────────────────
tree <- read.newick(tree_file)

# Root on the American anchor outgroup (midpoint fallback)
outgroup_name <- "ABlue-winged_TealSouth_CarolinaUSDA-000345-0032021_EPI_ISL_18133416__american_anchor/USA/2021-12-30"
if (outgroup_name %in% tree$tip.label) {
  tree <- root(tree, outgroup = outgroup_name, resolve.root = TRUE)
} else {
  message("Outgroup not found – falling back to midpoint root.")
  tree <- midpoint.root(tree)
}

# ── Load metadata ──────────────────────────────────────────────────────────
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
# Ensure all tips are present; fill unknowns gracefully
tip_df <- tibble(label = tree$tip.label) |>
  left_join(meta, by = "label") |>
  mutate(
    type_c = if_else(is.na(type_c), "unknown", type_c)
  )

# ── Colour / shape palette (mirrors iTOL logic) ────────────────────────────
#   flu_epi_isl  → coastal Ecuador   (#FF0000, red)
#   flu_sierra   → sierra Ecuador    (#00008B, dark-blue)
#   american_anchor                  (#800080, purple)
#   regional_context                 (#4CAF50, green)
#   unknown                          (#999999, grey)

type_colors <- c(
  flu_epi_isl      = "#FF0000",
  flu_sierra       = "#00008B",
  american_anchor  = "#800080",
  regional_context = "#4CAF50",
  unknown          = "#999999"
)

type_shapes <- c(
  flu_epi_isl      = 17L,   # filled triangle
  flu_sierra       = 17L,
  american_anchor  = 15L,   # filled square
  regional_context = 16L,   # filled circle
  unknown          = 16L
)

type_sizes <- c(
  flu_epi_isl      = 2.5,
  flu_sierra       = 2.5,
  american_anchor  = 2.0,
  regional_context = 1.2,
  unknown          = 1.0
)

type_alphas <- c(
  flu_epi_isl      = 1.0,
  flu_sierra       = 1.0,
  american_anchor  = 0.85,
  regional_context = 0.55,
  unknown          = 0.4
)

# Friendly legend labels
type_labels <- c(
  flu_epi_isl      = "Ecuador (coastal)",
  flu_sierra       = "Ecuador (sierra)",
  american_anchor  = "American anchor",
  regional_context = "Regional context",
  unknown          = "Unknown"
)

# Only include levels that actually appear in this tree
present_types <- unique(tip_df$type_c)
type_colors  <- type_colors[present_types]
type_shapes  <- type_shapes[present_types]
type_sizes   <- type_sizes[present_types]
type_alphas  <- type_alphas[present_types]
type_labels  <- type_labels[present_types]

# ── Internal-node markers (Flu clade / Sister clade) ──────────────────────
flu_tips <- tip_df |>
  filter(str_starts(type_c, "flu_")) |>
  pull(label)

green_nodes <- integer(0)
blue_nodes  <- integer(0)

if (length(flu_tips) > 1) {
  flu_mrca  <- getMRCA(tree, flu_tips)
  green_nodes <- flu_mrca

  # Walk up 3 levels marking sister branches
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

  # Coastal (red) sub-clade → additional green marker
  red_flu <- tip_df |>
    filter(type_c == "flu_epi_isl") |>
    pull(label)
  if (length(red_flu) > 1) {
    red_mrca <- getMRCA(tree, red_flu)
    if (!red_mrca %in% green_nodes)
      green_nodes <- c(green_nodes, red_mrca)
  }
}

if (length(green_nodes) + length(blue_nodes) > 0) {
  marker_df <- tibble(
    node   = c(green_nodes, blue_nodes),
    marker = c(rep("Flu clade",    length(green_nodes)),
               rep("Sister clade", length(blue_nodes)))
  )
} else {
  marker_df <- tibble(node = integer(0), marker = character(0))
}

mrk_fill <- c("Flu clade" = "#2ca02c", "Sister clade" = "#1f77b4")

# ── Build ggtree plot ──────────────────────────────────────────────────────
p <- ggtree(tree, color = "grey40", linewidth = 0.3) %<+% tip_df %<+% marker_df

# Internal-node clade markers
if (nrow(marker_df) > 0) {
  p <- p +
    geom_point2(
      aes(subset = !is.na(marker), fill = marker),
      shape = 22, size = 6, color = "grey20", stroke = 1.2
    ) +
    scale_fill_manual(values = mrk_fill, name = "Clade markers")
}

# Tip points – one geom per type so sizes / alphas can vary
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
  scale_shape_manual(values = type_shapes, labels = type_labels, name = "Samples")

# Theme
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

# ── Save ───────────────────────────────────────────────────────────────────
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
fig_h <- max(6, Ntip(tree) * 0.10)
ggsave(output_file, plot = p, width = 14, height = fig_h, dpi = 300, limitsize = FALSE)
cat(paste("Saved tree plot →", output_file, "\n"))
