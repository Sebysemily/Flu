#!/usr/bin/env Rscript
# plot_ggtree.R
# -------------
# Plots a phylogenetic tree with tip colours derived from the metadata.
#
# CLI args (supplied by Snakemake):
#   1  tree_file     – path to Newick tree
#   2  metadata_file – path to metadata/H5N1_context.csv
#   3  ecuador_metadata_file - path to config/flu_filtrado.csv
#   4  output_file   – path for the PNG output
#   5  title         – (optional) figure title

suppressPackageStartupMessages({
  library(phytools)   # read.newick
  library(tidyverse)  # dplyr, tibble, etc.
  library(ggtree)     # ggtree, geom_tippoint, …
  library(treeio)
})

# ── CLI args ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_ggtree.R <tree_file> <metadata_csv> <ecuador_metadata_csv> <output_file> [title]")
}

tree_file             <- args[1]
metadata_file         <- args[2]
ecuador_metadata_file <- args[3]
output_file           <- args[4]
title                 <- if (length(args) >= 5) args[5] else "Phylogenetic Tree"

if (!file.exists(tree_file))             stop(paste("Tree file not found:", tree_file))
if (!file.exists(metadata_file))         stop(paste("Metadata file not found:", metadata_file))
if (!file.exists(ecuador_metadata_file)) stop(paste("Ecuador metadata file not found:", ecuador_metadata_file))

# ── Load tree ──────────────────────────────────────────────────────────────
tree <- read.newick(tree_file)

# ── Load and concatenate metadata ──────────────────────────────────────────
# A. Read Ecuador metadata
ecuador_meta <- read_csv(ecuador_metadata_file, show_col_types = FALSE) |>
  select(
    file_name = EPI_ISL,
    collection_date = `Fecha colección`,
    sample_id = `Código USFQ`
  ) |>
  mutate(
    country = "Ecuador",
    expected_role = if_else(sample_id == "Flu-0406" | file_name %in% c("EPI_ISL_17973443", "EPI_ISL_17973458", "EPI_ISL_18137671"), "flu_costa", "flu_sierra")
  ) |>
  select(file_name, collection_date, country, expected_role)

# B. Read GISAID context metadata
context_meta <- read_csv(metadata_file, show_col_types = FALSE) |>
  select(file_name, collection_date, country, expected_role)

# C. Concatenate
metadata <- bind_rows(ecuador_meta, context_meta)

# D. Drop extra metadata rows not in the tree
extra_in_metadata <- setdiff(metadata$file_name, tree$tip.label)
if(length(extra_in_metadata) > 0) {
  message("Dropping metadata rows not in tree: ", 
          paste(extra_in_metadata, collapse=", "))
}
metadata <- metadata[metadata$file_name %in% tree$tip.label, ]

# E. Prepare tip dataframe for plotting
tip_df <- tibble(label = tree$tip.label) |>
  left_join(metadata, by = c("label" = "file_name"))

# Validate that all tree tips have metadata
missing_metadata <- tip_df |> filter(is.na(expected_role))
if (nrow(missing_metadata) > 0) {
  stop(paste("ERROR: Metadata not found for tree tip labels: ",
             paste(missing_metadata$label, collapse=", ")))
}

tip_df <- tip_df |>
  mutate(type_c = expected_role)

# ── Colour / shape palette (mirrors iTOL logic) ────────────────────────────
#   flu_costa        → coastal Ecuador   (#FF0000, red)
#   flu_sierra       → sierra Ecuador    (#00008B, dark-blue)
#   american_anchor                  (#800080, purple)
#   regional_context                 (#4CAF50, green)

type_colors <- c(
  flu_costa        = "#FF0000",
  flu_sierra       = "#00008B",
  american_anchor  = "#800080",
  regional_context = "#4CAF50"
)

type_shapes <- c(
  flu_costa        = 17L,   # filled triangle
  flu_sierra       = 17L,
  american_anchor  = 15L,   # filled square
  regional_context = 16L    # filled circle
)

type_sizes <- c(
  flu_costa        = 2.5,
  flu_sierra       = 2.5,
  american_anchor  = 2.0,
  regional_context = 1.2
)

type_alphas <- c(
  flu_costa        = 1.0,
  flu_sierra       = 1.0,
  american_anchor  = 0.85,
  regional_context = 0.55
)

# Friendly legend labels
type_labels <- c(
  flu_costa        = "Ecuador (coastal)",
  flu_sierra       = "Ecuador (sierra)",
  american_anchor  = "American anchor",
  regional_context = "Regional context"
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
    filter(type_c == "flu_costa") |>
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
