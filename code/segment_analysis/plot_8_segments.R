#!/usr/bin/env Rscript
# plot_8_segments.R
# Plots 8 segment trees as individual panels arranged in 2 rows x 4 columns.
# Each panel: rooted tree, background clades collapsed as triangles,
# flu samples coloured by role, context samples in grey.
#
# Usage:
#   Rscript plot_8_segments.R <meta.csv> <outgroup> <out.png> <tree1> ... <tree8>

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(stringr)
  library(patchwork)
})

source("code/segment_analysis/tree_aesthetics.R")

# ---------------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11) {
  stop("Usage: Rscript plot_8_segments.R <meta.csv> <outgroup> <out.png> <tree1> ... <tree8>")
}

metadata_path  <- args[1]
outgroup_sample <- args[2]
output_png     <- args[3]
tree_paths     <- args[4:11]

# Derive segment names from directory component (e.g. results/iq-tree/HA/HA.treefile -> HA)
segments <- sapply(tree_paths, function(x) {
  parts <- unlist(strsplit(x, "/"))
  parts[length(parts) - 1]
})

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
meta          <- read.csv(metadata_path, stringsAsFactors = FALSE)
flu_tips_all  <- meta$file_name[grepl("^flu_", meta$expected_role)]

# Tip role lookup: file_name -> expected_role
role_lookup <- setNames(meta$expected_role, meta$file_name)

# ---------------------------------------------------------------------------
# Helper: build one ggplot panel for a segment tree
# ---------------------------------------------------------------------------
build_segment_panel <- function(tree_path, segment_name, flu_tips_all,
                                 role_lookup, outgroup_sample) {

  tree <- read.tree(tree_path)

  # Root
  if (!is.na(outgroup_sample) && outgroup_sample %in% tree$tip.label) {
    tree <- ape::root(tree, outgroup = outgroup_sample, resolve.root = TRUE)
  } else {
    tree <- phytools::midpoint.root(tree)
  }
  tree <- ape::ladderize(tree, right = FALSE)

  flu_tips_in_tree <- intersect(tree$tip.label, flu_tips_all)
  n_tips           <- Ntip(tree)

  # Fortify for ggplot
  d <- ggtree::fortify(tree)

  # ------------------------------------------------------------------
  # Identify background clades to collapse (triangle polygons)
  # ------------------------------------------------------------------
  poly_list  <- list()
  poly_id    <- 1
  nodes_to_drop <- c()

  if (length(flu_tips_in_tree) > 0) {
    desc_tips <- phangorn::Descendants(tree, type = "tips")
    desc_all  <- phangorn::Descendants(tree, type = "all")

    for (node in (n_tips + 1):(n_tips + Nnode(tree))) {
      tips_in_node <- tree$tip.label[desc_tips[[node]]]
      if (length(intersect(tips_in_node, flu_tips_in_tree)) == 0) {
        parent <- tree$edge[tree$edge[, 2] == node, 1]
        if (length(parent) > 0) {
          parent_tips <- tree$tip.label[desc_tips[[parent]]]
          if (length(intersect(parent_tips, flu_tips_in_tree)) > 0) {
            desc_nodes <- desc_all[[node]]
            x_root     <- d$x[d$node == node]
            y_root     <- d$y[d$node == node]
            desc_data  <- d[d$node %in% desc_nodes, ]
            if (nrow(desc_data) > 0) {
              poly_list[[poly_id]] <- data.frame(
                x     = c(x_root, max(desc_data$x), max(desc_data$x)),
                y     = c(y_root, max(desc_data$y), min(desc_data$y)),
                group = paste0("poly_", poly_id)
              )
              poly_id       <- poly_id + 1
              nodes_to_drop <- c(nodes_to_drop, desc_nodes)
            }
          }
        }
      }
    }
  }

  # Remove collapsed nodes from tree data
  if (length(nodes_to_drop) > 0) {
    d <- d[!d$node %in% nodes_to_drop, ]
  }

  # ------------------------------------------------------------------
  # Tip annotation: colour & shape by role
  # ------------------------------------------------------------------
  tip_data <- d[d$isTip, ] |>
    left_join(meta, by = c("label" = "file_name")) |>
    mutate(
      is_flu = label %in% flu_tips_in_tree,
      role   = ifelse(is_flu, expected_role, "context")
    )

  # ------------------------------------------------------------------
  # Build ggplot
  # ------------------------------------------------------------------
  p <- ggplot() +
    geom_tree(data = d, color = "grey35", linewidth = 0.25) +
    theme_tree2() +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5,
                                     margin = margin(b = 4)),
      axis.line.x     = element_line(color = "grey40", linewidth = 0.35),
      axis.ticks.x    = element_line(color = "grey40", linewidth = 0.3),
      axis.text.x     = element_text(size = 6, color = "grey30"),
      plot.margin     = margin(6, 8, 6, 4),
      legend.position = "none"
    ) +
    ggtitle(segment_name)

  # Collapsed-clade triangles
  if (length(poly_list) > 0) {
    poly_df <- bind_rows(poly_list)
    p <- p + geom_polygon(
      data  = poly_df,
      aes(x = x, y = y, group = group),
      fill  = "grey82",
      color = "grey55",
      alpha = 0.75,
      linewidth = 0.25
    )
  }

  # Context tips (grey dots, small)
  ctx <- tip_data[!tip_data$is_flu, ]
  if (nrow(ctx) > 0) {
    p <- p + geom_point(
      data  = ctx,
      aes(x = x, y = y),
      color = "grey70",
      shape = 16,
      size  = 0.8,
      alpha = 0.7
    )
  }

  # Ecuador flu tips (coloured triangles, larger)
  flu_df <- tip_data[tip_data$is_flu, ]
  if (nrow(flu_df) > 0) {
    flu_df$role <- normalize_role_vector(setNames(flu_df$expected_role, flu_df$label))
    p <- p + geom_point(
      data  = flu_df,
      aes(x = x, y = y, color = role),
      shape = 17,
      size  = 2.2,
      alpha = 0.95
    ) +
    scale_color_manual(values = panel_type_colors, na.value = "grey70")
  }

  p
}

# ---------------------------------------------------------------------------
# Build all 8 panels
# ---------------------------------------------------------------------------
panels <- vector("list", length(tree_paths))
for (i in seq_along(tree_paths)) {
  cat("Processing", segments[i], "...\n")
  panels[[i]] <- build_segment_panel(
    tree_path       = tree_paths[i],
    segment_name    = segments[i],
    flu_tips_all    = flu_tips_all,
    role_lookup     = role_lookup,
    outgroup_sample = outgroup_sample
  )
}

# ---------------------------------------------------------------------------
# Shared legend
# ---------------------------------------------------------------------------
legend_roles <- names(panel_type_colors)   # flu_costa, flu_andine, flu_amazonia, american_anchor, regional_context
legend_shapes <- c(flu_costa = 17, flu_andine = 17, flu_amazonia = 17,
                   american_anchor = 16, regional_context = 16)

legend_data <- data.frame(
  x    = seq_along(legend_roles),
  y    = 1,
  role = legend_roles,
  stringsAsFactors = FALSE
)
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = role, shape = role)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = panel_type_colors,
    labels = panel_type_labels,
    name   = NULL,
    na.value = "grey70"
  ) +
  scale_shape_manual(
    values = legend_shapes,
    labels = panel_type_labels,
    name   = NULL
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.55, "cm"),
    legend.box       = "horizontal",
    legend.margin    = margin(4, 4, 4, 4)
  )

shared_legend <- cowplot::get_legend(legend_plot)

# ---------------------------------------------------------------------------
# Compose 2x4 grid with patchwork + shared legend via cowplot
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(cowplot))

grid_plot <- wrap_plots(panels, nrow = 2, ncol = 4) +
  plot_annotation(
    title = "H5N1 Phylogeny — All Segments",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                margin = margin(b = 6))
    )
  )

final_plot <- plot_grid(
  grid_plot,
  shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.06)
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
dir.create(dirname(output_png), showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename  = output_png,
  plot      = final_plot,
  width     = 24,
  height    = 14,
  dpi       = 300,
  bg        = "white",
  limitsize = FALSE
)
cat("Saved 2x4 panel plot to", output_png, "\n")
