#!/usr/bin/env Rscript
# plot_8_segments.R
# Plots 8 segment trees as individual panels arranged in 2 figures (4 segments each).
# Each segment has:
# - Top row: rooted tree, background clades collapsed as triangles,
#   flu samples coloured by role, context samples in grey, with a box around the flu cluster.
# - Bottom row: zoomed-in subtree of the flu cluster, tips coloured by country.
#
# Usage:
#   Rscript plot_8_segments.R <meta.csv> <outgroup> <out.png> <tree1> ... <tree8>
#   (out.png will be used to generate out_1.png and out_2.png)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(stringr)
  library(patchwork)
  library(cowplot)
})

source("code/segment_analysis/tree_aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11) {
  stop("Usage: Rscript plot_8_segments.R <meta.csv> <outgroup> <out.png> <tree1> ... <tree8>")
}

metadata_path  <- args[1]
outgroup_sample <- args[2]
output_png     <- args[3]
tree_paths     <- args[4:11]

# Derive segment names from directory component
segments <- sapply(tree_paths, function(x) {
  parts <- unlist(strsplit(x, "/"))
  parts[length(parts) - 2]
})

meta          <- read.csv(metadata_path, stringsAsFactors = FALSE)
flu_tips_all  <- meta$file_name[grepl("^flu_", meta$expected_role)]
role_lookup   <- setNames(meta$expected_role, meta$file_name)

# Get top N countries for consistent coloring
top_countries <- meta %>% count(country, sort=TRUE) %>% head(10) %>% pull(country)
country_palette <- scales::hue_pal()(length(top_countries))
names(country_palette) <- top_countries

build_segment_panels <- function(tree_path, segment_name, flu_tips_all, role_lookup, outgroup_sample) {
  tree <- read.tree(tree_path)
  
  if (!is.na(outgroup_sample) && outgroup_sample %in% tree$tip.label) {
    tree <- ape::root(tree, outgroup = outgroup_sample, resolve.root = TRUE)
  } else {
    tree <- phytools::midpoint.root(tree)
  }
  tree <- ape::ladderize(tree, right = FALSE)

  flu_tips_in_tree <- intersect(tree$tip.label, flu_tips_all)
  n_tips <- Ntip(tree)

  d <- ggtree::fortify(tree)

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

  if (length(nodes_to_drop) > 0) {
    d <- d[!d$node %in% nodes_to_drop, ]
  }

  tip_data <- d[d$isTip, ] |>
    left_join(meta, by = c("label" = "file_name")) |>
    mutate(
      is_flu = label %in% flu_tips_in_tree,
      role   = ifelse(is_flu, expected_role, "context")
    )

  p_top <- ggplot() +
    geom_tree(data = d, color = "grey35", linewidth = 0.25) +
    theme_tree2() +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 4)),
      axis.line.x     = element_line(color = "grey40", linewidth = 0.35),
      axis.ticks.x    = element_line(color = "grey40", linewidth = 0.3),
      axis.text.x     = element_text(size = 6, color = "grey30"),
      plot.margin     = margin(6, 8, 6, 4),
      legend.position = "none"
    ) +
    ggtitle(segment_name)

  if (length(poly_list) > 0) {
    poly_df <- bind_rows(poly_list)
    p_top <- p_top + geom_polygon(
      data = poly_df, aes(x = x, y = y, group = group),
      fill = "grey82", color = "grey55", alpha = 0.75, linewidth = 0.25
    )
  }

  ctx <- tip_data[!tip_data$is_flu, ]
  if (nrow(ctx) > 0) {
    p_top <- p_top + geom_point(
      data = ctx, aes(x = x, y = y),
      color = "grey70", shape = 16, size = 0.8, alpha = 0.7
    )
  }

  flu_df <- tip_data[tip_data$is_flu, ]
  if (nrow(flu_df) > 0) {
    flu_df$role <- normalize_role_vector(setNames(flu_df$expected_role, flu_df$label))
    p_top <- p_top + geom_point(
      data = flu_df, aes(x = x, y = y, color = role),
      shape = 17, size = 2.2, alpha = 0.95
    ) + scale_color_manual(values = panel_type_colors, na.value = "grey70")
  }

  # Draw rectangle and extract subtree
  p_bot <- NULL
  if (length(flu_tips_in_tree) > 0) {
    if (length(flu_tips_in_tree) > 1) {
      mrca_node <- ape::getMRCA(tree, flu_tips_in_tree)
      desc_mrca <- phangorn::Descendants(tree, mrca_node, type = "all")
      mrca_nodes_all <- c(mrca_node, desc_mrca)
    } else {
      mrca_node <- which(tree$tip.label == flu_tips_in_tree[1])
      mrca_nodes_all <- c(mrca_node)
    }

    mrca_data <- d[d$node %in% mrca_nodes_all, ]
    if (nrow(mrca_data) > 0) {
      xmin <- min(mrca_data$x)
      xmax <- max(mrca_data$x)
      ymin <- min(mrca_data$y)
      ymax <- max(mrca_data$y)
      x_range <- max(d$x, na.rm=TRUE) - min(d$x, na.rm=TRUE)
      
      p_top <- p_top + geom_rect(
        aes(xmin = xmin - 0.02 * x_range, xmax = xmax + 0.03 * x_range,
            ymin = ymin - 0.6, ymax = ymax + 0.6),
        fill = NA, color = "red", linetype = "dashed", linewidth = 0.5
      )
    }

    # Extract and plot subtree
    if (length(flu_tips_in_tree) > 1) {
      subtree <- ape::extract.clade(tree, mrca_node)
    } else {
      subtree <- ape::keep.tip(tree, flu_tips_in_tree[1])
    }
    
    d_sub <- ggtree::fortify(subtree)
    tip_data_sub <- d_sub[d_sub$isTip, ] |>
      left_join(meta, by = c("label" = "file_name")) |>
      mutate(
        is_flu = label %in% flu_tips_in_tree,
        role   = ifelse(is_flu, expected_role, "context")
      )
    
    # ensure country colors are mostly covered
    tip_data_sub$country_plot <- ifelse(tip_data_sub$country %in% names(country_palette), tip_data_sub$country, "Other")

    # Split into context and flu
    ctx_sub <- tip_data_sub[!tip_data_sub$is_flu, ]
    flu_sub <- tip_data_sub[tip_data_sub$is_flu, ]

    p_bot <- ggplot() +
      geom_tree(data = d_sub, color = "grey35", linewidth = 0.3) +
      theme_tree2() +
      theme(
        axis.line.x     = element_line(color = "grey40", linewidth = 0.35),
        axis.ticks.x    = element_line(color = "grey40", linewidth = 0.3),
        axis.text.x     = element_text(size = 6, color = "grey30"),
        plot.margin     = margin(6, 8, 6, 4),
        legend.position = "none"
      )

    if (nrow(ctx_sub) > 0) {
      p_bot <- p_bot + geom_point(
        data = ctx_sub,
        aes(x = x, y = y, fill = country_plot),
        shape = 21, size = 2.5, color = "white", stroke = 0.2, alpha = 0.5
      )
    }

    if (nrow(flu_sub) > 0) {
      flu_sub$role <- normalize_role_vector(setNames(flu_sub$expected_role, flu_sub$label))
      p_bot <- p_bot + geom_point(
        data = flu_sub,
        aes(x = x, y = y, color = role),
        shape = 17, size = 6, alpha = 1.0
      ) + scale_color_manual(values = panel_type_colors, na.value = "grey70")
    }
    
    # We must add scale_fill_manual for country colors
    p_bot <- p_bot + scale_fill_manual(values = c(country_palette, "Other" = "grey70"))
  }

  list(top = p_top, bot = p_bot)
}

panels <- vector("list", length(tree_paths))
for (i in seq_along(tree_paths)) {
  cat("Processing", segments[i], "...\n")
  panels[[i]] <- build_segment_panels(
    tree_path       = tree_paths[i],
    segment_name    = segments[i],
    flu_tips_all    = flu_tips_all,
    role_lookup     = role_lookup,
    outgroup_sample = outgroup_sample
  )
}

# Build Legends
legend_roles <- names(panel_type_colors)
legend_shapes <- c(flu_costa = 17, flu_andine = 17, flu_amazonia = 17, american_anchor = 16, regional_context = 16)
legend_data_role <- data.frame(x=seq_along(legend_roles), y=1, role=legend_roles, stringsAsFactors=FALSE)

p_leg_role <- ggplot(legend_data_role, aes(x=x, y=y, color=role, shape=role)) +
  geom_point(size=3) +
  scale_color_manual(values=panel_type_colors, labels=panel_type_labels, name="Role") +
  scale_shape_manual(values=legend_shapes, labels=panel_type_labels, name="Role") +
  theme_void() + theme(legend.position="bottom", legend.box="horizontal")
leg_role <- cowplot::get_legend(p_leg_role)

country_levels <- c(names(country_palette), "Other")
country_colors <- c(country_palette, "Other" = "grey70")
p_leg_country <- ggplot(data.frame(country=country_levels), aes(x=1,y=1,fill=country)) +
  geom_point(shape=21, size=3) +
  scale_fill_manual(values=country_colors, name="Country") +
  theme_void() + theme(legend.position="bottom", legend.box="horizontal")
leg_country <- cowplot::get_legend(p_leg_country)

combined_legends <- plot_grid(leg_role, leg_country, ncol=1)

# Compose Graphic 1 (first 4 segments)
layout_1 <- plot_grid(
  plot_grid(panels[[1]]$top, panels[[2]]$top, panels[[3]]$top, panels[[4]]$top, ncol=4),
  plot_grid(panels[[1]]$bot, panels[[2]]$bot, panels[[3]]$bot, panels[[4]]$bot, ncol=4),
  ncol = 1, rel_heights = c(1, 1)
)
final_plot_1 <- plot_grid(layout_1, combined_legends, ncol=1, rel_heights=c(1, 0.1))

# Compose Graphic 2 (last 4 segments)
layout_2 <- plot_grid(
  plot_grid(panels[[5]]$top, panels[[6]]$top, panels[[7]]$top, panels[[8]]$top, ncol=4),
  plot_grid(panels[[5]]$bot, panels[[6]]$bot, panels[[7]]$bot, panels[[8]]$bot, ncol=4),
  ncol = 1, rel_heights = c(1, 1)
)
final_plot_2 <- plot_grid(layout_2, combined_legends, ncol=1, rel_heights=c(1, 0.1))

out1 <- sub(".png$", "_1.png", output_png)
out2 <- sub(".png$", "_2.png", output_png)

dir.create(dirname(output_png), showWarnings = FALSE, recursive = TRUE)

ggsave(filename=out1, plot=final_plot_1, width=24, height=14, dpi=300, bg="white", limitsize=FALSE)
cat("Saved plot 1 to", out1, "\n")

ggsave(filename=out2, plot=final_plot_2, width=24, height=14, dpi=300, bg="white", limitsize=FALSE)
cat("Saved plot 2 to", out2, "\n")

# Create a dummy file for Snakemake's expected output
file.create(output_png)
cat("Created dummy file", output_png, "for Snakemake compatibility\n")
