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

# Define outliers per segment based on user input
outliers_map <- list(
  "EPI_ISL_20450124" = c("MP", "NA", "NP", "NS"),
  "EPI_ISL_20450123" = c("MP"),
  "EPI_ISL_20450128" = c("MP", "NS"),
  "EPI_ISL_20450130" = c("MP"),
  "EPI_ISL_20450131" = c("MP"),
  "EPI_ISL_20450132" = c("MP"),
  "EPI_ISL_20450133" = c("MP"),
  "EPI_ISL_18137671" = c("NS", "NA")
)

role_lookup   <- setNames(meta$expected_role, meta$file_name)

# Get top N countries for consistent coloring
top_countries <- meta %>% count(country, sort=TRUE) %>% head(10) %>% pull(country)
country_palette <- scales::hue_pal()(length(top_countries))
names(country_palette) <- top_countries

build_segment_panels <- function(tree_path, segment_name, flu_tips_all, outliers_map, role_lookup, outgroup_sample) {
  tree <- read.tree(tree_path)
  
  if (!is.na(outgroup_sample) && outgroup_sample %in% tree$tip.label) {
    tree <- ape::root(tree, outgroup = outgroup_sample, resolve.root = TRUE)
  } else {
    tree <- phytools::midpoint.root(tree)
  }
  tree <- ape::ladderize(tree, right = FALSE)

  flu_tips_in_tree <- intersect(tree$tip.label, flu_tips_all)
  
  # Determine which tips to use as anchors for this segment's MRCA
  outliers_this_seg <- names(outliers_map)[sapply(outliers_map, function(x) segment_name %in% x)]
  flu_anchors_in_tree <- setdiff(flu_tips_in_tree, outliers_this_seg)
  
  n_tips <- Ntip(tree)

  d <- ggtree::fortify(tree)

  poly_list  <- list()
  poly_id    <- 1
  nodes_to_drop <- c()

  if (length(flu_tips_in_tree) > 0) {
    desc_tips <- phangorn::Descendants(tree, type = "tips")
    desc_all  <- phangorn::Descendants(tree, type = "all")

    for (node in (n_tips + 1):(n_tips + Nnode(tree))) {
      if (node %in% nodes_to_drop) next
      tips_in_node <- tree$tip.label[desc_tips[[node]]]
      if (length(intersect(tips_in_node, flu_tips_in_tree)) == 0) {
        roles_node <- meta$expected_role[match(tips_in_node, meta$file_name)]
        roles_node <- roles_node[!is.na(roles_node) & roles_node != ""]
        if (length(unique(roles_node)) > 1) next
        
        desc_nodes <- desc_all[[node]]
        x_root     <- d$x[d$node == node]
        y_root     <- d$y[d$node == node]
        desc_data  <- d[d$node %in% desc_nodes, ]
        if (nrow(desc_data) > 0) {
          majority_role <- if (length(roles_node) > 0) unique(roles_node)[1] else "unknown"
          majority_role <- normalize_ecuador_role(majority_role)

          poly_list[[poly_id]] <- data.frame(
            x     = c(x_root, max(desc_data$x), max(desc_data$x)),
            y     = c(y_root, max(desc_data$y), min(desc_data$y)),
            group = paste0("poly_", poly_id),
            poly_role = majority_role
          )
          poly_id       <- poly_id + 1
          nodes_to_drop <- c(nodes_to_drop, desc_nodes)
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
      role   = expected_role
    )

  p_top <- ggplot() +
    geom_tree(data = d, color = "grey35", linewidth = 0.25) +
    theme_tree() +
    theme(
      plot.title      = element_text(face = "bold", size = 32, hjust = 0.5, margin = margin(b = 10)),
      plot.margin     = margin(6, 8, 20, 4),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off") +
    ggtitle(segment_name)

  # Add top scale dynamically (bottom right)
  x_max <- max(d$x, na.rm=TRUE)
  y_max <- max(d$y, na.rm=TRUE)
  scale_y <- y_max * -0.015
  scale_len_top <- signif(x_max * 0.15, 1) # 15% of the total width
  scale_x_start_top <- x_max - scale_len_top - (x_max * 0.05) # Shifted a bit left
  p_top <- p_top +
    annotate("segment", x = scale_x_start_top, xend = scale_x_start_top + scale_len_top, y = scale_y, yend = scale_y, linewidth = 1.2, color = "grey30") +
    annotate("text", x = scale_x_start_top + (scale_len_top / 2), y = scale_y - (y_max * 0.035), label = paste(scale_len_top, "subs/site"), size = 6, color = "grey30")





  if (length(poly_list) > 0) {
    poly_df <- bind_rows(poly_list)
    p_top <- p_top + geom_polygon(
      data = poly_df, aes(x = x, y = y, group = group, fill = poly_role, color = poly_role),
      alpha = 0.25, linewidth = 0.25
    ) + scale_fill_manual(values = panel_type_colors, na.value = "grey82")
  }

  tip_data$role <- normalize_role_vector(setNames(tip_data$expected_role, tip_data$label))
  tip_data$pt_size  <- ifelse(tip_data$is_flu, 2.2, 0.8)
  tip_data$pt_alpha <- ifelse(tip_data$is_flu, 0.95, 0.7)

  if (nrow(tip_data) > 0) {
    p_top <- p_top + geom_point(
      data = tip_data, aes(x = x, y = y, color = role, shape = host_type,
                           size = pt_size, alpha = pt_alpha),
      stroke = 0.2
    ) + scale_color_manual(values = panel_type_colors, na.value = "grey70") +
      scale_shape_manual(values = c("domesticated bird"=15, "wild bird"=16, "domesticated mammal"=17, "wild mammal"=18, "?"=3)) +
      scale_size_identity() + scale_alpha_identity()
  }

  # Draw rectangle and extract subtree
  p_bot <- NULL
  if (length(flu_anchors_in_tree) > 0) {
    if (length(flu_anchors_in_tree) > 1) {
      mrca_node <- ape::getMRCA(tree, flu_anchors_in_tree)
      desc_mrca <- phangorn::Descendants(tree, mrca_node, type = "all")
      mrca_nodes_all <- c(mrca_node, desc_mrca)
    } else {
      mrca_node <- which(tree$tip.label == flu_anchors_in_tree[1])
      mrca_nodes_all <- c(mrca_node)
    }

    mrca_data <- d[d$node %in% mrca_nodes_all, ]
    if (nrow(mrca_data) > 0) {
      mrca_tips_core <- mrca_data[mrca_data$isTip & !(mrca_data$label %in% outliers_this_seg), ]
      xmin <- d$x[d$node == mrca_node]
      xmax <- max(mrca_tips_core$x, na.rm = TRUE)
      ymin <- min(mrca_data$y, na.rm=TRUE)
      ymax <- max(mrca_data$y, na.rm=TRUE)
      x_range <- max(d$x, na.rm=TRUE) - min(d$x, na.rm=TRUE)
      
      p_top <- p_top + geom_rect(
        aes(xmin = xmin - 0.02 * x_range, xmax = xmax + 0.03 * x_range,
            ymin = ymin - 0.6, ymax = ymax + 0.6),
        fill = NA, color = "red", linetype = "dashed", linewidth = 2.5
      )
    }

    # Extract and plot subtree
    if (length(flu_tips_in_tree) > 1) {
      subtree <- ape::extract.clade(tree, mrca_node)
    } else {
      subtree <- ape::keep.tip(tree, flu_tips_in_tree[1])
    }
    # Drop extreme outliers from the zoomed subtree so they do not distort the scale
    # INSTEAD OF DROPPING: We will keep them in the tree but restrict the plot limits
    d_sub <- ggtree::fortify(subtree)
    
    outliers_in_sub <- intersect(outliers_this_seg, subtree$tip.label)
    d_sub_core <- d_sub[!(d_sub$label %in% outliers_in_sub), ]
    x_max_core <- max(d_sub_core$x, na.rm=TRUE)
    x_min_core <- min(d_sub_core$x, na.rm=TRUE)
    x_range_core <- x_max_core - x_min_core

    tip_data_sub <- d_sub[d_sub$isTip, ] |>
      left_join(meta, by = c("label" = "file_name")) |>
      mutate(
        is_flu = label %in% flu_tips_in_tree,
        role   = expected_role
      )
    
    tip_data_sub$role <- normalize_role_vector(setNames(tip_data_sub$expected_role, tip_data_sub$label))
    tip_data_sub$pt_size  <- ifelse(tip_data_sub$is_flu, 2, 2.5)
    tip_data_sub$pt_alpha <- ifelse(tip_data_sub$is_flu, 1.0, 0.6)

    p_bot <- ggplot() +
      geom_tree(data = d_sub, color = "grey35", linewidth = 0.3) +
      theme_tree() +
      theme(
        plot.margin     = margin(10, 10, 20, 10),
        legend.position = "none",
        panel.border    = element_rect(color = "red", linetype = "dashed", fill = NA, linewidth = 0.5)
      ) +
      coord_cartesian(xlim = c(x_min_core, x_max_core + (x_range_core * 0.05)), clip = "on")

    scale_len_bot <- signif(x_range_core * 0.2, 1) # ~20% of the core subtree width
    scale_x_start <- x_max_core - scale_len_bot - (x_range_core * 0.05) # Shift left a bit
    y_max_sub <- max(d_sub$y, na.rm=TRUE)
    
    p_bot <- p_bot +
      annotate("segment", x = scale_x_start, xend = scale_x_start + scale_len_bot, y = -(y_max_sub * 0.03), yend = -(y_max_sub * 0.03), linewidth = 1, color = "grey30") +
      annotate("text", x = scale_x_start + (scale_len_bot / 2), y = -(y_max_sub * 0.06), label = paste(scale_len_bot, "subs/site"), size = 5, color = "grey30")

    if (nrow(tip_data_sub) > 0) {
      p_bot <- p_bot + geom_point(
        data = tip_data_sub,
        aes(x = x, y = y, color = role, shape = host_type,
            size = pt_size, alpha = pt_alpha),
        stroke = 0.2
      ) + scale_color_manual(values = panel_type_colors, na.value = "grey70") +
        scale_shape_manual(values = c("domesticated bird"=15, "wild bird"=16, "domesticated mammal"=17, "wild mammal"=18, "?"=3)) +
        scale_size_identity() + scale_alpha_identity()
    }
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
    outliers_map    = outliers_map,
    role_lookup     = role_lookup,
    outgroup_sample = outgroup_sample
  )
}

# Build Legends
legend_roles <- PANEL_TYPE_RIBBON_ORDER
legend_data_role <- data.frame(x=seq_along(legend_roles), y=1, role=legend_roles, stringsAsFactors=FALSE)

p_leg_role <- ggplot(legend_data_role, aes(x=x, y=y, color=role)) +
  geom_point(size=12, shape=15) +
  scale_color_manual(values=panel_type_colors, labels=panel_type_labels, name="Geographic groups") +
  theme_void() + theme(
    legend.position="bottom", 
    legend.box="horizontal",
    legend.margin = margin(t=0, b=0, r=0, l=0),
    legend.text = element_text(margin = margin(r = 15), size=18),
    legend.title = element_text(size=20, face="bold", vjust=0.5, margin=margin(r=20))
  ) +
  guides(color = guide_legend(override.aes = list(size = 9, alpha = c(1, 1, 1, 0.4, 0.4)), nrow = 1, title.position = "left", title.vjust = 0.5, keywidth = unit(2, "cm")))
leg_role <- cowplot::get_legend(p_leg_role)

host_shape_mapping <- c(
  "domesticated bird" = 15,
  "wild bird" = 16,
  "domesticated mammal" = 17,
  "wild mammal" = 18,
  "Unknown" = 3
)

p_leg_shape <- ggplot(data.frame(host=names(host_shape_mapping)), aes(x=1,y=1,shape=host)) +
  geom_point(size=12, color="grey30") +
  scale_shape_manual(values=host_shape_mapping, name="Host Type") +
  theme_void() + theme(
    legend.position="bottom", 
    legend.box="horizontal",
    legend.margin = margin(t=0, b=0, r=0, l=0),
    legend.text = element_text(margin = margin(r = 40), size=18),
    legend.title = element_text(size=20, face="bold", vjust=0.5, margin=margin(r=20))
  ) +
  guides(shape = guide_legend(override.aes = list(size = 9), nrow = 1, title.position = "left", title.vjust = 0.5, keywidth = unit(2, "cm")))
leg_shape <- cowplot::get_legend(p_leg_shape)

combined_legends <- plot_grid(leg_role, leg_shape, ncol=1, rel_heights=c(1, 1))

# Reorder trees to: PB1, PB2, NP, NS and PA, HA, NA, MP
desired_order <- c("PB1", "PB2", "NP", "NS", "PA", "HA", "NA", "MP")
match_idx <- match(desired_order, segments)
panels_reordered <- panels[match_idx]

# Compose Graphic 1 (first 4 segments: PB1, PB2, NP, NS)
layout_1 <- plot_grid(
  plot_grid(panels_reordered[[1]]$top, panels_reordered[[2]]$top, panels_reordered[[3]]$top, panels_reordered[[4]]$top, ncol=4),
  plot_grid(panels_reordered[[1]]$bot, panels_reordered[[2]]$bot, panels_reordered[[3]]$bot, panels_reordered[[4]]$bot, ncol=4),
  ncol = 1, rel_heights = c(1, 1)
)
final_plot_1 <- plot_grid(layout_1, combined_legends, ncol=1, rel_heights=c(1, 0.12), labels=c("A", ""), label_size=60, label_fontfamily="sans")

# Compose Graphic 2 (last 4 segments: PA, HA, NA, MP)
layout_2 <- plot_grid(
  plot_grid(panels_reordered[[5]]$top, panels_reordered[[6]]$top, panels_reordered[[7]]$top, panels_reordered[[8]]$top, ncol=4),
  plot_grid(panels_reordered[[5]]$bot, panels_reordered[[6]]$bot, panels_reordered[[7]]$bot, panels_reordered[[8]]$bot, ncol=4),
  ncol = 1, rel_heights = c(1, 1)
)
final_plot_2 <- plot_grid(layout_2, combined_legends, ncol=1, rel_heights=c(1, 0.12), labels=c("B", ""), label_size=60, label_fontfamily="sans")


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
