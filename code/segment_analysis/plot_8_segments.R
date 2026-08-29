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

# Load genoflu details
genoflu_path <- "metadata/genoflu_results.tsv"
if (file.exists(genoflu_path)) {
  gf <- read.delim(genoflu_path, stringsAsFactors = FALSE, check.names = FALSE)
  gf_cols <- c("Strain", "Genotype List Used, >=98.0%", "Genotype Sample Title List", "Genotype Percent Match List")
  gf_cols_present <- intersect(gf_cols, names(gf))
  gf_sub <- gf[, gf_cols_present, drop = FALSE]
  names(gf_sub)[names(gf_sub) == "Genotype List Used, >=98.0%"] <- "Genotype_List"
  names(gf_sub)[names(gf_sub) == "Genotype Sample Title List"] <- "Genotype_Top_Hits"
  names(gf_sub)[names(gf_sub) == "Genotype Percent Match List"] <- "Genotype_Percent_Match"
  meta <- merge(meta, gf_sub, by.x = "file_name", by.y = "Strain", all.x = TRUE)
} else {
  meta$Genotype_List <- NA
  meta$Genotype_Top_Hits <- NA
  meta$Genotype_Percent_Match <- NA
}

# Define outliers per segment based on user input
outliers_map <- list(
  "EPI_ISL_18137671" = c("NS")
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
  
  # Identify outliers for this segment
  outliers_this_seg <- names(outliers_map)[sapply(outliers_map, function(x) segment_name %in% x)]
  
  flu_tips_in_tree <- intersect(tree$tip.label, flu_tips_all)
  flu_anchors_in_tree <- setdiff(flu_tips_in_tree, outliers_this_seg)
  
  n_tips <- Ntip(tree)

  d <- ggtree::fortify(tree)


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
  scale_len_top <- 0.005 # Fixed standard scale for top plot
  scale_x_start_top <- x_max - scale_len_top - (x_max * 0.05) # Shifted a bit left
  p_top <- p_top +
    annotate("segment", x = scale_x_start_top, xend = scale_x_start_top + scale_len_top, y = scale_y, yend = scale_y, linewidth = 1.2, color = "grey30") +
    annotate("text", x = scale_x_start_top + (scale_len_top / 2), y = scale_y - (y_max * 0.035), label = paste(scale_len_top, "subs/site"), size = 6, color = "grey30")






  tip_data$role <- normalize_role_vector(setNames(tip_data$expected_role, tip_data$label))
  tip_data$pt_size  <- ifelse(tip_data$is_flu, 2.2, 0.8)
  tip_data$pt_alpha <- ifelse(tip_data$is_flu, 1.0, 
                              ifelse(tip_data$role == "american_anchor", 0.3, 0.7))

  if (nrow(tip_data) > 0) {
    p_top <- p_top + geom_point(
      data = tip_data, aes(x = x, y = y, color = role, shape = host_type,
                           size = pt_size, alpha = pt_alpha),
      stroke = 0.2
    ) + scale_color_manual(values = panel_type_colors, na.value = "grey70") +
      scale_shape_manual(values = c("domesticated bird"=15, "wild bird"=16, "domesticated mammal"=17, "wild mammal"=18, "?"=3)) +
      scale_size_identity() + scale_alpha_identity()
      
    booby_data_top <- tip_data[tip_data$label == "EPI_ISL_18137671", ]
    if (nrow(booby_data_top) > 0) {
      p_top <- p_top + geom_text(data = booby_data_top, aes(x = x, y = y), label = "*", color = "red", size = 12, fontface = "bold", hjust = -0.2, vjust = 0.75)
    }
    p_top <- p_top + coord_cartesian(clip = "off")
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
    d_sub_core <- d_sub[d_sub$isTip & !(d_sub$label %in% outliers_in_sub), ]
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
    tip_data_sub$pt_alpha <- ifelse(tip_data_sub$is_flu, 1.0, 
                                    ifelse(tip_data_sub$role == "american_anchor", 0.3, 0.6))

    p_bot <- ggplot() +
      geom_tree(data = d_sub, color = "grey35", linewidth = 0.3) +
      theme_tree() +
      theme(
        plot.margin     = margin(10, 10, 20, 10),
        legend.position = "none",
        panel.border    = element_rect(color = "red", linetype = "dashed", fill = NA, linewidth = 0.5)
      ) +
      coord_cartesian(xlim = c(x_min_core, x_max_core + (x_range_core * 0.05)), clip = "off")

    scale_len_bot <- 0.003 # Fixed scale bar for all zoomed segments
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
        
      booby_data_bot <- tip_data_sub[tip_data_sub$label == "EPI_ISL_18137671", ]
      if (nrow(booby_data_bot) > 0) {
        p_bot <- p_bot + geom_text(data = booby_data_bot, aes(x = x, y = y), label = "*", color = "red", size = 12, fontface = "bold", hjust = -0.2, vjust = 0.75)
      }
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
    legend.justification = "left",
    legend.text = element_text(margin = margin(r = 15), size=18),
    legend.title = element_text(size=20, face="bold", vjust=0.5, margin=margin(r=10))
  ) +
  guides(color = guide_legend(override.aes = list(size = 9, alpha = c(1, 1, 1, 0.3, 0.4)), nrow = 1, title.position = "left", title.vjust = 0.5, keywidth = unit(2, "cm")))
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
    legend.margin = margin(t=0, b=0, r=0, l=15),
    legend.justification = "left",
    legend.text = element_text(margin = margin(r = 40), size=18),
    legend.title = element_text(size=20, face="bold", vjust=0.5, margin=margin(r=10))
  ) +
  guides(shape = guide_legend(override.aes = list(size = 9), nrow = 1, title.position = "left", title.vjust = 0.5, keywidth = unit(2, "cm")))
leg_shape <- cowplot::get_legend(p_leg_shape)

p_leg_reassortant <- ggplot(data.frame(x=1, y=1, type="NS outlier"), aes(x,y, color=type)) + 
  geom_point(shape="*", size=16, stroke=1.5) +
  scale_color_manual(values=c("NS outlier"="red"), name=NULL) +
  theme_void() + theme(
    legend.position="bottom", 
    legend.margin = margin(t=0, b=0, r=0, l=0),
    legend.justification = "left",
    legend.text = element_text(size=18, margin=margin(r=10))
  ) + guides(color = guide_legend(override.aes = list(size = 14), nrow = 1))
leg_reassortant <- cowplot::get_legend(p_leg_reassortant)

grid_2x2 <- plot_grid(
  leg_reassortant, leg_role,
  leg_shape, NULL,
  nrow = 2,
  rel_widths = c(0.25, 1.5),
  rel_heights = c(1, 1)
)
combined_legends <- plot_grid(NULL, grid_2x2, nrow=1, rel_widths=c(0.12, 1))

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


# -------------------------------------------------------------
# Export ignored_in_per_segment.csv based on outliers_map
# -------------------------------------------------------------
out_csv <- "results/ignored_in_per_segment.csv"

b32_constellation <- list(
  "PB2" = "am2.1",
  "PB1" = "am1.2",
  "PA"  = "ea1",
  "HA"  = "ea1",
  "NP"  = "am1.4.1",
  "NA"  = "ea1",
  "MP"  = "ea1",
  "NS"  = "am1.1"
)

# Extract information from meta to match the original CSV structure
ignored_list <- lapply(names(outliers_map), function(sid) {
  row <- meta[meta$file_name == sid, ]
  if (nrow(row) == 0) return(NULL)
  
  segs <- paste(sort(outliers_map[[sid]]), collapse=", ")
  geno_raw <- as.character(row$genotype[1])
  if (is.na(geno_raw)) geno_raw <- ""
  
  geno_list <- as.character(row$Genotype_List[1])
  if (is.na(geno_list)) geno_list <- ""
  
  incompatible_segs <- c()
  
  if (grepl("B3.2", geno_raw) && !grepl("Not assigned", geno_raw)) {
    comp <- "Yes"
  } else {
    comp <- "No"
  }
  
  top_hits <- as.character(row$Genotype_Top_Hits[1])
  if (is.na(top_hits)) top_hits <- ""
  pcts <- as.character(row$Genotype_Percent_Match[1])
  if (is.na(pcts)) pcts <- ""
  
  ignored_segs_list <- outliers_map[[sid]]
  
  ignored_info <- c()
  other_info <- c()
  
  if (top_hits != "") {
    hits_vec <- strsplit(top_hits, ",\\s*")[[1]]
    pcts_vec <- strsplit(pcts, ",\\s*")[[1]]
    
    for (k in seq_along(hits_vec)) {
      h <- hits_vec[k]
      p <- if (k <= length(pcts_vec)) pcts_vec[k] else ""
      parts <- strsplit(h, ":")[[1]]
      lin <- parts[1]
      seg <- parts[length(parts)]
      
      info_str <- sprintf("%s (%s, %s)", seg, lin, p)
      if (seg %in% ignored_segs_list) {
        ignored_info <- c(ignored_info, info_str)
      } else {
        other_info <- c(other_info, info_str)
      }
    }
  }
  
  data.frame(
    Sample_ID = sid,
    Segments_Ignored_In_Plots = segs,
    ">98% match to B3.2 genotype" = comp,
    Host = as.character(row$host_type[1]),
    Collection_Date = as.character(row$collection_date[1]),
    Location = paste0(row$country[1], " (", row$province[1], ")"),
    Role = as.character(row$expected_role[1]),
    Ignored_Segments_Top_Hits = paste(ignored_info, collapse="; "),
    Other_Segments_Top_Hits = paste(other_info, collapse="; "),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})


ignored_df <- do.call(rbind, ignored_list)

# Additional Ecuadorian sequences retained in trees with sub-threshold (<98%) GenoFLU hits
ec_retained <- meta[grepl("^flu_", meta$expected_role) & !(meta$file_name %in% names(outliers_map)), ]
retained_divergent_list <- list()

if (nrow(ec_retained) > 0) {
  for (i in seq_len(nrow(ec_retained))) {
    row_ec <- ec_retained[i, ]
    top_hits <- as.character(row_ec$Genotype_Top_Hits)
    pcts <- as.character(row_ec$Genotype_Percent_Match)
    if (!is.na(top_hits) && top_hits != "") {
      hits_vec <- strsplit(top_hits, ",\\s*")[[1]]
      pcts_vec <- strsplit(pcts, ",\\s*")[[1]]
      
      for (j in seq_along(hits_vec)) {
        h <- hits_vec[j]
        p_str <- if (j <= length(pcts_vec)) pcts_vec[j] else ""
        p_val <- as.numeric(gsub("%", "", p_str))
        parts <- strsplit(h, ":")[[1]]
        lin <- parts[1]
        seg <- parts[length(parts)]
        canon <- if (seg %in% names(b32_constellation)) b32_constellation[[seg]] else "unknown"
        
        # Only include if the top-hit lineage diverges from the canonical B3.2 constellation
        if (seg %in% names(b32_constellation) && lin != canon) {
          retained_divergent_list[[length(retained_divergent_list) + 1]] <- data.frame(
            Sample_ID = row_ec$file_name,
            Segment = seg,
            GenoFLU_Top_Hit_Lineage = lin,
            Canonical_B3.2_Lineage = canon,
            Percent_Identity = p_str,
            Full_Top_Hit_Reference = h,
            Host = row_ec$host_type,
            Location = paste0(row_ec$country, " (", row_ec$province, ")"),
            Collection_Date = row_ec$collection_date,
            Role = row_ec$expected_role,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
retained_divergent_df <- if (length(retained_divergent_list) > 0) do.call(rbind, retained_divergent_list) else NULL

if (!is.null(ignored_df)) {
  write.csv(ignored_df, out_csv, row.names = FALSE, quote = TRUE)
  
  if (!is.null(retained_divergent_df) && nrow(retained_divergent_df) > 0) {
    cat("\n", file=out_csv, append=TRUE)
    cat("\"--- Additional Ecuadorian Sequences Retained in Trees with Non-Canonical B3.2 GenoFLU Hits ---\"\n", file=out_csv, append=TRUE)
    write.table(retained_divergent_df, out_csv, row.names = FALSE, quote = TRUE, sep = ",", append = TRUE)
  }
  
  # Append the official B3.2 constellation at the end of the CSV
  cat("\n", file=out_csv, append=TRUE)
  cat("\"Official B3.2 Constellation\",\"Lineage\"\n", file=out_csv, append=TRUE)
  for (seg in names(b32_constellation)) {
    cat(sprintf("\"%s\",\"%s\"\n", seg, b32_constellation[[seg]]), file=out_csv, append=TRUE)
  }
}
# -------------------------------------------------------------

