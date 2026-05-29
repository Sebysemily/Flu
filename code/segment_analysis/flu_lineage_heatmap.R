#!/usr/bin/env Rscript
# flu_lineage_heatmap.R
# GenoFLU constellation heatmap (Ecuador panel).
#
# Usage:
#   Rscript code/segment_analysis/flu_lineage_heatmap.R <metadata_csv> [output_png]
#
# Expects metadata/H5N1_context.csv (Ecuador rows with expected_role + segment lineages).
# When sourced, exposes build_flu_lineage_heatmap() only (no side effects).

source("code/segment_analysis/lineage_palette.R")

ECUADOR_PANEL_ROLES <- c("flu_costa", "flu_andine", "flu_sierra", "flu_amazonia")

compose_heatmap_with_legend <- function(p, costa_star_color) {
  legend_panel <- cowplot::get_legend(p)
  coastal_note <- cowplot::ggdraw() +
    cowplot::draw_label(
      label = "*",
      x = 0.22,
      y = 0.90,
      hjust = 0.5,
      vjust = 0.5,
      color = costa_star_color,
      fontface = "bold",
      size = 20
    ) +
    cowplot::draw_label(
      label = "coastal",
      x = 0.22,
      y = 0.60,
      hjust = 0.5,
      vjust = 0.5,
      color = "black",
      size = 6
    )

  legend_column <- cowplot::plot_grid(
    legend_panel,
    coastal_note,
    ncol = 1,
    rel_heights = c(1, 0.10)
  )

  cowplot::plot_grid(
    p + ggplot2::theme(legend.position = "none"),
    legend_column,
    ncol = 2,
    rel_widths = c(1, 0.20),
    align = "h"
  )
}

build_flu_lineage_heatmap <- function(
    metadata_path,
    sample_order = NULL,
    title = "Avian Influenza Genotype Constellation"
) {
  flu_data <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!"expected_role" %in% colnames(flu_data)) {
    stop("metadata must contain expected_role (use metadata/H5N1_context.csv)")
  }

  flu_data <- flu_data[flu_data$expected_role %in% ECUADOR_PANEL_ROLES, , drop = FALSE]
  if (nrow(flu_data) == 0) {
    stop("no Ecuador panel samples (flu_costa / flu_andine / flu_sierra / flu_amazonia) in metadata")
  }

  key_col <- if ("file_name" %in% colnames(flu_data)) {
    "file_name"
  } else if ("EPI_ISL" %in% colnames(flu_data)) {
    "EPI_ISL"
  } else if (ncol(flu_data) >= 2) {
    colnames(flu_data)[2]
  } else {
    colnames(flu_data)[1]
  }

  segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
  unknown_label <- "unknown"
  costa_marker_fill <- ".flu_costa_panel"
  costa_star_color <- "#FF0000"
  legend_title <- "Sub-lineage"
  segment_levels <- c(rev(segments), "Panel")

  segment_flags <- flu_data %>%
    tidyr::pivot_longer(
      cols = all_of(segments),
      names_to = "Segment",
      values_to = "Assembled"
    ) %>%
    dplyr::mutate(
      Sample_ID = .data[[key_col]],
      Assembled = toupper(stringr::str_trim(Assembled)) == "SI"
    )

  lineages <- flu_data %>%
    tidyr::pivot_longer(
      cols = ends_with("_lineage"),
      names_to = "Segment",
      values_to = "Lineage_Number"
    ) %>%
    dplyr::mutate(
      Segment = stringr::str_replace(Segment, "_lineage", ""),
      Sample_ID = .data[[key_col]],
      Lineage_Number = stringr::str_trim(Lineage_Number)
    )

  clean_metadata <- segment_flags %>%
    dplyr::filter(Assembled) %>%
    dplyr::left_join(
      lineages %>% dplyr::select(Sample_ID, Segment, Lineage_Number),
      by = c("Sample_ID", "Segment")
    ) %>%
    dplyr::mutate(
      Lineage_Number = dplyr::if_else(
        !is.na(Lineage_Number) & Lineage_Number != "",
        Lineage_Number,
        unknown_label
      )
    )

  sample_roles <- flu_data %>%
    dplyr::transmute(
      Sample_ID = .data[[key_col]],
      expected_role = .data[["expected_role"]]
    ) %>%
    dplyr::distinct()

  plotted_samples <- clean_metadata %>%
    dplyr::distinct(Sample_ID) %>%
    dplyr::left_join(
      sample_roles %>% dplyr::select(Sample_ID, expected_role),
      by = "Sample_ID"
    )

  if (!is.null(sample_order)) {
    sample_order <- intersect(as.character(sample_order), plotted_samples$Sample_ID)
    plotted_samples <- plotted_samples %>%
      dplyr::filter(Sample_ID %in% sample_order)
    clean_metadata <- clean_metadata %>%
      dplyr::filter(Sample_ID %in% sample_order)
  }

  if (is.null(sample_order) || length(sample_order) == 0) {
    sample_order <- plotted_samples %>%
      dplyr::arrange(expected_role, Sample_ID) %>%
      dplyr::pull(Sample_ID)
  }

  clean_metadata <- clean_metadata %>%
    dplyr::mutate(
      Sample_ID = factor(Sample_ID, levels = sample_order),
      Segment = factor(Segment, levels = segment_levels)
    )

  costa_panel_markers <- plotted_samples %>%
    dplyr::filter(expected_role == "flu_costa") %>%
    dplyr::transmute(
      Sample_ID = factor(Sample_ID, levels = sample_order),
      Segment = factor("Panel", levels = segment_levels),
      Lineage_Number = costa_marker_fill,
      star = "*"
    )

  unique_lineages <- unique(clean_metadata$Lineage_Number)
  palette <- build_lineage_palette(
    unique_lineages,
    unknown_label = unknown_label,
    costa_marker_fill = costa_marker_fill
  )

  # Build legend with individual sub-lineages (am1.1, am1.2, ea1, etc.)
  # ordered am* first (light→dark orange), then ea* (light→dark blue)
  ordered_lins <- order_lineages_for_legend(unique_lineages, unknown_label = unknown_label)
  ordered_lins <- ordered_lins[
    tolower(ordered_lins) != unknown_label & ordered_lins != "" &
    !grepl("^\\.flu_costa", ordered_lins)
  ]

  fill_breaks <- ordered_lins
  fill_labels <- ordered_lins

  heatmap_legend_key_size <- grid::unit(1.1, "cm")

  draw_flu_fill_key <- function(data, params, size) {
    grid::rectGrob(
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      gp = grid::gpar(fill = data$fill, col = "white", linewidth = 0.4)
    )
  }

  plot_data <- dplyr::bind_rows(clean_metadata, costa_panel_markers)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample_ID, y = Segment, fill = Lineage_Number)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      data = costa_panel_markers,
      ggplot2::aes(label = star),
      color = costa_star_color,
      size = 5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      breaks = fill_breaks,
      labels = fill_labels,
      name = legend_title,
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(
      labels = c(setNames(rev(segments), rev(segments)), Panel = "")
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(
      title = title,
      x = "Sample (EPI_ISL)",
      y = "Segment"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        size = c(4, 4, 4, 4, 12, 4, 4, 12, 4),
        face = c("plain", "plain", "plain", "plain", "bold", "plain", "plain", "bold", "plain")
      ),
      axis.title.x = ggplot2::element_blank(),
      axis.title.x.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 8, face = "bold"),
      legend.title = ggplot2::element_text(size = 8, face = "bold", angle = 90),
      legend.text = ggplot2::element_text(size = 8),
      legend.position = "right",
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title.position = "left",
        title.hjust = 0.5,
        title.vjust = 0.5,
        keywidth = heatmap_legend_key_size,
        keyheight = heatmap_legend_key_size,
        draw.key = draw_flu_fill_key
      )
    )

  p <- compose_heatmap_with_legend(p, costa_star_color)

  list(
    plot = p,
    palette = palette,
    sample_order = sample_order,
    costa_marker_fill = costa_marker_fill
  )
}

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(RColorBrewer)
    library(cowplot)
  })

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop(
      "Usage: Rscript flu_lineage_heatmap.R ",
      "<metadata_csv> [output_png]"
    )
  }

  metadata_path <- args[1]
  output_path <- if (length(args) >= 2) args[2] else "flu_lineage.png"

  heatmap_obj <- build_flu_lineage_heatmap(metadata_path = metadata_path)

  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(
    filename = output_path,
    plot = heatmap_obj$plot,
    width = 12,
    height = 3.2,
    dpi = 300,
    bg = "white"
  )

  rds_path <- sub("\\.[^.]+$", ".rds", output_path)
  saveRDS(heatmap_obj, rds_path)

  cat("Plot saved to", output_path, "\n")
  cat("Heatmap object saved to", rds_path, "\n")
}
