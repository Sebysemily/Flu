# flu_lineage_heatmap.R
# Shared GenoFLU constellation heatmap builder (Ecuador panel).

source("code/segment_analysis/lineage_palette.R")

build_flu_lineage_heatmap <- function(
    metadata_path,
    sample_order = NULL,
    title = "Avian Influenza Genotype Constellation"
) {
  flu_data <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

  key_col <- if ("EPI_ISL" %in% colnames(flu_data)) {
    "EPI_ISL"
  } else if (ncol(flu_data) >= 2) {
    colnames(flu_data)[2]
  } else {
    colnames(flu_data)[1]
  }

  segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
  unknown_label <- "unknown"
  coastal_epi_isls <- c(
    "EPI_ISL_17973443", "EPI_ISL_17973458", "EPI_ISL_18137671"
  )
  costa_marker_fill <- ".flu_costa_panel"
  legend_title <- "Sub‑Lineage / Status"
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

  usfq_col <- if ("Código USFQ" %in% colnames(flu_data)) {
    "Código USFQ"
  } else if ("Codigo USFQ" %in% colnames(flu_data)) {
    "Codigo USFQ"
  } else {
    NA_character_
  }

  sample_roles <- flu_data %>%
    dplyr::transmute(
      Sample_ID = .data[[key_col]],
      usfq = if (!is.na(usfq_col)) .data[[usfq_col]] else NA_character_,
      expected_role = if ("expected_role" %in% colnames(flu_data)) {
        .data[["expected_role"]]
      } else {
        NA_character_
      }
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      expected_role = dplyr::if_else(
        !is.na(expected_role) & expected_role != "",
        expected_role,
        dplyr::if_else(
          usfq == "Flu-0406" | Sample_ID %in% coastal_epi_isls,
          "flu_costa",
          "flu_sierra"
        )
      )
    )

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
  fill_legend_breaks <- setdiff(names(palette), costa_marker_fill)
  fill_breaks <- c(fill_legend_breaks, costa_marker_fill)
  fill_labels <- c(fill_legend_breaks, "* Ecuador (coastal)")

  plot_data <- dplyr::bind_rows(clean_metadata, costa_panel_markers)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample_ID, y = Segment, fill = Lineage_Number)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      data = costa_panel_markers,
      ggplot2::aes(label = star),
      color = "#FF0000",
      size = 10,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      breaks = fill_breaks,
      labels = fill_labels,
      name = legend_title
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
      axis.text.y = ggplot2::element_text(size = 6),
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
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
        title.vjust = 0.5
      )
    )

  list(
    plot = p,
    palette = palette,
    sample_order = sample_order,
    costa_marker_fill = costa_marker_fill
  )
}
