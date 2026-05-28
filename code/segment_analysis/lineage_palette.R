# lineage_palette.R
# Shared GenoFLU sub-lineage colors: am* = Oranges, ea* = Blues, other = Greens.

build_lineage_palette <- function(
    lineages,
    unknown_label = "unknown",
    costa_marker_fill = NULL,
    costa_marker_color = "white"
) {
  lineages <- unique(as.character(lineages))
  lineages <- lineages[!is.na(lineages) & lineages != ""]

  unknown_lineages <- lineages[tolower(lineages) == unknown_label]
  am_lineages <- sort(lineages[grepl("^am", lineages, ignore.case = TRUE)])
  ea_lineages <- sort(lineages[grepl("^ea", lineages, ignore.case = TRUE)])
  other_lineages <- sort(setdiff(
    lineages,
    c(unknown_lineages, am_lineages, ea_lineages)
  ))

  palette <- c()
  if (length(unknown_lineages) > 0) {
    unknown_colors <- rep("grey85", length(unknown_lineages))
    names(unknown_colors) <- unknown_lineages
    palette <- c(palette, unknown_colors)
  }
  if (length(am_lineages) > 0) {
    am_base <- RColorBrewer::brewer.pal(
      max(3, min(9, length(am_lineages))),
      "Oranges"
    )
    am_colors <- grDevices::colorRampPalette(am_base)(length(am_lineages))
    names(am_colors) <- am_lineages
    palette <- c(palette, am_colors)
  }
  if (length(ea_lineages) > 0) {
    ea_base <- RColorBrewer::brewer.pal(
      max(3, min(9, length(ea_lineages))),
      "Blues"
    )
    ea_colors <- grDevices::colorRampPalette(ea_base)(length(ea_lineages))
    names(ea_colors) <- ea_lineages
    palette <- c(palette, ea_colors)
  }
  if (length(other_lineages) > 0) {
    other_base <- RColorBrewer::brewer.pal(
      max(3, min(8, length(other_lineages))),
      "Greens"
    )
    other_colors <- grDevices::colorRampPalette(other_base)(length(other_lineages))
    names(other_colors) <- other_lineages
    palette <- c(palette, other_colors)
  }

  if (!is.null(costa_marker_fill) && nzchar(costa_marker_fill)) {
    palette[costa_marker_fill] <- costa_marker_color
  }

  palette
}

order_lineages_for_legend <- function(lineages, unknown_label = "unknown") {
  lineages <- unique(as.character(lineages))
  lineages <- lineages[!is.na(lineages) & lineages != ""]

  unknown <- lineages[tolower(lineages) == unknown_label]
  am <- sort(lineages[grepl("^am", lineages, ignore.case = TRUE)])
  ea <- sort(lineages[grepl("^ea", lineages, ignore.case = TRUE)])
  other <- sort(setdiff(lineages, c(unknown, am, ea)))

  c(am, ea, other, unknown)
}

extend_lineage_palette <- function(
    base_palette,
    extra_lineages,
    costa_marker_fill = ".flu_costa_panel"
) {
  base_names <- setdiff(names(base_palette), costa_marker_fill)
  extra <- unique(as.character(extra_lineages))
  extra <- extra[!is.na(extra) & extra != ""]
  extra <- setdiff(extra, costa_marker_fill)

  all_lineages <- unique(c(base_names, extra))
  full_palette <- build_lineage_palette(all_lineages)

  shared <- intersect(base_names, names(full_palette))
  full_palette[shared] <- unname(base_palette[shared])

  if (costa_marker_fill %in% names(base_palette)) {
    full_palette[costa_marker_fill] <- unname(base_palette[costa_marker_fill])
  }

  full_palette
}
