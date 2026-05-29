# lineage_palette.R
# Shared GenoFLU sub-lineage colors: am* = Oranges, ea* = Blues, other = Greens.

# Canonical colours for simplified tanglegram / composite display.
COLLAPSED_AM_COLOR <- "#FDAE6B"   # am1.2 orange (reference for all am* tips)
CANONICAL_EA1_COLOR <- "#4292C6"  # lighter blue for all ea* tips (vs andine navy)
COLLAPSED_EA_COLOR <- CANONICAL_EA1_COLOR
AM21_RIBBON_COLOR <- "#E6550D"    # am2.1 — american-anchor ribbon colour

LEGEND_AM_KEY <- ".legend_am"
LEGEND_EA_KEY <- ".legend_ea"

LINEAGE_FAMILY_LABELS <- c(am = "AM", ea = "EA")
LINEAGE_FAMILY_COLORS <- c(
  am = COLLAPSED_AM_COLOR,
  ea = COLLAPSED_EA_COLOR
)

LINEAGE_COLOR_OVERRIDES <- c(
  ea1 = CANONICAL_EA1_COLOR
)

apply_lineage_color_overrides <- function(palette) {
  for (lin in names(LINEAGE_COLOR_OVERRIDES)) {
    if (lin %in% names(palette)) {
      palette[lin] <- LINEAGE_COLOR_OVERRIDES[lin]
    }
  }

  # Keep sub-lineages with distinct shades instead of collapsing
  # am_names <- names(palette)[grepl("^am", names(palette), ignore.case = TRUE)]
  # if (length(am_names) > 0) {
  #   palette[am_names] <- COLLAPSED_AM_COLOR
  # }

  # ea_names <- names(palette)[grepl("^ea", names(palette), ignore.case = TRUE)]
  # if (length(ea_names) > 0) {
  #   palette[ea_names] <- COLLAPSED_EA_COLOR
  # }

  palette
}

lineage_families_present <- function(lineages, unknown_label = "unknown") {
  lineages <- unique(as.character(lineages))
  lineages <- lineages[!is.na(lineages) & lineages != ""]
  lineages <- lineages[tolower(lineages) != unknown_label]

  fams <- character()
  if (any(grepl("^am", lineages, ignore.case = TRUE))) {
    fams <- c(fams, "am")
  }
  if (any(grepl("^ea", lineages, ignore.case = TRUE))) {
    fams <- c(fams, "ea")
  }
  fams
}

lineage_family_legend_breaks <- function(lineages = NULL, unknown_label = "unknown") {
  lineages_vec <- if (is.null(lineages)) {
    character()
  } else {
    vals <- unique(as.character(lineages))
    vals[!is.na(vals) & vals != "" & tolower(vals) != unknown_label]
  }

  fams <- if (length(lineages_vec) == 0) {
    c("am", "ea")
  } else {
    lineage_families_present(lineages_vec, unknown_label = unknown_label)
  }

  keys <- character()
  labels <- character()
  colors <- character()

  if ("am" %in% fams) {
    am_key <- if (length(lineages_vec) > 0) {
      sort(lineages_vec[grepl("^am", lineages_vec, ignore.case = TRUE)])[1]
    } else {
      LEGEND_AM_KEY
    }
    keys <- c(keys, am_key)
    labels <- c(labels, LINEAGE_FAMILY_LABELS["am"])
    colors <- c(colors, LINEAGE_FAMILY_COLORS["am"])
  }

  if ("ea" %in% fams) {
    ea_key <- if (length(lineages_vec) > 0) {
      sort(lineages_vec[grepl("^ea", lineages_vec, ignore.case = TRUE)])[1]
    } else {
      LEGEND_EA_KEY
    }
    keys <- c(keys, ea_key)
    labels <- c(labels, LINEAGE_FAMILY_LABELS["ea"])
    colors <- c(colors, LINEAGE_FAMILY_COLORS["ea"])
  }

  list(
    keys = unname(keys),
    labels = unname(labels),
    colors = unname(colors)
  )
}

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

  apply_lineage_color_overrides(palette)
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

  # For the tanglegram, collapse all am* to one orange and all ea* to one blue
  am_names <- names(full_palette)[grepl("^am", names(full_palette), ignore.case = TRUE)]
  if (length(am_names) > 0) full_palette[am_names] <- COLLAPSED_AM_COLOR

  ea_names <- names(full_palette)[grepl("^ea", names(full_palette), ignore.case = TRUE)]
  if (length(ea_names) > 0) full_palette[ea_names] <- COLLAPSED_EA_COLOR

  full_palette
}
