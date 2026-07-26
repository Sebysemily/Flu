# tanglegram_legend.R
# Legend panel for composite tanglegram (ribbons by type, tips by role/lineage).

source("code/segment_analysis/tree_aesthetics.R")
source("code/segment_analysis/lineage_palette.R")

LEGEND_ROW_SPACING <- 2.0
LEGEND_SECTION_GAP <- 2.5
LEGEND_TITLE_OFFSET <- 2.0
LEGEND_ECUADOR_TIP_SIZE <- 7
LEGEND_CONTEXT_TIP_SIZE <- LEGEND_ECUADOR_TIP_SIZE

lineage_color_for <- function(lineage, palette) {
  if (is.na(lineage) || lineage == "" || tolower(lineage) == "unknown") {
    return("grey80")
  }
  if (!is.null(palette) && lineage %in% names(palette)) {
    return(unname(palette[lineage]))
  }
  "grey60"
}

legend_items_df <- function(labels, colors, marker = c("ribbon", "triangle", "circle")) {
  if (length(marker) == 1) {
    marker <- match.arg(marker)
  }
  data.frame(
    label = labels,
    color = colors,
    marker = marker,
    stringsAsFactors = FALSE
  )
}

layout_legend_sections <- function(
    sections,
    row_spacing = LEGEND_ROW_SPACING,
    section_gap = LEGEND_SECTION_GAP,
    title_offset = LEGEND_TITLE_OFFSET
) {
  active <- sections[vapply(
    sections,
    function(s) !is.null(s$items) && nrow(s$items) > 0,
    logical(1)
  )]
  if (length(active) == 0) {
    return(list(sections = list(), y_top = 1, y_bottom = 0))
  }

  n_rows <- sum(vapply(active, function(s) nrow(s$items), integer(1)))
  n_titles <- length(active)
  y_top <- n_titles * (title_offset + section_gap) + n_rows * row_spacing

  y_cursor <- y_top
  laid_out <- vector("list", length(active))

  for (i in seq_along(active)) {
    section <- active[[i]]
    items <- section$items
    n_items <- nrow(items)

    title_y <- y_cursor
    item_ys <- title_y - title_offset - ((seq_len(n_items) - 1) * row_spacing)
    items$y <- item_ys

    laid_out[[i]] <- list(
      title = section$title,
      title_y = title_y,
      items = items
    )

    y_cursor <- min(item_ys) - section_gap
  }

  list(
    sections = laid_out,
    y_top = y_top,
    y_bottom = y_cursor
  )
}

draw_legend_markers <- function(items) {
  layers <- list()
  marker_types <- unique(items$marker)

  if ("ribbon" %in% marker_types) {
    ribbon_items <- items[items$marker == "ribbon", , drop = FALSE]
    layers <- c(
      layers,
      list(
        ggplot2::geom_segment(
          data = ribbon_items,
          ggplot2::aes(x = 0, xend = 1.5, y = y, yend = y),
          linewidth = 4.5,
          lineend = "round",
          color = ribbon_items$color
        )
      )
    )
  }

  if ("triangle" %in% marker_types) {
    triangle_items <- items[items$marker == "triangle", , drop = FALSE]
    layers <- c(
      layers,
      list(
        ggplot2::geom_point(
          data = triangle_items,
          ggplot2::aes(x = 0.5, y = y),
          shape = 17,
          size = LEGEND_ECUADOR_TIP_SIZE,
          color = triangle_items$color
        )
      )
    )
  }

  if ("circle" %in% marker_types) {
    circle_items <- items[items$marker == "circle", , drop = FALSE]
    layers <- c(
      layers,
      list(
        ggplot2::geom_point(
          data = circle_items,
          ggplot2::aes(x = 0.5, y = y),
          shape = CONTEXT_TIP_SHAPE,
          size = LEGEND_CONTEXT_TIP_SIZE,
          color = circle_items$color
        )
      )
    )
  }

  layers
}

draw_legend_section <- function(section) {
  items <- section$items
  c(
    list(
      ggplot2::annotate(
        "text",
        x = 0,
        y = section$title_y,
        label = section$title,
        hjust = 0,
        fontface = "bold",
        size = 6.0
      ),
      ggplot2::geom_text(
        data = items,
        ggplot2::aes(x = 1.8, y = y, label = label),
        hjust = 0,
        size = 5.0
      )
    ),
    draw_legend_markers(items)
  )
}

build_tanglegram_legend <- function(
    ribbon_roles = NULL,
    lineage_col_label = "lineage",
    context_lineages = NULL,
    lineage_palette = NULL,
    flu_alpha = 1,
    context_alpha = 0.30
) {
  if (is.null(ribbon_roles)) {
    show_roles <- PANEL_TYPE_RIBBON_ORDER
  } else {
    show_roles <- PANEL_TYPE_RIBBON_ORDER[
      PANEL_TYPE_RIBBON_ORDER %in% unique(vapply(
        ribbon_roles,
        normalize_ecuador_role,
        character(1)
      ))
    ]
  }

  ribbon_items <- legend_items_df(
    labels = unname(panel_type_labels[show_roles]),
    colors = vapply(
      show_roles,
      function(role) ribbon_color_for_type(role, flu_alpha, context_alpha),
      character(1)
    ),
    marker = "ribbon"
  )

  # Build the tree tips legend to show genotype colors instead of geographic ones.
  b41_color <- lineage_color_for("B4.1", lineage_palette)
  
  labels <- c(
    "B3.2", 
    "B2.2", 
    "B1.3", 
    "B4.1", 
    "Others", 
    "unknown"
  )
  colors <- c(
    scales::alpha("#0000FF", 0.7),      # B3.2 - alpha 0.7 (Blue)
    scales::alpha("#FF0000", 0.7),      # B2.2 - alpha 0.7 (Red)
    scales::alpha("#FF8C00", 0.7),      # B1.3 - alpha 0.7 (Orange)
    scales::alpha(b41_color, 0.7),     # B4.1 - alpha 0.7
    scales::alpha("#78C67A", 0.7),     # Others (B2.1 pastel green) - alpha 0.7
    scales::alpha("grey85", 0.7)       # unknown - alpha 0.7
  )
  markers <- c(
    "circle", 
    "circle", 
    "circle", 
    "circle", 
    "circle", 
    "circle"
  )
  
  tip_items <- legend_items_df(
    labels = labels,
    colors = colors,
    marker = markers
  )

  sections <- list(
    list(title = "Ribbons (geographic groups)", items = ribbon_items),
    list(title = "Tree tips", items = tip_items)
  )

  layout <- layout_legend_sections(sections)
  p <- ggplot2::ggplot()

  for (section in layout$sections) {
    for (layer in draw_legend_section(section)) {
      p <- p + layer
    }
  }

  p +
    ggplot2::coord_cartesian(
      xlim = c(0, 18.0),
      ylim = c(layout$y_bottom - 0.5, layout$y_top + 0.8),
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}
