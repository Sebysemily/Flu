# plot_cophylo_tanglegram.R
# phytools::cophylo tanglegram renderer for HA/PB2 composite figures.

set.seed(3980)
source("code/segment_analysis/tree_aesthetics.R")
source("code/segment_analysis/tanglegram_legend.R")

tip_style_for_phylo <- function(tip_labels, roles, lineages, palette) {
  n <- length(tip_labels)
  pch <- rep(16L, n)
  col <- rep("grey80", n)
  cex <- rep(COPHYLO_CONTEXT_TIP_CEX, n)

  for (i in seq_len(n)) {
    tip <- tip_labels[i]
    role <- normalize_ecuador_role(roles[tip])
    is_ec <- is_ecuador_tip(roles[tip])
    
    genotype <- lineages[tip]
    style <- get_custom_genotype_style(genotype, is_ec, palette)

    # Determine transparency (alpha) and size (cex)
    if (is_ec) {
      col[i] <- scales::alpha(style$color, 1.0) # Ecuador has no alpha
      cex[i] <- COPHYLO_FLU_TIP_CEX
    } else {
      col[i] <- scales::alpha(style$color, 0.7) # Rest have alpha 0.7
      cex[i] <- COPHYLO_CONTEXT_TIP_CEX
    }

    pch[i] <- style$pch
  }

  list(pch = pch, col = col, cex = cex)
}

draw_cophylo_tip_markers <- function(plot_side, style) {
  n <- plot_side$Ntip
  points(
    plot_side$xx[seq_len(n)],
    plot_side$yy[seq_len(n)],
    pch = style$pch,
    col = style$col,
    cex = style$cex
  )
}

render_cophylo_tanglegram <- function(
    tree_ha,
    tree_pb2,
    tip_roles,
    tip_lineages_ha,
    tip_lineages_pb2,
    lineage_palette,
    link_lwd = 4,
    edge_lwd = 1.2
) {
  obj <- phytools::cophylo(tree_ha, tree_pb2, rotate = FALSE)

  link_col <- vapply(
    obj$assoc[, 1],
    function(tip) ribbon_color_for_type(tip_roles[tip]),
    character(1)
  )

  left_style <- tip_style_for_phylo(
    obj$trees[[1]]$tip.label,
    tip_roles,
    tip_lineages_ha,
    lineage_palette
  )
  right_style <- tip_style_for_phylo(
    obj$trees[[2]]$tip.label,
    tip_roles,
    tip_lineages_pb2,
    lineage_palette
  )

  phytools::plot.cophylo(
    obj,
    link.type = "curved",
    link.lwd = link_lwd,
    link.lty = 1,
    link.col = link_col,
    tip.lty = 0, # Turn off default black dotted tip lines
    fsize = 0,   # Remove gap for labels
    tip.len = 0, # Remove default tip extension gap
    ftype = "off",
    pts = FALSE,
    frame = FALSE,
    edge.width = edge_lwd,
    scale.bar = c(0, 0)
  )

  plot_obj <- get("last_plot.cophylo", envir = ape::.PlotPhyloEnv)

  # Function to reduce alpha of a vector of colors by a factor
  reduce_alpha <- function(cols, factor = 0.5) {
    unname(vapply(cols, function(c) {
      rgba <- col2rgb(c, alpha = TRUE)
      rgb(rgba[1,1], rgba[2,1], rgba[3,1], rgba[4,1] * factor, maxColorValue = 255)
    }, character(1)))
  }

  # Extract the start points of the ribbons so we can draw colored lines to them
  # In phytools plot.cophylo, the tip lines go from xx[i] to a specific x coordinate.
  # Instead of calculating the complex phytools internal offset, we can extract the bounds 
  # from the tree and draw horizontal segments up to the ribbon start.
  # The ribbon starts are symmetric around 0.
  # We know the rightmost point of the left tree is max(plot_obj$left$xx)
  # but phytools adds a gap. The simplest way to perfectly match the color is to 
  # draw a horizontal line from the tip's X to a computed boundary.
  
  # Wait, an easier and more robust way is to compute the ribbon boundary exactly 
  # as phytools does, or just find where the tip label space ends.
  # By default part=0.4, tip.len=0.1, fsize=1.0. 
  # Actually, we can just trace the bounding box, or just draw to x=0 for left and x=0 for right, 
  # but we draw BEFORE the tip markers and UNDER the ribbons? 
  # No, we just overdraw. But the ribbons have alpha=0.3! 
  # If we draw to X=0, it overlaps the ribbon and causes darker lines!
  # So we must draw exactly up to the ribbon edge.
  # For left tree, the edge is at x = -0.1 (approx).
  
  # Let's compute the phytools gap:
  # In phytools, x_left_ribbon is usually around -0.05 or something.
  # Let's find the max X of the left tree tips:
  # x_left_max <- max(plot_obj$left$xx)
  # The ribbon starts exactly at x_left_max + gap.
  # Actually, we can just use the link_col vector. 
  # Note: link_col is in the order of obj$assoc. 
  # So we need to map obj$assoc to the left and right tips.
  
  left_tip_order <- match(obj$trees[[1]]$tip.label, obj$assoc[, 1])
  left_line_cols <- link_col[left_tip_order]
  
  right_tip_order <- match(obj$trees[[2]]$tip.label, obj$assoc[, 2])
  right_line_cols <- link_col[right_tip_order]
  
  # By setting fsize=0 and tip.len=0 in plot.cophylo, the ribbons now start EXACTLY 
  # at max(plot_obj$left$xx) for the left tree and min(plot_obj$right$xx) for the right tree.
  x_ribbon_left <- max(plot_obj$left$xx)
  x_ribbon_right <- min(plot_obj$right$xx)
  
  # Use lower thickness (lwd = 1.0) to reduce clutter
  segments(
    x0 = plot_obj$left$xx, 
    y0 = plot_obj$left$yy, 
    x1 = x_ribbon_left, 
    y1 = plot_obj$left$yy, 
    col = left_line_cols, 
    lty = "solid", lwd = 1.0
  )
  
  segments(
    x0 = plot_obj$right$xx, 
    y0 = plot_obj$right$yy, 
    x1 = x_ribbon_right, 
    y1 = plot_obj$right$yy, 
    col = right_line_cols, 
    lty = "solid", lwd = 1.0
  )

  draw_cophylo_tip_markers(plot_obj$left, left_style)
  draw_cophylo_tip_markers(plot_obj$right, right_style)

  list(
    obj = obj,
    tip_order = obj$assoc[, 1]
  )
}
