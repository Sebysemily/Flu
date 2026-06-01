# plot_cophylo_tanglegram.R
# phytools::cophylo tanglegram renderer for HA/PB2 composite figures.

set.seed(3980)
source("code/segment_analysis/tree_aesthetics.R")
source("code/segment_analysis/tanglegram_legend.R")

tip_style_for_phylo <- function(tip_labels, roles, lineages, palette) {
  n <- length(tip_labels)
  pch <- rep(CONTEXT_TIP_SHAPE, n)
  col <- rep("grey80", n)
  cex <- rep(COPHYLO_CONTEXT_TIP_CEX, n)

  for (i in seq_len(n)) {
    tip <- tip_labels[i]
    role <- normalize_ecuador_role(roles[tip])
    if (is_ecuador_tip(roles[tip])) {
      pch[i] <- flu_tip_shapes[role]
      col[i] <- unname(flu_tip_colors[role])
      cex[i] <- COPHYLO_FLU_TIP_CEX
    } else {
      col[i] <- lineage_color_for(lineages[tip], palette)
    }
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
    ftype = "off",
    pts = FALSE,
    frame = FALSE,
    edge.width = edge_lwd,
    scale.bar = c(0, 0)
  )

  plot_obj <- get("last_plot.cophylo", envir = ape::.PlotPhyloEnv)
  draw_cophylo_tip_markers(plot_obj$left, left_style)
  draw_cophylo_tip_markers(plot_obj$right, right_style)

  list(
    obj = obj,
    tip_order = obj$assoc[, 1]
  )
}
