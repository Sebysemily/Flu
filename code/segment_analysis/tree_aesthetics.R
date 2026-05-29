# tree_aesthetics.R
# Tip and ribbon styling for Ecuador + context samples in tanglegram figures.

source("code/segment_analysis/lineage_palette.R")

flu_tip_colors <- c(
  flu_costa    = "#FF0000",
  flu_andine   = "#00008B",
  flu_amazonia = "#008000"
)

flu_tip_shapes <- c(
  flu_costa    = 17L,
  flu_andine   = 17L,
  flu_amazonia = 17L
)

flu_tip_labels <- c(
  flu_costa    = "Ecuador (coastal)",
  flu_andine   = "Ecuador (andine)",
  flu_amazonia = "Ecuador (amazon)"
)

FLU_TIP_ROLES <- c("flu_costa", "flu_andine", "flu_sierra", "flu_amazonia")
FLU_TIP_DISPLAY_ROLES <- names(flu_tip_colors)

FLU_TIP_CEX <- 7
CONTEXT_TIP_CEX <- 5
CONTEXT_TIP_SHAPE <- 16L

COPHYLO_FLU_TIP_CEX <- 2.2
COPHYLO_CONTEXT_TIP_CEX <- COPHYLO_FLU_TIP_CEX * CONTEXT_TIP_CEX / FLU_TIP_CEX

normalize_ecuador_role <- function(role) {
  role <- as.character(role)
  if (is.na(role) || role == "") {
    return(role)
  }
  if (role == "flu_sierra") {
    return("flu_andine")
  }
  role
}

normalize_role_vector <- function(roles) {
  normalized <- vapply(roles, normalize_ecuador_role, character(1))
  names(normalized) <- names(roles)
  normalized
}

is_ecuador_tip <- function(role) {
  !is.na(role) && role %in% FLU_TIP_ROLES
}

panel_type_colors <- c(
  flu_costa        = "#FF0000",
  flu_andine       = "#00008B",
  flu_amazonia     = "#008000",
  american_anchor  = AM21_RIBBON_COLOR,
  regional_context = "#1B9E77"
)

panel_type_labels <- c(
  flu_costa        = "Ecuador (coastal)",
  flu_andine       = "Ecuador (andine)",
  flu_amazonia     = "Ecuador (amazon)",
  american_anchor  = "American anchor",
  regional_context = "Regional context"
)

PANEL_TYPE_ROLES <- names(panel_type_colors)

PANEL_TYPE_RIBBON_ORDER <- c(
  "flu_costa",
  "flu_andine",
  "flu_amazonia",
  "american_anchor",
  "regional_context"
)

ribbon_color_for_type <- function(role, flu_alpha = 1, context_alpha = 0.30) {
  role <- normalize_ecuador_role(as.character(role))
  if (is.na(role) || role == "" || !(role %in% PANEL_TYPE_ROLES)) {
    return(scales::alpha("grey70", context_alpha))
  }
  base <- unname(panel_type_colors[role])
  if (role %in% FLU_TIP_DISPLAY_ROLES) {
    return(base)
  }
  scales::alpha(base, context_alpha)
}
