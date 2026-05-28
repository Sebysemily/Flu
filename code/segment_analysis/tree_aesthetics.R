# tree_aesthetics.R
# Tip and ribbon styling for Ecuador + context samples in tanglegram figures.

flu_tip_colors <- c(
  flu_costa  = "#FF0000",
  flu_sierra = "#00008B"
)

flu_tip_shapes <- c(
  flu_costa  = 17L,
  flu_sierra = 17L
)

flu_tip_labels <- c(
  flu_costa  = "Ecuador (coastal)",
  flu_sierra = "Ecuador (sierra)"
)

FLU_TIP_ROLES <- names(flu_tip_colors)

FLU_TIP_CEX <- 7
CONTEXT_TIP_CEX <- 5
CONTEXT_TIP_SHAPE <- 16L

is_ecuador_tip <- function(role) {
  !is.na(role) && role %in% FLU_TIP_ROLES
}

panel_type_colors <- c(
  flu_costa        = "#FF0000",
  flu_sierra       = "#00008B",
  american_anchor  = "#D95F02",
  regional_context = "#1B9E77"
)

panel_type_labels <- c(
  flu_costa        = "Ecuador (coastal)",
  flu_sierra       = "Ecuador (sierra)",
  american_anchor  = "American anchor",
  regional_context = "Regional context"
)

PANEL_TYPE_ROLES <- names(panel_type_colors)

PANEL_TYPE_RIBBON_ORDER <- c(
  "flu_costa",
  "flu_sierra",
  "american_anchor",
  "regional_context"
)

ribbon_color_for_type <- function(role, flu_alpha = 0.90, context_alpha = 0.32) {
  role <- as.character(role)
  if (is.na(role) || role == "" || !(role %in% PANEL_TYPE_ROLES)) {
    return(scales::alpha("grey70", context_alpha))
  }
  base <- unname(panel_type_colors[role])
  if (role %in% FLU_TIP_ROLES) {
    return(scales::alpha(base, flu_alpha))
  }
  scales::alpha(base, context_alpha)
}
