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
  flu_costa        = "#8B0000", # Dark Red
  flu_andine       = "#00008B",
  flu_amazonia     = "#008000",
  american_anchor  = "#FFA500", # Bright Orange
  regional_context = "#1B9E77",
  # Diverse palette for countries
  Peru             = "#A34700", # Naranja muy oscuro/opaco
  Colombia         = "#E91E63", # Pink-ish (rosado ish)
  Venezuela        = "#8D6E63", # Mocha / Brown
  Brazil           = "#FDD835", # Bright Yellow
  Bolivia          = "#5E35B1", # Deep Purple
  Chile            = "#1E88E5", # Strong Blue
  Argentina        = "#29B6F6", # Azul celeste vibrante (Argentina)
  Uruguay          = "#00897B", # Teal
  USA              = "grey70",
  "Falkland Islands"= "grey80",
  Antarctica       = "#B0BEC5"  # Gris Hielo (Ice Gray) para separarlo de los azules
)

panel_type_labels <- c(
  flu_costa        = "Ecuador (coastal)",
  flu_andine       = "Ecuador (andine)",
  flu_amazonia     = "Ecuador (amazon)",
  american_anchor  = "U.S.A",
  regional_context = "South America"
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


get_custom_genotype_style <- function(genotype, is_ec, palette = NULL) {
  if (is.na(genotype) || genotype == "" || tolower(genotype) == "none" || grepl("^not assigned", tolower(trimws(genotype)))) {
    genotype <- "unknown"
  }
  
  geno_clean <- tolower(trimws(genotype))
  
  # Default values
  pch <- 16L
  label <- "Others"
  base_color <- "#78C67A" # Default B2.1 pastel green
  
  if (geno_clean == "b3.2") {
    base_color <- "#0000FF" # Blue
    pch <- 16L
    label <- "B3.2"
  } else if (geno_clean == "b2.2") {
    base_color <- "#FF0000" # Red
    label <- "B2.2"
  } else if (geno_clean == "b1.3") {
    base_color <- "#FF8C00" # Orange/Gold
    label <- "B1.3"
  } else if (geno_clean == "b4.1") {
    base_color <- "#00BFFF" # Celeste / Deep Sky Blue
    label <- "B4.1"
  } else if (geno_clean == "unknown" || geno_clean == "partial genome" || geno_clean == "unassigned" || geno_clean == "partial / unassigned") {
    base_color <- "grey85"
    label <- "Partial / Unassigned"
  }
  
  list(color = base_color, pch = pch, label = label)
}

# --- GLOBAL FONT SETTINGS ---
# The journal requires a non-serif font (Arial/Helvetica). "sans" is the default non-serif.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_get() + ggplot2::theme(text = ggplot2::element_text(family = "sans")))
  ggplot2::update_geom_defaults("text", list(family = "sans"))
  ggplot2::update_geom_defaults("label", list(family = "sans"))
}
