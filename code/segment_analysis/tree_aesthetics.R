# tree_aesthetics.R
# Shared tip styling for pipeline ggtree figures and future manual overlays.

type_colors <- c(
  flu_costa        = "#FF0000",
  flu_sierra       = "#00008B",
  flu_epi_isl      = "#FF0000",
  american_anchor  = "#800080",
  regional_context = "#4CAF50",
  unknown          = "#999999"
)

type_shapes <- c(
  flu_costa        = 17L,
  flu_sierra       = 17L,
  flu_epi_isl      = 17L,
  american_anchor  = 15L,
  regional_context = 16L,
  unknown          = 16L
)

type_sizes <- c(
  flu_costa        = 2.5,
  flu_sierra       = 2.5,
  flu_epi_isl      = 2.5,
  american_anchor  = 2.0,
  regional_context = 1.2,
  unknown          = 1.0
)

type_alphas <- c(
  flu_costa        = 1.0,
  flu_sierra       = 1.0,
  flu_epi_isl      = 1.0,
  american_anchor  = 0.85,
  regional_context = 0.55,
  unknown          = 0.4
)

type_labels <- c(
  flu_costa        = "Ecuador (coastal)",
  flu_sierra       = "Ecuador (sierra)",
  flu_epi_isl      = "Ecuador (coastal EPI)",
  american_anchor  = "American anchor",
  regional_context = "Regional context",
  unknown          = "Unknown"
)

mrk_fill <- c("Flu clade" = "#2ca02c", "Sister clade" = "#1f77b4")

filter_aesthetics_for_types <- function(present_types) {
  list(
    colors = type_colors[present_types],
    shapes = type_shapes[present_types],
    sizes  = type_sizes[present_types],
    alphas = type_alphas[present_types],
    labels = type_labels[present_types]
  )
}
