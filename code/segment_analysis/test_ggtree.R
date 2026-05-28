cat << 'EOF' > code/segment_analysis/test_ggtree.R
library(ape)
library(ggplot2)
library(ggtree)

args <- c("results/phylogeny/iq-tree/PA/PA_final.treefile", "results/phylogeny/itol_styles/PA/tree_colors.txt", "results/figures/PA_tree_test.png", "PA Tree")

tree_file <- args[1]
colors_file <- args[2]
output_file <- args[3]
title <- args[4]

tree <- read.tree(tree_file)
lines <- readLines(colors_file)
data_idx <- which(lines == "DATA")
data_lines <- lines[(data_idx + 1):length(lines)]

taxon_map <- list()
for (line in data_lines) {
  line <- trimws(line)
  if (line == "" || grepl("^#", line)) next
  parts <- strsplit(line, ",")[[1]]
  if (length(parts) >= 3 && parts[2] == "label") {
    taxon_map[[parts[1]]] <- parts[3]
  }
}

tips <- tree$tip.label
tip_colors <- sapply(tips, function(x) {
  col <- taxon_map[[x]]
  if (is.null(col)) "#555555" else col
})

flu_tips <- tips[grepl("^Flu-", tips)]
if (length(flu_tips) > 1) {
  flu_mrca <- getMRCA(tree, flu_tips)
  flu_clade_tips <- extract.clade(tree, flu_mrca)$tip.label
} else if (length(flu_tips) == 1) {
  flu_clade_tips <- flu_tips
} else {
  flu_clade_tips <- character(0)
}

shape_vec <- rep(16, length(tips))
names(shape_vec) <- tips
shape_vec[tips %in% flu_clade_tips] <- 2
shape_vec[tips %in% flu_tips] <- 17

annotation <- data.frame(
  label = tips,
  color = tip_colors,
  shape = shape_vec,
  stringsAsFactors = FALSE
)

p <- ggtree(tree, color="transparent") %<+% annotation +
  geom_tippoint(aes(color = color, shape = shape), size = 6) +
  scale_color_identity() +
  scale_shape_identity() +
  theme_void() +
  ggtitle(title) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
fig_height <- max(6, length(tips) * 0.25)
ggsave(output_file, plot = p, width = 12, height = fig_height, dpi = 300)
cat(paste("Successfully saved tree plot to", output_file, "\n"))
EOF
