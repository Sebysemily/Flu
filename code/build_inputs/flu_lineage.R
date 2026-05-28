#!/usr/bin/env Rscript

# Flu Lineage Constellation Plot
# Generates a heatmap-like tile plot of segment lineages per sample.
# Usage: Rscript flu_lineage.R <metadata_csv> [output_png]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
})

source("code/segment_analysis/flu_lineage_heatmap.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide path to metadata CSV (first argument). Optional second argument: output PNG path.")
}
metadata_path <- args[1]
output_path <- ifelse(length(args) >= 2, args[2], "flu_lineage.png")

heatmap_obj <- build_flu_lineage_heatmap(metadata_path = metadata_path)
p <- heatmap_obj$plot

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = output_path,
  plot = p,
  width = 12,
  height = 3.2,
  dpi = 300,
  bg = "white"
)

rds_path <- sub("\\.[^.]+$", ".rds", output_path)
saveRDS(heatmap_obj, rds_path)

cat("Plot successfully saved to", output_path, "\n")
cat("Heatmap object saved to", rds_path, "\n")
