lines <- readLines("code/segment_analysis/plot_8_segments.R")

start_idx <- grep("bottom_row <- plot_grid", lines)
end_idx <- grep("combined_legends <- plot_grid", lines)

new_code <- c(
  "grid_2x2 <- plot_grid(",
  "  leg_reassortant, leg_role,",
  "  leg_shape, NULL,",
  "  nrow = 2,",
  "  rel_widths = c(0.4, 1.2),",
  "  rel_heights = c(1, 1)",
  ")",
  "combined_legends <- plot_grid(NULL, grid_2x2, nrow=1, rel_widths=c(0.12, 1))"
)

new_lines <- c(lines[1:(start_idx-1)], new_code, lines[(end_idx+1):length(lines)])
writeLines(new_lines, "code/segment_analysis/plot_8_segments.R")
