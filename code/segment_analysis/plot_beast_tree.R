#!/usr/bin/env Rscript
# plot_beast_tree.R
# Plots a BEAST MCC tree with tip colors by geography (expected_role),
# similar to the 8-segments collapsed figures, with a time scale.
# Usage: Rscript plot_beast_tree.R <mcc.tree> <metadata.tsv> <output.png> [title]

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(dplyr)
  library(readr)
  library(phytools)
  library(readr)
  library(maps)
})

source("code/segment_analysis/tree_aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_beast_tree.R <mcc.tree> <metadata.tsv> <output.png> [title]")
}

tree_file     <- args[1]
metadata_file <- args[2]
output_file   <- args[3]
title         <- if (length(args) >= 4) args[4] else "BEAST Time-Scaled Tree"

tree <- treeio::read.beast(tree_file)
tree_phylo <- tree@phylo

if (any(tree_phylo$edge.length < 0, na.rm = TRUE)) {
  tree_phylo$edge.length[tree_phylo$edge.length < 0] <- 0.00001
  tree@phylo <- tree_phylo
}

meta <- read_tsv(metadata_file, show_col_types = FALSE)

tip_data <- data.frame(label = tree_phylo$tip.label, stringsAsFactors = FALSE) |>
  left_join(meta, by = c("label" = "Taxon")) |>
  rename(host_type = Host) |>
  mutate(
    # Fallback in case expected_role is missing in TSV (we rely on Location for now)
    role = Location
  )

# Fix role mapping for Ecuador
tip_data$role <- ifelse(tip_data$role == "ecuador_coastal", "flu_costa", tip_data$role)
tip_data$role <- ifelse(tip_data$role == "ecuador_andine", "flu_andine", tip_data$role)
tip_data$role <- ifelse(tip_data$role == "ecuador_amazonia", "flu_amazonia", tip_data$role)

tip_data$is_flu <- grepl("^flu_", tip_data$role)
tip_data$pt_size <- 8.0
tip_data$pt_alpha <- case_when(
  tip_data$is_flu ~ 1.0,
  tip_data$role == "USA" ~ 0.2, # Extremely low for USA (the anchor)
  TRUE ~ 0.9  # High alpha for warm country colors
)

valid_dates <- as.Date(tip_data$Date, format="%Y-%m-%d")
valid_dates <- valid_dates[!is.na(valid_dates)]
mrsd_date <- if (length(valid_dates) > 0) max(valid_dates) else NULL

p <- ggtree(tree, mrsd = mrsd_date, color = "grey35", linewidth = 1.5) %<+% tip_data

# --- COMPRESS USA TIPS VERTICALLY ---
usa_nodes <- which(p$data$role == "USA")
if (length(usa_nodes) > 0) {
  min_usa <- min(p$data$y[usa_nodes], na.rm=TRUE)
  max_usa <- max(p$data$y[usa_nodes], na.rm=TRUE)
  # Compress by 80%
  comp_ratio <- 0.2
  affected <- which(p$data$y >= min_usa & p$data$y <= max_usa)
  p$data$y[affected] <- min_usa + (p$data$y[affected] - min_usa) * comp_ratio
  # Shift the rest of the tree to close the gap
  if (min_usa == 1) { # USA is at the bottom
    shift_amt <- (max_usa - min_usa) * (1 - comp_ratio)
    above <- which(p$data$y > max_usa)
    p$data$y[above] <- p$data$y[above] - shift_amt
  }
}

diff_x <- max(p$data$x, na.rm = TRUE) - min(p$data$x, na.rm = TRUE)
tip_offset <- diff_x * 0.0015

if (nrow(tip_data) > 0) {
  p <- p + geom_point(
    aes(x = x + tip_offset, y = y, color = role, shape = host_type, size = pt_size, alpha = pt_alpha),
    stroke = 1.0
  ) +
    scale_color_manual(values = panel_type_colors, guide = "none", na.value = "grey70", na.translate = FALSE) +
    scale_shape_manual(values = c("domesticated bird"=15, "wild bird"=16, "domesticated mammal"=17, "wild mammal"=18, "?"=3), name = "Host Type", na.value=16, na.translate = FALSE) +
    scale_size_identity() + scale_alpha_identity() + scale_linewidth_identity()
}

if (!is.null(tree@data) && "posterior" %in% names(tree@data)) {
  # Calculate dynamic node size based on clade vertical separation
  node_y <- setNames(p$data$y, as.character(p$data$node))
  p$data$y_span <- sapply(p$data$node, function(n) {
    children <- p$data$node[p$data$parent == n & p$data$node != n]
    if (length(children) > 0) {
      abs(diff(range(node_y[as.character(children)])))
    } else {
      0
    }
  })
  max_span <- max(p$data$y_span, na.rm=TRUE)
  p$data$node_size <- if(max_span > 0) 1.5 + (p$data$y_span / max_span) * 4.5 else 3.0

  p$data$support_class <- factor(
    case_when(
      p$data$posterior >= 0.90 ~ "≥ 90%",
      p$data$posterior >= 0.80 & p$data$posterior < 0.90 ~ "80% - 89%",
      TRUE ~ NA_character_
    ),
    levels = c("≥ 90%", "80% - 89%")
  )
  p <- p + 
    geom_nodepoint(aes(fill = support_class, size = node_size), shape = 21, color = "black", stroke = 1.2, na.rm = TRUE) +
    scale_fill_manual(values = c("≥ 90%" = "black", "80% - 89%" = "white"), name = "Posterior Support", na.translate = FALSE)
}

p <- p + theme_tree() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.0, 0.95),
    legend.justification = c(0, 1),
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.direction = "vertical",
    legend.title = element_text(face="bold", size=58),
    legend.text = element_text(size=50),
    legend.key.size = unit(4.0, "cm"),
    plot.margin = margin(t=5, r=5, b=40, l=35)
  ) +
  guides(
    shape = guide_legend(ncol = 2, override.aes = list(size = 28)),
    fill = guide_legend(ncol = 2, override.aes = list(size = 24))
  )

# --- ADD SOUTH AMERICA MAP AND CONNECTING LINES ---

# --- ADD SOUTH AMERICA MAP AND CONNECTING LINES ---

# 1. Get South America map data
map_countries <- c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", "Falkland Islands", "Peru", "Uruguay", "Venezuela")
line_countries <- c(map_countries, "Antarctica")
world_map <- map_data("world")
sa_map <- world_map %>% filter(region %in% map_countries)

# 2. Calculate scaling and perfectly proportional dimensions
min_tree_x <- min(p$data$x, na.rm = TRUE)
max_tree_x <- max(p$data$x, na.rm = TRUE)
tree_y_max <- max(p$data$y, na.rm = TRUE)
diff_x <- max_tree_x - min_tree_x

# Aspect ratio compensation removed since fig_w auto-adjusts
min_lat <- min(sa_map$lat)
max_lat <- max(sa_map$lat)
min_lon <- min(sa_map$long)
max_lon <- max(sa_map$long)

map_target_height <- tree_y_max * 0.80
scale_y <- map_target_height / (max_lat - min_lat)

map_target_width <- diff_x * 0.605
scale_x <- map_target_width / (max_lon - min_lon)

# 3. Transform map coordinates
# Map sits deeply inside the right edge of the tree for maximum compactness
x_offset_map <- max_tree_x - diff_x * 0.07

# Elevate the map entirely so it sits perfectly relative to the axis
map_y_shift <- tree_y_max * 0.06

sa_map <- sa_map %>% mutate(
  trans_y = (lat - min_lat) * scale_y + map_y_shift,
  trans_x = (long - min_lon) * scale_x + x_offset_map
)

# 4. Calculate country centroids on the transformed map
# Filter out Galapagos (long < -85) so Ecuador's centroid doesn't fall in the ocean
sa_map_mainland <- sa_map %>% filter(!(region == "Ecuador" & long < -85))

centroids <- sa_map_mainland %>% group_by(region) %>% summarize(
  c_x = mean(trans_x),
  c_y = mean(trans_y)
)

# Manually add a fake centroid for Antarctica below South America
ant_c_x <- mean(sa_map_mainland$trans_x, na.rm=TRUE)
ant_c_y <- min(sa_map_mainland$trans_y, na.rm=TRUE) - (tree_y_max * 0.08)

centroids <- bind_rows(centroids, data.frame(
  region = "Antarctica",
  c_x = ant_c_x,
  c_y = ant_c_y
))

# Create a fake Antarctica with horizontal gradient using vertical segments
# Make the right side much steeper/shorter to free up ocean space for the geographic legend
ant_radius_x_left <- (max(sa_map_mainland$trans_x, na.rm=T) - min(sa_map_mainland$trans_x, na.rm=T)) * 0.45
ant_radius_x_right <- (max(sa_map_mainland$trans_x, na.rm=T) - min(sa_map_mainland$trans_x, na.rm=T)) * 0.15

x_vals <- seq(ant_c_x - ant_radius_x_left, ant_c_x + ant_radius_x_right, length.out = 400)
base_y <- ant_c_y - tree_y_max * 0.05
peak_y <- ant_c_y + tree_y_max * 0.015 # Ensure the point ant_c_y is embedded inside the ice
noise <- sin(seq(0, 8*pi, length.out=400)) * (tree_y_max * 0.003)

y_vals <- ifelse(x_vals < ant_c_x,
                 base_y + (peak_y - base_y) * (1 - (abs(x_vals - ant_c_x)/ant_radius_x_left)^2) + noise,
                 base_y + (peak_y - base_y) * (1 - (abs(x_vals - ant_c_x)/ant_radius_x_right)^2) + noise)

alpha_vals <- ifelse(x_vals < ant_c_x,
                     1 - (abs(x_vals - ant_c_x)/ant_radius_x_left)^1.5,
                     1 - (abs(x_vals - ant_c_x)/ant_radius_x_right)^1.5)

fake_antarctica_grad <- data.frame(
  x = x_vals,
  y = base_y,
  yend = y_vals,
  alpha = alpha_vals
)

# Map our custom location names to the map regions
tip_coords <- p$data %>% 
  filter(isTip == TRUE, role %in% line_countries | role %in% c("flu_costa", "flu_andine", "flu_amazonia")) %>%
  mutate(
    map_region = case_when(
      role %in% c("flu_costa", "flu_andine", "flu_amazonia") ~ "Ecuador",
      TRUE ~ role
    )
  ) %>%
  left_join(centroids, by = c("map_region" = "region"))

# Introduce geographical offsets for Ecuador sub-regions so lines don't converge to a single point
tip_coords <- tip_coords %>% mutate(
  c_y = case_when(
    role == "flu_costa" ~ c_y + (map_target_height * 0.015), # Shift Coast up (North) slightly
    role == "flu_andine" ~ c_y - (map_target_height * 0.015), # Shift Andine down (South) slightly
    TRUE ~ c_y
  ),
  c_x = case_when(
    role == "flu_costa" ~ c_x - (map_target_width * 0.015), # Shift Coast left (West) slightly so it hits land, not ocean
    role == "flu_andine" ~ c_x + (map_target_width * 0.005), # Shift Andine right slightly
    role == "Peru" ~ c_x - (map_target_width * 0.06), # Shift Peru heavily to the left (West)
    role == "Brazil" ~ c_x + (map_target_width * 0.15), # Mover Brasil más a la derecha
    role == "Argentina" ~ c_x + (map_target_width * 0.05), # Mover Argentina levemente a la derecha
    TRUE ~ c_x
  ),
  line_alpha = ifelse(grepl("^flu_", role), 1.0, 0.75),
  line_width = ifelse(grepl("^flu_", role), 2.0, 0.45)
) %>%
  filter(!is.na(c_x))

# 5. Add Map, Fake Antarctica, and Lines to the plot
p <- p + 
  # Fake Antarctica landmass with gradient
  geom_segment(
    data = fake_antarctica_grad,
    aes(x = x, xend = x, y = y, yend = yend, alpha = alpha),
    color = "grey75", linewidth = 0.6, inherit.aes = FALSE
  ) +
  # 1st layer: Base Map
  geom_polygon(
    data = sa_map,
    aes(x = trans_x, y = trans_y, group = group),
    fill = "grey90", color = "white", linewidth = 0.4, inherit.aes = FALSE
  ) +
  # 2nd layer: Add horizontal connecting lines (Tip to edge)
  # The line starts extremely close to the tip shape to visually connect to it
  geom_segment(
    data = tip_coords,
    aes(x = x + (diff_x * 0.003), y = y, xend = max_tree_x, yend = y, color = role, alpha = line_alpha, linewidth = line_width),
    inherit.aes = FALSE
  ) +
  # 4th layer: Add angled connecting lines (Edge to map)
  geom_segment(
    data = tip_coords,
    aes(x = max_tree_x, y = y, xend = c_x, yend = c_y, color = role, alpha = line_alpha, linewidth = line_width),
    inherit.aes = FALSE
  )

# Add custom colored text legend for countries
# Moved to the top-left, just below the metadata legend
roles_present <- sort(unique(tip_coords$role))
custom_legend <- data.frame(
  role = roles_present,
  y = seq(tree_y_max * 0.68, by = -tree_y_max * 0.035, length.out = length(roles_present))
)
p <- p + geom_text(
  data = custom_legend,
  aes(x = -Inf, y = y, label = role, color = role),
  fontface = "bold", size = 26, hjust = 0, inherit.aes = FALSE
)

# Force the X-axis to only appear under the tree (not the map) manually
b_ticks <- scales::breaks_pretty(n = 6)(c(min_tree_x, max_tree_x))
b_ticks <- b_ticks[b_ticks <= max_tree_x & b_ticks >= min_tree_x]

axis_y <- -tree_y_max * 0.02
p <- p + 
  # Main axis line from the exact start to exact end of the tree
  geom_segment(aes(x = min_tree_x, xend = max_tree_x, y = axis_y, yend = axis_y), color = "black", linewidth = 1.5, inherit.aes = FALSE) +
  # Ticks
  geom_segment(data = data.frame(x = b_ticks), aes(x = x, xend = x, y = axis_y, yend = axis_y - tree_y_max * 0.015), color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  # Tick Labels
  geom_text(data = data.frame(x = b_ticks), aes(x = x, y = axis_y - tree_y_max * 0.04), label = b_ticks, size = 14, color = "black", inherit.aes = FALSE)

# Add Floating Title to compact the layout vertically
p <- p + annotate(
  "text",
  x = min_tree_x + (max_tree_x - min_tree_x) * 0.65,
  y = tree_y_max * 0.98,
  label = title,
  size = 24, fontface = "bold", hjust = 0.5
)

# 6. Save the final combined plot
dir.create(dirname(output_file), showWarnings=FALSE, recursive=TRUE)

# Calculate mathematically exact fig_w to prevent map from squishing horizontally
fig_h <- min(35, max(8, Ntip(tree_phylo) * 0.10))
diff_x_total <- diff_x * 1.1 + map_target_width
# 1 degree lat = scale_y * (fig_h / tree_y_max) inches
# 1 degree lon = scale_x * (fig_w / diff_x_total) inches
# We set them equal to maintain aspect ratio 1:1
fig_w <- (scale_y * fig_h / tree_y_max) * diff_x_total / scale_x

ggsave(output_file, plot=p, width=fig_w, height=fig_h, dpi=300, limitsize=FALSE, bg = "white")
