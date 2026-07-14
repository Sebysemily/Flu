library(ggplot2)
library(dplyr)
library(maps)

# Simulate tree bounds
max_tree_x <- 2024
tree_y_max <- 367

world_map <- map_data("world")
sa_countries <- c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", "Peru", "Uruguay", "Venezuela", "Falkland Islands")
sa_map <- world_map %>% filter(region %in% sa_countries)

# Bounding box of SA
min_lat <- min(sa_map$lat)
max_lat <- max(sa_map$lat)
min_lon <- min(sa_map$long)
max_lon <- max(sa_map$long)

# Scale SA map to fit tree Y axis
scale_factor <- tree_y_max / (max_lat - min_lat) * 0.8  # 80% of tree height

# Transform coordinates
sa_map <- sa_map %>% mutate(
  trans_y = (lat - min_lat) * scale_factor + (tree_y_max * 0.1),
  trans_x = (long - min_lon) * scale_factor + max_tree_x + 5
)

# Centroids
centroids <- sa_map %>% group_by(region) %>% summarize(
  c_x = mean(trans_x),
  c_y = mean(trans_y)
)

print(head(centroids))
