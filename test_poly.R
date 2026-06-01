library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)

set.seed(123)
t1 <- rtree(50)

d <- fortify(t1)

# Let's say we collapse node 60
collapse_node <- 60
descendants <- phangorn::Descendants(t1, collapse_node, type="all")

# Polygon coords
x_root <- d$x[d$node == collapse_node]
y_root <- d$y[d$node == collapse_node]

desc_data <- d[d$node %in% descendants, ]
x_max <- max(desc_data$x)
y_min <- min(desc_data$y)
y_max <- max(desc_data$y)

poly <- data.frame(
    x = c(x_root, x_max, x_max),
    y = c(y_root, y_max, y_min),
    group = collapse_node
)

# Remove descendants from d
d_pruned <- d[!d$node %in% descendants, ]

p <- ggplot() +
    geom_tree(data = d_pruned) +
    geom_polygon(data = poly, aes(x=x, y=y, group=group), fill="grey80", color="grey60") +
    theme_tree2()

ggsave("scratch_poly.png", p, width=5, height=5)
