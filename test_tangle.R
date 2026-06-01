library(ggtree)
library(ape)
library(ggplot2)

set.seed(123)
t1 <- rtree(10)
t2 <- rtree(10)

p1 <- ggtree(t1)
p1 <- collapse(p1, node=12)

d1 <- p1$data
d2 <- fortify(t2)
d2$x <- d2$x + max(d1$x) + 1

# Try to plot
p_combined <- p1 + geom_tree(data=d2)
ggsave("test_tangle.png", p_combined, width=10, height=5)
