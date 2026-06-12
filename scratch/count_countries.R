library(ape)
library(dplyr)
library(phangorn)

meta <- read.csv("metadata/H5N1_context.csv", stringsAsFactors = FALSE)
flu_tips_all <- meta$file_name[grepl("^flu_", meta$expected_role)]

tree_paths <- c(
  "results/phylogeny/iq-tree/PB2/lineage/H5N1_PB2_fast.treefile",
  "results/phylogeny/iq-tree/PB1/lineage/H5N1_PB1_fast.treefile",
  "results/phylogeny/iq-tree/PA/lineage/H5N1_PA_fast.treefile",
  "results/phylogeny/iq-tree/HA/lineage/H5N1_HA_fast.treefile",
  "results/phylogeny/iq-tree/NP/lineage/H5N1_NP_fast.treefile",
  "results/phylogeny/iq-tree/NA/lineage/H5N1_NA_fast.treefile",
  "results/phylogeny/iq-tree/MP/lineage/H5N1_MP_fast.treefile",
  "results/phylogeny/iq-tree/NS/lineage/H5N1_NS_fast.treefile"
)

cat("Segment\tCountry\tCount\n")
for (p in tree_paths) {
  if (file.exists(p)) {
    seg <- unlist(strsplit(p, "/"))
    seg <- seg[length(seg) - 2]
    tree <- read.tree(p)
    flu <- intersect(tree$tip.label, flu_tips_all)
    if (length(flu) > 1) {
      mrca <- getMRCA(tree, flu)
      desc <- tree$tip.label[phangorn::Descendants(tree, mrca, "tips")[[1]]]
    } else if (length(flu) == 1) {
      desc <- flu
    } else {
      desc <- c()
    }
    
    if (length(desc) > 0) {
      meta_sub <- meta[meta$file_name %in% desc, ]
      counts <- table(meta_sub$country)
      for (c_name in names(counts)) {
        cat(sprintf("%s\t%s\t%d\n", seg, c_name, counts[c_name]))
      }
    }
  }
}
