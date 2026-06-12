library(ape)
meta <- read.csv("metadata/H5N1_context.csv", stringsAsFactors = FALSE)
flu_tips_all <- meta$file_name[grepl("^flu_", meta$expected_role)]

for (seg in c("HA", "PB2")) {
  tree <- read.tree(paste0("results/phylogeny/iq-tree/", seg, "/lineage/H5N1_", seg, "_fast.treefile"))
  flu <- intersect(tree$tip.label, flu_tips_all)
  if(length(flu) > 1) {
    mrca <- getMRCA(tree, flu)
    desc <- length(phangorn::Descendants(tree, mrca, "tips")[[1]])
    cat("Seg:", seg, "Flu tips:", length(flu), "Tips in MRCA clade:", desc, "Total tree tips:", Ntip(tree), "\n")
  }
}
