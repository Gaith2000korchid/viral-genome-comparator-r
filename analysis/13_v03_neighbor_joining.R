# Exploratory Neighbor-Joining tree from the v0.3 whole-genome distance matrix.
# This is a distance tree, not a maximum-likelihood phylogeny.

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Install package 'ape' before running this analysis.")
}

manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)
d <- as.matrix(
  read.csv(
    "results/v0.3_p_distance_matrix.csv",
    row.names = 1,
    check.names = FALSE
  )
)
storage.mode(d) <- "numeric"

tree <- ape::nj(as.dist(d))
ape::write.tree(tree, file = "results/v0.3_neighbor_joining_tree.nwk")

label_map <- setNames(
  paste0(manifest$virus_name, " [", manifest$current_genus, "]"),
  manifest$accession
)
tree_display <- tree
tree_display$tip.label <- unname(label_map[tree$tip.label])

png("results/v0.3_neighbor_joining_tree.png", width = 1800, height = 1800, res = 150)
plot(
  tree_display,
  type = "phylogram",
  cex = 0.65,
  main = "Exploratory Neighbor-Joining tree from whole-genome distance"
)
dev.off()
