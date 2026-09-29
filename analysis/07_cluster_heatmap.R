# Visualize the v0.2 uncorrected p-distance matrix.
# This script uses base R only so the visualization remains lightweight.

matrix_path <- "results/p_distance_matrix.csv"
manifest_path <- "inst/extdata/genome_manifest.csv"

if (!file.exists(matrix_path)) {
  stop("Missing distance matrix. Run analysis/06_pairwise_metrics.R first.")
}

distance_matrix <- as.matrix(
  read.csv(matrix_path, row.names = 1, check.names = FALSE)
)
storage.mode(distance_matrix) <- "numeric"

if (nrow(distance_matrix) != ncol(distance_matrix)) {
  stop("Distance matrix must be square.")
}
if (!identical(rownames(distance_matrix), colnames(distance_matrix))) {
  stop("Distance matrix row and column names must match.")
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
label_map <- setNames(
  paste0(manifest$virus_name, " [", manifest$current_genus, "]"),
  manifest$accession
)

missing_labels <- setdiff(rownames(distance_matrix), names(label_map))
if (length(missing_labels) > 0) {
  stop("Accessions missing from manifest: ", paste(missing_labels, collapse = ", "))
}

hc <- hclust(as.dist(distance_matrix), method = "average")
order_idx <- hc$order
ordered <- distance_matrix[order_idx, order_idx, drop = FALSE]
pretty_labels <- unname(label_map[rownames(ordered)])

dir.create("results", showWarnings = FALSE)

write.csv(
  data.frame(
    cluster_order = seq_along(pretty_labels),
    accession = rownames(ordered),
    label = pretty_labels,
    stringsAsFactors = FALSE
  ),
  "results/cluster_order.csv",
  row.names = FALSE
)

png("results/p_distance_dendrogram.png", width = 1600, height = 1000, res = 150)
plot(
  hc,
  labels = unname(label_map[hc$labels]),
  main = "Pilot bacteriophage clustering from whole-genome p-distance",
  xlab = "",
  sub = "Average-linkage hierarchical clustering"
)
dev.off()

png("results/p_distance_heatmap.png", width = 1800, height = 1600, res = 150)
par(mar = c(13, 13, 4, 5))

n <- nrow(ordered)
image(
  x = seq_len(n),
  y = seq_len(n),
  z = t(ordered[n:1, , drop = FALSE]),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Whole-genome uncorrected p-distance"
)

axis(1, at = seq_len(n), labels = pretty_labels, las = 2, cex.axis = 0.75)
axis(2, at = seq_len(n), labels = rev(pretty_labels), las = 2, cex.axis = 0.75)

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    text(
      x = j,
      y = n - i + 1,
      labels = sprintf("%.3f", ordered[i, j]),
      cex = 0.75
    )
  }
}

box()
dev.off()

message("Wrote:")
message("  results/cluster_order.csv")
message("  results/p_distance_dendrogram.png")
message("  results/p_distance_heatmap.png")
