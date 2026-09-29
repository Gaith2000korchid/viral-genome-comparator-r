# Compare distance-derived clustering with curated genus labels.

matrix_path <- "results/p_distance_matrix.csv"
manifest_path <- "inst/extdata/genome_manifest.csv"

if (!file.exists(matrix_path)) {
  stop("Missing distance matrix. Run analysis/06_pairwise_metrics.R first.")
}

distance_matrix <- as.matrix(
  read.csv(matrix_path, row.names = 1, check.names = FALSE)
)
storage.mode(distance_matrix) <- "numeric"

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
manifest <- manifest[match(rownames(distance_matrix), manifest$accession), ]

if (anyNA(manifest$accession)) {
  stop("Distance-matrix accessions are missing from the manifest.")
}

hc <- hclust(as.dist(distance_matrix), method = "average")
cluster_labels <- cutree(hc, k = length(unique(manifest$current_genus)))
ari <- adjusted_rand_index(cluster_labels, manifest$current_genus)

comparison <- data.frame(
  accession = manifest$accession,
  virus_name = manifest$virus_name,
  current_genus = manifest$current_genus,
  distance_cluster = unname(cluster_labels),
  stringsAsFactors = FALSE
)

summary <- data.frame(
  n_genomes = nrow(comparison),
  n_taxonomic_genera = length(unique(manifest$current_genus)),
  adjusted_rand_index = ari,
  exact_genus_recovery = ari == 1,
  stringsAsFactors = FALSE
)

dir.create("results", showWarnings = FALSE)
write.csv(comparison, "results/taxonomy_cluster_comparison.csv", row.names = FALSE)
write.csv(summary, "results/taxonomy_concordance_summary.csv", row.names = FALSE)

print(comparison)
print(summary)
