# Diagnose the apparent disagreement between local genus structure and
# a two-cluster hierarchical cut in the v0.3 distance matrix.

manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)
d <- as.matrix(
  read.csv(
    "results/v0.3_p_distance_matrix.csv",
    row.names = 1,
    check.names = FALSE
  )
)
storage.mode(d) <- "numeric"

manifest <- manifest[match(rownames(d), manifest$accession), ]
if (anyNA(manifest$current_genus)) {
  stop("Distance-matrix accession missing from v0.3 manifest.")
}

groups <- manifest$current_genus
hc <- hclust(as.dist(d), method = "average")
cluster_k2 <- cutree(hc, k = 2)

assignments <- data.frame(
  accession = rownames(d),
  virus_name = manifest$virus_name,
  current_genus = groups,
  cluster_k2 = unname(cluster_k2),
  stringsAsFactors = FALSE
)

pair_rows <- list()
k <- 1L
for (i in seq_len(nrow(d) - 1L)) {
  for (j in (i + 1L):nrow(d)) {
    pair_rows[[k]] <- data.frame(
      genus_1 = groups[i],
      genus_2 = groups[j],
      same_genus = groups[i] == groups[j],
      distance = d[i, j],
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
pairs <- do.call(rbind, pair_rows)

within_by_genus <- do.call(
  rbind,
  lapply(sort(unique(groups)), function(g) {
    values <- pairs$distance[
      pairs$same_genus & pairs$genus_1 == g & pairs$genus_2 == g
    ]
    data.frame(
      comparison = paste0("within_", g),
      n_pairs = length(values),
      mean_distance = mean(values),
      median_distance = median(values),
      min_distance = min(values),
      max_distance = max(values),
      stringsAsFactors = FALSE
    )
  })
)

between <- pairs$distance[!pairs$same_genus]
distance_summary <- rbind(
  within_by_genus,
  data.frame(
    comparison = "between_genera",
    n_pairs = length(between),
    mean_distance = mean(between),
    median_distance = median(between),
    min_distance = min(between),
    max_distance = max(between),
    stringsAsFactors = FALSE
  )
)

cluster_sizes <- as.data.frame(table(cluster_k2), stringsAsFactors = FALSE)
names(cluster_sizes) <- c("cluster_k2", "n_genomes")

dir.create("results", showWarnings = FALSE)
write.csv(assignments, "results/v0.3_cluster_assignments_k2.csv", row.names = FALSE)
write.csv(distance_summary, "results/v0.3_genus_distance_summary.csv", row.names = FALSE)
write.csv(cluster_sizes, "results/v0.3_cluster_sizes_k2.csv", row.names = FALSE)

print(assignments[order(assignments$cluster_k2, assignments$current_genus), ])
print(distance_summary)
print(cluster_sizes)
