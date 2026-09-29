# Quantify how well sequence-derived distances recover the two genus labels.

manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)
distance_path <- "results/v0.3_p_distance_matrix.csv"

if (!file.exists(distance_path)) {
  stop("Missing distance matrix. Run analysis/10_v03_align_distance.R first.")
}

d <- as.matrix(read.csv(distance_path, row.names = 1, check.names = FALSE))
storage.mode(d) <- "numeric"

if (!identical(rownames(d), colnames(d))) {
  stop("Distance matrix row and column names must match.")
}

manifest <- manifest[match(rownames(d), manifest$accession), ]
if (anyNA(manifest$current_genus)) {
  stop("Distance-matrix accession missing from v0.3 manifest.")
}

groups <- manifest$current_genus
pair_mask <- upper.tri(d)
same_genus <- outer(groups, groups, "==")

within <- d[pair_mask & same_genus]
between <- d[pair_mask & !same_genus]

hc <- hclust(as.dist(d), method = "average")
clusters <- cutree(hc, k = length(unique(groups)))
ari <- adjusted_rand_index(groups, clusters)

d_nn <- d
diag(d_nn) <- Inf
nearest_index <- apply(d_nn, 1, which.min)
nearest_accession <- colnames(d_nn)[nearest_index]
nearest_group <- groups[nearest_index]
nn_same <- nearest_group == groups

set.seed(20260929)
n_perm <- 999L
observed_delta <- mean(between) - mean(within)
perm_delta <- numeric(n_perm)

for (b in seq_len(n_perm)) {
  permuted <- sample(groups, replace = FALSE)
  same_perm <- outer(permuted, permuted, "==")
  within_perm <- d[pair_mask & same_perm]
  between_perm <- d[pair_mask & !same_perm]
  perm_delta[b] <- mean(between_perm) - mean(within_perm)
}

permutation_p <- (1 + sum(perm_delta >= observed_delta)) / (n_perm + 1)

summary <- data.frame(
  n_genomes = nrow(d),
  n_teseptimavirus = sum(groups == "Teseptimavirus"),
  n_teetrevirus = sum(groups == "Teetrevirus"),
  mean_within_genus_distance = mean(within),
  mean_between_genus_distance = mean(between),
  between_minus_within = observed_delta,
  nearest_neighbor_genus_accuracy = mean(nn_same),
  adjusted_rand_index_k2 = ari,
  permutation_p_value = permutation_p,
  permutations = n_perm,
  stringsAsFactors = FALSE
)

nearest <- data.frame(
  accession = rownames(d),
  virus_name = manifest$virus_name,
  genus = groups,
  nearest_accession = nearest_accession,
  nearest_virus_name = manifest$virus_name[match(nearest_accession, manifest$accession)],
  nearest_genus = nearest_group,
  nearest_distance = d[cbind(seq_len(nrow(d)), nearest_index)],
  same_genus = nn_same,
  stringsAsFactors = FALSE
)

dir.create("results", showWarnings = FALSE)
write.csv(summary, "results/v0.3_taxonomy_concordance.csv", row.names = FALSE)
write.csv(nearest, "results/v0.3_nearest_neighbors.csv", row.names = FALSE)

png("results/v0.3_within_between_distance.png", width = 1400, height = 1000, res = 150)
boxplot(
  list(
    "Within genus" = within,
    "Between genera" = between
  ),
  ylab = "Uncorrected DECIPHER distance",
  main = "Within-genus versus between-genus whole-genome distance"
)
dev.off()

print(summary)
