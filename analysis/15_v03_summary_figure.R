# Regenerate the compact figure from a versioned numerical summary.
summary <- read.csv("results/v0.3_taxonomy_concordance.csv")
svg("results/v0.3_taxonomy_summary.svg", width = 7, height = 4.5)
par(mar = c(5, 4, 4, 2) + 0.1)
barplot(
  c("Within genus" = summary$mean_within_genus_distance, "Between genera" = summary$mean_between_genus_distance),
  col = c("#2563eb", "#64748b"), ylim = c(0, 0.6),
  ylab = "Mean uncorrected DECIPHER distance",
  main = "29-genome panel: local signal, imperfect global clustering"
)
mtext(sprintf("Nearest-neighbor genus agreement: %.0f%% | forced k=2 ARI: %.3f",
              100 * summary$nearest_neighbor_genus_accuracy, summary$adjusted_rand_index_k2), side = 1, line = 3.5, cex = 0.8)
dev.off()
