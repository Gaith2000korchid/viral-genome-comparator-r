# Quantify pairwise identity, coverage and uncorrected p-distance
# from the v0.2 multiple alignment.

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings before running this analysis.")
}
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  stop("Install DECIPHER before running this analysis.")
}

alignment_path <- "results/pilot_alignment.fasta"
if (!file.exists(alignment_path)) {
  stop("Missing alignment. Run analysis/05_align_genomes.R first.")
}

aligned <- Biostrings::readDNAStringSet(alignment_path)
metrics <- pairwise_alignment_metrics(aligned)

p_distance <- DECIPHER::DistanceMatrix(
  aligned,
  method = "overlap",
  includeTerminalGaps = FALSE,
  penalizeGapLetterMatches = TRUE,
  correction = "none",
  processors = 1,
  verbose = FALSE
)

write.csv(metrics, "results/pairwise_alignment_metrics.csv", row.names = FALSE)
write.csv(
  as.matrix(p_distance),
  "results/p_distance_matrix.csv",
  row.names = TRUE
)

print(metrics)
print(round(as.matrix(p_distance), 4))
