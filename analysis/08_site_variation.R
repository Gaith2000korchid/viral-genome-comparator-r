# Site-by-site variation and conservation across the pilot alignment.

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings before running this analysis.")
}

alignment_path <- "results/pilot_alignment.fasta"
if (!file.exists(alignment_path)) {
  stop("Missing alignment. Run analysis/05_align_genomes.R first.")
}

aligned <- Biostrings::readDNAStringSet(alignment_path)
site_metrics <- alignment_site_metrics(aligned)

dir.create("results", showWarnings = FALSE)
write.csv(site_metrics, "results/site_metrics.csv", row.names = FALSE)

summary_table <- data.frame(
  alignment_columns = nrow(site_metrics),
  fully_covered_columns = sum(site_metrics$coverage_fraction == 1),
  variable_columns = sum(site_metrics$variable_site & site_metrics$coverage_fraction == 1),
  invariant_columns = sum(!site_metrics$variable_site & site_metrics$coverage_fraction == 1),
  mean_consensus_frequency_full_coverage = mean(
    site_metrics$consensus_frequency[site_metrics$coverage_fraction == 1],
    na.rm = TRUE
  ),
  mean_entropy_full_coverage = mean(
    site_metrics$shannon_entropy[site_metrics$coverage_fraction == 1],
    na.rm = TRUE
  )
)

write.csv(summary_table, "results/site_metrics_summary.csv", row.names = FALSE)

png("results/conservation_profile.png", width = 1800, height = 900, res = 150)
plot(
  site_metrics$alignment_position,
  site_metrics$consensus_frequency,
  type = "l",
  xlab = "Alignment position",
  ylab = "Consensus frequency",
  main = "Genome-wide conservation profile across the six-phage pilot",
  ylim = c(0, 1)
)
abline(h = 0.9, lty = 2)
dev.off()

png("results/entropy_profile.png", width = 1800, height = 900, res = 150)
plot(
  site_metrics$alignment_position,
  site_metrics$shannon_entropy,
  type = "l",
  xlab = "Alignment position",
  ylab = "Shannon entropy (bits)",
  main = "Genome-wide site variability across the six-phage pilot"
)
dev.off()

print(summary_table)
