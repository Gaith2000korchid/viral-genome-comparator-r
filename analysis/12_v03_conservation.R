# Site-level and windowed conservation analysis of the v0.3 alignment.

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings before running this analysis.")
}

alignment_path <- "results/v0.3_alignment.fasta"
if (!file.exists(alignment_path)) {
  stop("Missing v0.3 alignment. Run analysis/10_v03_align_distance.R first.")
}

aligned <- Biostrings::readDNAStringSet(alignment_path)
manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)
groups <- manifest$current_genus[match(names(aligned), manifest$accession)]

if (anyNA(groups)) {
  stop("Alignment accession missing from v0.3 manifest.")
}

summarize_windows <- function(site_summary, label, window_size = 500L) {
  window <- ceiling(site_summary$aligned_position / window_size)
  pieces <- split(site_summary, window)

  do.call(rbind, lapply(seq_along(pieces), function(i) {
    x <- pieces[[i]]
    high <- x$occupancy_percent >= 80

    data.frame(
      panel = label,
      window = i,
      start_aligned_position = min(x$aligned_position),
      end_aligned_position = max(x$aligned_position),
      mean_occupancy_percent = mean(x$occupancy_percent),
      mean_entropy_high_occupancy = if (any(high)) mean(x$entropy_bits[high]) else NA_real_,
      conserved_fraction_high_occupancy = if (any(high)) mean(x$site_class[high] == "conserved") else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

panel_summaries <- list()
window_summaries <- list()

analyze_panel <- function(x, label) {
  site <- alignment_site_summary(
    x,
    min_occupancy = 0.8,
    conserved_threshold = 0.95
  )

  high <- site$occupancy_percent >= 80
  summary <- data.frame(
    panel = label,
    n_sequences = length(x),
    aligned_width = nrow(site),
    high_occupancy_sites = sum(high),
    conserved_sites = sum(site$site_class == "conserved"),
    variable_sites = sum(site$site_class == "variable"),
    low_occupancy_sites = sum(site$site_class == "low_occupancy"),
    mean_entropy_high_occupancy = if (any(high)) mean(site$entropy_bits[high]) else NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    site = site,
    summary = summary,
    windows = summarize_windows(site, label)
  )
}

all_result <- analyze_panel(aligned, "all")
write.csv(
  all_result$site,
  "results/v0.3_site_conservation.csv",
  row.names = FALSE
)
panel_summaries[[1]] <- all_result$summary
window_summaries[[1]] <- all_result$windows

for (g in sort(unique(groups))) {
  result <- analyze_panel(aligned[groups == g], g)
  panel_summaries[[length(panel_summaries) + 1]] <- result$summary
  window_summaries[[length(window_summaries) + 1]] <- result$windows
}

summary_table <- do.call(rbind, panel_summaries)
window_table <- do.call(rbind, window_summaries)

write.csv(
  summary_table,
  "results/v0.3_conservation_summary.csv",
  row.names = FALSE
)
write.csv(
  window_table,
  "results/v0.3_conservation_windows.csv",
  row.names = FALSE
)

panels <- unique(window_table$panel)
line_types <- seq_along(panels)

png("results/v0.3_entropy_windows.png", width = 1800, height = 1000, res = 150)
plot(
  NA,
  xlim = range(window_table$start_aligned_position),
  ylim = range(window_table$mean_entropy_high_occupancy, na.rm = TRUE),
  xlab = "Aligned genome position",
  ylab = "Mean Shannon entropy (bits)",
  main = "Windowed nucleotide diversity across the v0.3 alignment"
)
for (i in seq_along(panels)) {
  x <- window_table[window_table$panel == panels[i], ]
  lines(
    x$start_aligned_position,
    x$mean_entropy_high_occupancy,
    lty = line_types[i],
    lwd = 2
  )
}
legend(
  "topright",
  legend = panels,
  lty = line_types,
  lwd = 2,
  bty = "n"
)
dev.off()

print(summary_table)
