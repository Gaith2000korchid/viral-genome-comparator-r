# Run from the repository root inside the versioned environment.
library(viralGenomeComparator)

steps <- c(
  "08_download_v03_panel.R", "09_v03_genome_qc.R",
  "10_v03_align_distance.R", "11_v03_taxonomy_concordance.R",
  "12_v03_conservation.R", "13_v03_neighbor_joining.R",
  "14_v03_cluster_diagnostics.R", "15_v03_summary_figure.R"
)

dir.create("results", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "results/v0.3_session_info.txt")
for (step in steps) {
  message("Running ", step)
  source(file.path("analysis", step))
}
writeLines(
  capture.output(sessionInfo()),
  "results/v0.3_session_info.txt"
)
