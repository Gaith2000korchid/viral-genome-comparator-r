# Orchestrate the six-genome pilot. Load the package before sourcing scripts
# that call its functions.
library(viralGenomeComparator)

steps <- c(
  "03_download_genomes.R",
  "04_genome_qc.R",
  "05_align_genomes.R",
  "06_pairwise_metrics.R",
  "07_cluster_heatmap.R"
)
for (step in steps) {
  message("Running ", step)
  source(file.path("analysis", step))
}
