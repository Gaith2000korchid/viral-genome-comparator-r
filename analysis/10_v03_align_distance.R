# Align the expanded v0.3 panel and compute pairwise sequence metrics.

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings before running this analysis.")
}
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  stop("Install DECIPHER before running this analysis.")
}

manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)
paths <- file.path("data/raw_v03", paste0(manifest$accession, ".fasta"))

missing <- paths[!file.exists(paths)]
if (length(missing) > 0) {
  stop(
    "Missing downloaded FASTA files: ",
    paste(missing, collapse = ", "),
    ". Run analysis/08_download_v03_panel.R first."
  )
}

genomes <- do.call(c, lapply(paths, Biostrings::readDNAStringSet))
names(genomes) <- manifest$accession

aligned <- DECIPHER::AlignSeqs(
  genomes,
  processors = 2,
  verbose = TRUE
)

dir.create("results", showWarnings = FALSE)
Biostrings::writeXStringSet(
  aligned,
  filepath = "results/v0.3_alignment.fasta",
  format = "fasta"
)

metrics <- pairwise_alignment_metrics(aligned)
write.csv(
  metrics,
  "results/v0.3_pairwise_alignment_metrics.csv",
  row.names = FALSE
)

p_distance <- DECIPHER::DistanceMatrix(
  aligned,
  method = "overlap",
  includeTerminalGaps = FALSE,
  penalizeGapLetterMatches = TRUE,
  correction = NA,
  processors = 2,
  verbose = FALSE
)

write.csv(
  as.matrix(p_distance),
  "results/v0.3_p_distance_matrix.csv",
  row.names = TRUE
)

message("Aligned width: ", unique(Biostrings::width(aligned)), " columns")
