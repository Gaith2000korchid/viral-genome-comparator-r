# Multiple alignment of the six-genome pilot panel.
# DECIPHER is used as the production aligner; the educational Needleman-Wunsch
# implementation in R/needleman_wunsch.R is intentionally not used here.

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings before running this analysis.")
}
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  stop("Install DECIPHER before running this analysis.")
}

manifest <- read.csv("inst/extdata/genome_manifest.csv", stringsAsFactors = FALSE)
paths <- file.path("data/raw", paste0(manifest$accession, ".fasta"))

missing <- paths[!file.exists(paths)]
if (length(missing) > 0) {
  stop(
    "Missing downloaded FASTA files: ",
    paste(missing, collapse = ", "),
    ". Run analysis/03_download_genomes.R first."
  )
}

genomes <- do.call(
  c,
  lapply(paths, Biostrings::readDNAStringSet)
)
names(genomes) <- manifest$accession

aligned <- DECIPHER::AlignSeqs(
  genomes,
  processors = 1,
  verbose = TRUE
)

dir.create("results", showWarnings = FALSE)
Biostrings::writeXStringSet(
  aligned,
  filepath = "results/pilot_alignment.fasta",
  format = "fasta"
)
