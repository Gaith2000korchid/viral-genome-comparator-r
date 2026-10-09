# Align the expanded v0.3 panel with MAFFT and compute pairwise metrics.
# MAFFT replaces DECIPHER::AlignSeqs here because the 29-genome profile
# alignment exceeds DECIPHER's maximum internal alignment size.

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Install Biostrings before running this analysis.")
}
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  stop("Install DECIPHER before running this analysis.")
}
if (!nzchar(Sys.which("mafft"))) {
  stop("MAFFT executable was not found on PATH.")
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

dir.create("results", showWarnings = FALSE)

input_path <- "results/v0.3_unaligned_panel.fasta"
alignment_path <- "results/v0.3_alignment.fasta"
mafft_log <- "results/v0.3_mafft.log"

Biostrings::writeXStringSet(
  genomes,
  filepath = input_path,
  format = "fasta"
)

status <- system2(
  "mafft",
  args = c("--auto", "--thread", "2", input_path),
  stdout = alignment_path,
  stderr = mafft_log
)

if (!identical(status, 0L)) {
  stop("MAFFT failed with exit status ", status, ". See ", mafft_log, ".")
}

aligned <- Biostrings::readDNAStringSet(alignment_path)
if (length(aligned) != nrow(manifest)) {
  stop(
    "Expected ", nrow(manifest),
    " aligned genomes but MAFFT returned ", length(aligned), "."
  )
}
if (length(unique(Biostrings::width(aligned))) != 1) {
  stop("MAFFT output is not a rectangular multiple sequence alignment.")
}
if (!setequal(names(aligned), manifest$accession)) {
  stop("MAFFT output sequence names do not match the frozen manifest.")
}

metrics <- pairwise_alignment_metrics(aligned)
write.csv(
  metrics,
  "results/v0.3_pairwise_alignment_metrics.csv",
  row.names = FALSE
)

# DECIPHER 3.2.0, pinned in environment-linux-64.explicit.txt, accepts
# correction = "none" for the uncorrected overlap distance. correction = NA
# is rejected as an invalid correction method in this release. Later manuals
# document NA as the uncorrected default; do not change this without retesting
# the frozen environment.
p_distance <- DECIPHER::DistanceMatrix(
  aligned,
  includeTerminalGaps = FALSE,
  penalizeGapLetterMatches = TRUE,
  correction = "none",
  processors = 2,
  verbose = FALSE
)

write.csv(
  as.matrix(p_distance),
  "results/v0.3_p_distance_matrix.csv",
  row.names = TRUE
)

versions <- c(
  paste("R", R.version.string),
  paste("viralGenomeComparator", as.character(utils::packageVersion("viralGenomeComparator"))),
  paste("Biostrings", as.character(utils::packageVersion("Biostrings"))),
  paste("DECIPHER", as.character(utils::packageVersion("DECIPHER"))),
  paste("MAFFT", paste(system2("mafft", "--version", stdout = TRUE, stderr = TRUE), collapse = " "))
)
writeLines(versions, "results/v0.3_software_versions.txt")

message("Aligned genomes: ", length(aligned))
message("Aligned width: ", unique(Biostrings::width(aligned)), " columns")
