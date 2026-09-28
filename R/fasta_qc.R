# FASTA-level QC helpers for complete DNA genomes.

fasta_qc <- function(fasta_path) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for FASTA QC.", call. = FALSE)
  }

  genomes <- Biostrings::readDNAStringSet(fasta_path, format = "fasta")
  if (length(genomes) == 0) {
    stop("FASTA file contains no sequences.", call. = FALSE)
  }

  alphabet <- Biostrings::alphabetFrequency(genomes, baseOnly = FALSE, collapse = FALSE)
  canonical <- rowSums(alphabet[, c("A", "C", "G", "T"), drop = FALSE])
  gc <- rowSums(alphabet[, c("G", "C"), drop = FALSE])
  total <- Biostrings::width(genomes)
  ambiguous <- total - canonical

  data.frame(
    record = names(genomes),
    length = total,
    gc_percent = ifelse(canonical == 0, NA_real_, 100 * gc / canonical),
    ambiguous_count = ambiguous,
    ambiguous_percent = 100 * ambiguous / total,
    stringsAsFactors = FALSE
  )
}
