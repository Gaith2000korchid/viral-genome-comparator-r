# Basic nucleotide-sequence utilities used in the v0.1 learning layer.

normalize_dna <- function(sequence) {
  sequence <- toupper(gsub("\\s+", "", sequence))
  if (!nzchar(sequence)) {
    stop("DNA sequence must not be empty.", call. = FALSE)
  }
  if (grepl("[^ACGTN]", sequence)) {
    stop("DNA sequence contains unsupported characters; allowed: A, C, G, T, N.", call. = FALSE)
  }
  sequence
}

reverse_complement <- function(sequence) {
  sequence <- normalize_dna(sequence)
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  complement <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  paste(rev(unname(complement[bases])), collapse = "")
}

sequence_qc <- function(sequence) {
  sequence <- normalize_dna(sequence)
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  counts <- table(factor(bases, levels = c("A", "C", "G", "T", "N")))
  canonical_length <- sum(counts[c("A", "C", "G", "T")])

  data.frame(
    length = nchar(sequence),
    gc_percent = if (canonical_length == 0) NA_real_ else
      100 * sum(counts[c("G", "C")]) / canonical_length,
    n_count = unname(counts["N"]),
    n_percent = 100 * unname(counts["N"]) / nchar(sequence),
    stringsAsFactors = FALSE
  )
}
