# Site-level conservation summaries from aligned DNA sequences.

alignment_site_summary <- function(aligned,
                                   min_occupancy = 0.8,
                                   conserved_threshold = 0.95) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for conservation analysis.", call. = FALSE)
  }
  if (!inherits(aligned, "DNAStringSet")) {
    stop("'aligned' must be a Biostrings DNAStringSet.", call. = FALSE)
  }
  if (length(aligned) < 2) {
    stop("At least two aligned sequences are required.", call. = FALSE)
  }
  if (length(unique(Biostrings::width(aligned))) != 1) {
    stop("All sequences must have the same aligned width.", call. = FALSE)
  }
  if (!is.numeric(min_occupancy) || length(min_occupancy) != 1 ||
      min_occupancy < 0 || min_occupancy > 1) {
    stop("'min_occupancy' must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(conserved_threshold) || length(conserved_threshold) != 1 ||
      conserved_threshold < 0 || conserved_threshold > 1) {
    stop("'conserved_threshold' must be a number between 0 and 1.", call. = FALSE)
  }

  raw_counts <- Biostrings::consensusMatrix(aligned, as.prob = FALSE)
  bases <- c("A", "C", "G", "T")

  counts <- matrix(
    0,
    nrow = length(bases),
    ncol = ncol(raw_counts),
    dimnames = list(bases, colnames(raw_counts))
  )
  present <- intersect(bases, rownames(raw_counts))
  counts[present, ] <- raw_counts[present, , drop = FALSE]

  canonical_total <- colSums(counts)
  occupancy <- canonical_total / length(aligned)
  major_count <- apply(counts, 2, max)
  major_index <- max.col(t(counts), ties.method = "first")
  major_base <- bases[major_index]
  major_base[canonical_total == 0] <- NA_character_

  major_frequency <- ifelse(
    canonical_total == 0,
    NA_real_,
    major_count / canonical_total
  )

  denom <- pmax(canonical_total, 1)
  probabilities <- sweep(counts, 2, denom, "/")
  entropy_bits <- -colSums(
    ifelse(probabilities > 0, probabilities * log2(probabilities), 0)
  )

  distinct_canonical_bases <- colSums(counts > 0)

  site_class <- ifelse(
    occupancy < min_occupancy,
    "low_occupancy",
    ifelse(
      major_frequency >= conserved_threshold,
      "conserved",
      "variable"
    )
  )

  data.frame(
    aligned_position = seq_len(ncol(counts)),
    occupancy_percent = 100 * occupancy,
    major_base = major_base,
    major_frequency_percent = 100 * major_frequency,
    entropy_bits = entropy_bits,
    distinct_canonical_bases = distinct_canonical_bases,
    site_class = site_class,
    stringsAsFactors = FALSE
  )
}
