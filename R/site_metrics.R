# Site-by-site summaries for aligned DNA sequences.

alignment_site_metrics <- function(aligned) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for site metrics.", call. = FALSE)
  }
  if (!inherits(aligned, "DNAStringSet")) {
    stop("'aligned' must be a Biostrings DNAStringSet.", call. = FALSE)
  }
  if (length(aligned) < 2) {
    stop("At least two aligned sequences are required.", call. = FALSE)
  }
  widths <- Biostrings::width(aligned)
  if (length(unique(widths)) != 1) {
    stop("All sequences must have the same aligned width.", call. = FALSE)
  }

  chars <- do.call(rbind, lapply(as.character(aligned), function(x) {
    strsplit(x, "", fixed = TRUE)[[1]]
  }))

  canonical <- c("A", "C", "G", "T")
  n_seq <- nrow(chars)
  n_pos <- ncol(chars)

  rows <- vector("list", n_pos)

  for (pos in seq_len(n_pos)) {
    column <- chars[, pos]
    canonical_values <- column[column %in% canonical]
    n_canonical <- length(canonical_values)
    counts <- table(factor(canonical_values, levels = canonical))
    positive <- counts[counts > 0]

    allele_count <- sum(counts > 0)
    consensus_base <- if (n_canonical == 0) NA_character_ else
      names(counts)[which.max(counts)]
    consensus_frequency <- if (n_canonical == 0) NA_real_ else
      max(counts) / n_canonical

    entropy <- if (n_canonical <= 1 || length(positive) <= 1) {
      0
    } else {
      p <- positive / sum(positive)
      -sum(p * log2(p))
    }

    rows[[pos]] <- data.frame(
      alignment_position = pos,
      canonical_count = n_canonical,
      coverage_fraction = n_canonical / n_seq,
      gap_or_ambiguous_count = n_seq - n_canonical,
      gap_or_ambiguous_fraction = (n_seq - n_canonical) / n_seq,
      allele_count = allele_count,
      consensus_base = consensus_base,
      consensus_frequency = consensus_frequency,
      shannon_entropy = entropy,
      variable_site = allele_count > 1,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}
