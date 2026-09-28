# Pairwise metrics derived from a multiple DNA alignment.

pairwise_alignment_metrics <- function(aligned) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for alignment metrics.", call. = FALSE)
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

  labels <- names(aligned)
  if (is.null(labels) || any(!nzchar(labels))) {
    labels <- paste0("sequence_", seq_along(aligned))
  }

  chars <- lapply(as.character(aligned), function(x) strsplit(x, "", fixed = TRUE)[[1]])
  canonical <- c("A", "C", "G", "T")
  rows <- list()
  k <- 1L

  for (i in seq_len(length(chars) - 1L)) {
    for (j in (i + 1L):length(chars)) {
      x <- chars[[i]]
      y <- chars[[j]]

      x_letter <- !x %in% c("-", ".")
      y_letter <- !y %in% c("-", ".")
      shared <- x_letter & y_letter
      comparable <- shared & x %in% canonical & y %in% canonical

      len_x <- sum(x_letter)
      len_y <- sum(y_letter)
      shared_bases <- sum(shared)
      comparable_bases <- sum(comparable)
      matches <- sum(x[comparable] == y[comparable])

      rows[[k]] <- data.frame(
        sequence_1 = labels[i],
        sequence_2 = labels[j],
        length_1 = len_x,
        length_2 = len_y,
        shared_bases = shared_bases,
        comparable_bases = comparable_bases,
        identity_percent = if (comparable_bases == 0) NA_real_ else
          100 * matches / comparable_bases,
        coverage_shorter_percent = 100 * shared_bases / min(len_x, len_y),
        coverage_longer_percent = 100 * shared_bases / max(len_x, len_y),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}
