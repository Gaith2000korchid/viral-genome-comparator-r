# Transparent Needleman-Wunsch implementation for short teaching examples.
# This is intentionally not intended for whole viral genomes.

needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -2) {
  seq1 <- normalize_dna(seq1)
  seq2 <- normalize_dna(seq2)
  x <- strsplit(seq1, "", fixed = TRUE)[[1]]
  y <- strsplit(seq2, "", fixed = TRUE)[[1]]
  n <- length(x)
  m <- length(y)

  score <- matrix(0, nrow = n + 1, ncol = m + 1)
  score[, 1] <- (0:n) * gap
  score[1, ] <- (0:m) * gap

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      diagonal <- score[i, j] + if (x[i] == y[j]) match else mismatch
      up <- score[i, j + 1] + gap
      left <- score[i + 1, j] + gap
      score[i + 1, j + 1] <- max(diagonal, up, left)
    }
  }

  aligned1 <- character()
  aligned2 <- character()
  i <- n
  j <- m

  while (i > 0 || j > 0) {
    if (i > 0 && j > 0) {
      diagonal_score <- score[i, j] + if (x[i] == y[j]) match else mismatch
      if (score[i + 1, j + 1] == diagonal_score) {
        aligned1 <- c(x[i], aligned1)
        aligned2 <- c(y[j], aligned2)
        i <- i - 1
        j <- j - 1
        next
      }
    }

    if (i > 0 && score[i + 1, j + 1] == score[i, j + 1] + gap) {
      aligned1 <- c(x[i], aligned1)
      aligned2 <- c("-", aligned2)
      i <- i - 1
    } else {
      aligned1 <- c("-", aligned1)
      aligned2 <- c(y[j], aligned2)
      j <- j - 1
    }
  }

  list(
    aligned_seq1 = paste(aligned1, collapse = ""),
    aligned_seq2 = paste(aligned2, collapse = ""),
    score = score[n + 1, m + 1],
    score_matrix = score
  )
}
