# Adjusted Rand Index for two categorical partitions.

adjusted_rand_index <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b)) {
    stop("Label vectors must have the same length.", call. = FALSE)
  }
  if (length(labels_a) < 2) {
    stop("At least two observations are required.", call. = FALSE)
  }
  if (anyNA(labels_a) || anyNA(labels_b)) {
    stop("Labels must not contain missing values.", call. = FALSE)
  }

  tab <- table(labels_a, labels_b)
  choose2 <- function(x) x * (x - 1) / 2

  index <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(sum(tab))

  expected <- row_pairs * col_pairs / total_pairs
  maximum <- (row_pairs + col_pairs) / 2
  denominator <- maximum - expected

  if (denominator == 0) {
    return(if (index == maximum) 1 else 0)
  }

  (index - expected) / denominator
}
