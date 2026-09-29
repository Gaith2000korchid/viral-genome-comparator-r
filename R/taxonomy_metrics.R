# Agreement metrics for comparing sequence-derived clusters with metadata groups.

adjusted_rand_index <- function(truth, predicted) {
  if (length(truth) != length(predicted)) {
    stop("'truth' and 'predicted' must have the same length.", call. = FALSE)
  }
  if (length(truth) < 2) {
    stop("At least two observations are required.", call. = FALSE)
  }
  if (anyNA(truth) || anyNA(predicted)) {
    stop("Cluster labels must not contain missing values.", call. = FALSE)
  }

  tab <- table(as.character(truth), as.character(predicted))
  choose2 <- function(x) x * (x - 1) / 2

  index <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(sum(tab))

  expected <- row_pairs * col_pairs / total_pairs
  max_index <- (row_pairs + col_pairs) / 2
  denominator <- max_index - expected

  if (denominator == 0) {
    return(if (index == max_index) 1 else 0)
  }

  (index - expected) / denominator
}
