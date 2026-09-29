test_that("adjusted_rand_index is one for identical partitions", {
  truth <- c("A", "A", "B", "B")
  predicted <- c(2, 2, 1, 1)
  expect_equal(adjusted_rand_index(truth, predicted), 1)
})

test_that("adjusted_rand_index is invariant to cluster label names", {
  truth <- c("A", "A", "B", "B", "C", "C")
  x <- c("x", "x", "y", "y", "z", "z")
  y <- c("p", "p", "r", "r", "q", "q")
  expect_equal(adjusted_rand_index(truth, x), adjusted_rand_index(truth, y))
})

test_that("adjusted_rand_index validates lengths", {
  expect_error(adjusted_rand_index(c("A", "B"), "A"), "same length")
})
