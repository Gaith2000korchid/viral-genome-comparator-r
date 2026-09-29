test_that("alignment site metrics distinguish conservation, variation and gaps", {
  skip_if_not_installed("Biostrings")

  aligned <- Biostrings::DNAStringSet(c(
    a = "ACGT",
    b = "AC-T",
    c = "ATGT"
  ))

  metrics <- alignment_site_metrics(aligned)

  expect_equal(nrow(metrics), 4)
  expect_equal(metrics$coverage_fraction[1], 1)
  expect_equal(metrics$consensus_frequency[1], 1)
  expect_false(metrics$variable_site[1])

  expect_equal(metrics$allele_count[2], 2)
  expect_true(metrics$variable_site[2])

  expect_equal(metrics$coverage_fraction[3], 2 / 3)
  expect_equal(metrics$gap_or_ambiguous_count[3], 1)
})

test_that("adjusted Rand index recognizes identical partitions", {
  expect_equal(
    adjusted_rand_index(c("A", "A", "B", "B"), c("x", "x", "y", "y")),
    1
  )
})

test_that("adjusted Rand index validates input", {
  expect_error(adjusted_rand_index(c("A", "B"), "A"), "same length")
})
