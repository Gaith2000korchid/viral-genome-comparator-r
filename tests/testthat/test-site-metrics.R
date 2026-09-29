test_that("alignment_site_summary distinguishes conservation, variation and gaps", {
  skip_if_not_installed("Biostrings")

  aligned <- Biostrings::DNAStringSet(c(
    a = "ACGT",
    b = "AC-T",
    c = "ATGT"
  ))

  metrics <- alignment_site_summary(
    aligned,
    min_occupancy = 0.8,
    conserved_threshold = 0.95
  )

  expect_equal(nrow(metrics), 4)
  expect_equal(metrics$occupancy_percent[1], 100)
  expect_equal(metrics$major_frequency_percent[1], 100)
  expect_equal(metrics$site_class[1], "conserved")

  expect_equal(metrics$distinct_canonical_bases[2], 2)
  expect_equal(metrics$site_class[2], "variable")

  expect_equal(metrics$occupancy_percent[3], 100 * 2 / 3)
  expect_equal(metrics$site_class[3], "low_occupancy")
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
