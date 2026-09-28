test_that("pairwise alignment metrics separate identity from coverage", {
  skip_if_not_installed("Biostrings")

  aligned <- Biostrings::DNAStringSet(c(
    genome_A = "ACGT-",
    genome_B = "AC-TT"
  ))

  metrics <- pairwise_alignment_metrics(aligned)

  expect_equal(nrow(metrics), 1)
  expect_equal(metrics$shared_bases, 3)
  expect_equal(metrics$comparable_bases, 3)
  expect_equal(metrics$identity_percent, 100)
  expect_equal(metrics$coverage_shorter_percent, 75)
  expect_equal(metrics$coverage_longer_percent, 75)
})

test_that("pairwise alignment metrics reject unaligned widths", {
  skip_if_not_installed("Biostrings")

  unaligned <- Biostrings::DNAStringSet(c("ACGT", "ACG"))
  expect_error(pairwise_alignment_metrics(unaligned), "same aligned width")
})
