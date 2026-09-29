test_that("alignment_site_summary classifies conserved and variable sites", {
  skip_if_not_installed("Biostrings")

  aligned <- Biostrings::DNAStringSet(c(
    a = "ACGT",
    b = "ACGT",
    c = "ATGT",
    d = "A-GT"
  ))

  summary <- alignment_site_summary(
    aligned,
    min_occupancy = 0.75,
    conserved_threshold = 0.95
  )

  expect_equal(nrow(summary), 4)
  expect_equal(summary$site_class[1], "conserved")
  expect_equal(summary$site_class[2], "variable")
  expect_equal(summary$occupancy_percent[2], 75)
  expect_true(summary$entropy_bits[2] > 0)
})

test_that("alignment_site_summary rejects unaligned sequences", {
  skip_if_not_installed("Biostrings")

  x <- Biostrings::DNAStringSet(c("ACGT", "ACG"))
  expect_error(alignment_site_summary(x), "same aligned width")
})
