test_that("fasta_qc summarizes a DNA FASTA record", {
  skip_if_not_installed("Biostrings")

  fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">toy", "ACGTNN"), fasta)

  qc <- fasta_qc(fasta)
  expect_equal(nrow(qc), 1)
  expect_equal(qc$length, 6)
  expect_equal(qc$gc_percent, 50)
  expect_equal(qc$ambiguous_count, 2)
  expect_equal(qc$ambiguous_percent, 100 * 2 / 6)
})
