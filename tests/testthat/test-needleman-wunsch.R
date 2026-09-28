test_that("Needleman-Wunsch aligns identical sequences", {
  result <- needleman_wunsch("ACGT", "ACGT")
  expect_equal(result$aligned_seq1, "ACGT")
  expect_equal(result$aligned_seq2, "ACGT")
  expect_equal(result$score, 4)
})

test_that("Needleman-Wunsch returns a valid global alignment with a gap", {
  result <- needleman_wunsch("ACGT", "AGT")
  expect_equal(gsub("-", "", result$aligned_seq1), "ACGT")
  expect_equal(gsub("-", "", result$aligned_seq2), "AGT")
  expect_equal(nchar(result$aligned_seq1), nchar(result$aligned_seq2))
  expect_equal(result$score, 1)
})
