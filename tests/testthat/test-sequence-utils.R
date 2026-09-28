test_that("reverse complement handles canonical bases and N", {
  expect_equal(reverse_complement("ATGC"), "GCAT")
  expect_equal(reverse_complement("A NGC"), "GCNT")
})

test_that("sequence QC reports length, GC and ambiguous bases", {
  qc <- sequence_qc("ACGTNN")
  expect_equal(qc$length, 6)
  expect_equal(qc$gc_percent, 50)
  expect_equal(qc$n_count, 2)
  expect_equal(qc$n_percent, 100 * 2 / 6)
})

test_that("invalid DNA is rejected", {
  expect_error(normalize_dna("ACGTX"), "unsupported")
})
