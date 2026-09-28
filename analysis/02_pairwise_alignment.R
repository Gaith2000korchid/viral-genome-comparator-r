# Needleman-Wunsch demonstration on short sequences.
# The educational implementation is O(n*m) in time and memory and therefore
# must not be used as the production whole-genome aligner.

example <- needleman_wunsch(
  "GATTACA",
  "GCATGCU"
)
# The second sequence deliberately contains U and will be rejected because
# v0.1 currently models DNA only. Replace U with T to run the example:
example <- needleman_wunsch("GATTACA", "GCATGCT")
example$aligned_seq1
example$aligned_seq2
example$score
