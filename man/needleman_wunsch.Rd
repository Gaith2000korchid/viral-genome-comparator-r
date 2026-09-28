\name{needleman_wunsch}
\alias{needleman_wunsch}
\title{Globally align two short DNA sequences with Needleman-Wunsch}
\usage{
needleman_wunsch(seq1, seq2, match = 1, mismatch = -1, gap = -2)
}
\arguments{
  \item{seq1}{First DNA character string.}
  \item{seq2}{Second DNA character string.}
  \item{match}{Score assigned to a nucleotide match.}
  \item{mismatch}{Score assigned to a nucleotide mismatch.}
  \item{gap}{Linear gap penalty.}
}
\value{
A list containing the two aligned sequences, the optimal global-alignment score, and the dynamic-programming score matrix.
}
\description{
Implements the Needleman-Wunsch dynamic-programming algorithm for educational use. The implementation uses O(n*m) time and memory and is not intended as a production whole-genome aligner.
}
\examples{
needleman_wunsch("GATTACA", "GCATGCT")
}
