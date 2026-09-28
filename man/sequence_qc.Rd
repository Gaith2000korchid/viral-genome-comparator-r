\name{sequence_qc}
\alias{sequence_qc}
\title{Compute basic DNA sequence quality-control statistics}
\usage{
sequence_qc(sequence)
}
\arguments{
  \item{sequence}{A DNA character string containing A, C, G, T, or N.}
}
\value{
A one-row data frame containing sequence length, GC percentage among canonical bases, N count, and N percentage.
}
\description{
Calculates simple quality-control statistics for a DNA sequence after normalization and validation.
}
\examples{
sequence_qc("ACGTNN")
}
