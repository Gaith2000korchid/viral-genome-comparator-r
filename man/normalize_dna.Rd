\name{normalize_dna}
\alias{normalize_dna}
\title{Normalize and validate a DNA sequence}
\usage{
normalize_dna(sequence)
}
\arguments{
  \item{sequence}{A character string containing DNA bases. Whitespace is removed and letters are converted to uppercase.}
}
\value{A normalized uppercase DNA character string.}
\description{
Removes whitespace, converts the input to uppercase, and validates that only A, C, G, T, and N are present.
}
\examples{
normalize_dna("ac gt n")
}
