\name{reverse_complement}
\alias{reverse_complement}
\title{Compute the reverse complement of a DNA sequence}
\usage{
reverse_complement(sequence)
}
\arguments{
  \item{sequence}{A DNA character string containing A, C, G, T, or N.}
}
\value{A character string containing the reverse-complement DNA sequence.}
\description{
Validates a DNA sequence, complements each nucleotide, and reverses the resulting sequence. Ambiguous N bases remain N.
}
\examples{
reverse_complement("ATGC")
}
