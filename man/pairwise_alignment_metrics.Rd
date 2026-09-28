\name{pairwise_alignment_metrics}
\alias{pairwise_alignment_metrics}
\title{Compute pairwise identity and coverage from a DNA alignment}
\usage{
pairwise_alignment_metrics(aligned)
}
\arguments{
  \item{aligned}{A Biostrings \code{DNAStringSet} whose sequences all have the same aligned width.}
}
\value{
A data frame with one row per sequence pair. It reports ungapped sequence lengths, shared aligned bases, comparable canonical A/C/G/T positions, nucleotide identity among comparable positions, and coverage relative to the shorter and longer genome.
}
\description{
Separates nucleotide identity from alignment coverage so that high identity over a short shared region is not mistaken for high whole-genome similarity. Gap characters are excluded from shared-base counts. Non-canonical letters are excluded from the identity denominator.
}
