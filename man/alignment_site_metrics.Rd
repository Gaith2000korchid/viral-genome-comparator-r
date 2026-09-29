\name{alignment_site_metrics}
\alias{alignment_site_metrics}
\title{Summarize variation and conservation across a DNA alignment}
\usage{
alignment_site_metrics(aligned)
}
\arguments{
  \item{aligned}{A Biostrings \code{DNAStringSet} whose sequences all have the same aligned width.}
}
\value{
A data frame with one row per alignment column, including canonical-base coverage, gap or ambiguity fraction, allele count, consensus base, consensus frequency, Shannon entropy, and a variable-site flag.
}
\description{
Computes transparent site-level summaries from a multiple DNA alignment. Conservation is represented by consensus frequency among canonical A/C/G/T observations, while coverage is reported separately so poorly covered columns are not mistaken for conserved sites.
}
