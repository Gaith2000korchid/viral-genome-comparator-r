\name{alignment_site_summary}
\alias{alignment_site_summary}
\title{Summarize conservation across a DNA multiple alignment}
\usage{
alignment_site_summary(
  aligned,
  min_occupancy = 0.8,
  conserved_threshold = 0.95
)
}
\arguments{
  \item{aligned}{A Biostrings \code{DNAStringSet} with equal aligned widths.}
  \item{min_occupancy}{Minimum fraction of sequences containing a canonical A/C/G/T base for a site to be considered well occupied.}
  \item{conserved_threshold}{Minimum major-base frequency among canonical bases for a well-occupied site to be classified as conserved.}
}
\value{
A data frame with aligned position, occupancy, major base, major-base frequency, Shannon entropy in bits, number of observed canonical bases, and site class.
}
\description{
Computes site-level conservation statistics from a DNA multiple alignment. Gaps and non-canonical symbols reduce occupancy but are excluded from the canonical-base entropy denominator.
}
