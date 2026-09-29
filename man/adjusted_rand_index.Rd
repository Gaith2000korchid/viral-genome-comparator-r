\name{adjusted_rand_index}
\alias{adjusted_rand_index}
\title{Compute the Adjusted Rand Index}
\usage{
adjusted_rand_index(truth, predicted)
}
\arguments{
  \item{truth}{Vector of reference categorical labels.}
  \item{predicted}{Vector of predicted or sequence-derived cluster labels.}
}
\value{A numeric Adjusted Rand Index. Identical non-degenerate partitions return 1.}
\description{
Compares two partitions while adjusting the Rand Index for agreement expected by chance. The v0.3 workflow uses this function to compare whole-genome distance clustering with curated genus labels.
}
