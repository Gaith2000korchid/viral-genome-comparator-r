\name{adjusted_rand_index}
\alias{adjusted_rand_index}
\title{Compute the Adjusted Rand Index}
\usage{
adjusted_rand_index(truth, predicted)
}
\arguments{
  \item{truth}{Reference group labels.}
  \item{predicted}{Cluster labels to compare with the reference grouping.}
}
\value{A numeric Adjusted Rand Index. A value of 1 indicates identical partitions.}
\description{
Compares two partitions while correcting the Rand Index for agreement expected by chance. Cluster label names themselves do not need to match.
}
