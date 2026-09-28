\name{fasta_qc}
\alias{fasta_qc}
\title{Compute QC statistics for DNA records in a FASTA file}
\usage{
fasta_qc(fasta_path)
}
\arguments{
  \item{fasta_path}{Path to a FASTA file containing DNA sequences.}
}
\value{
A data frame with one row per FASTA record and columns for record name, sequence length, GC percentage, ambiguous-base count, and ambiguous-base percentage.
}
\description{
Reads DNA sequences with Biostrings and computes basic genome-level quality-control statistics. GC percentage uses canonical A, C, G, and T bases as its denominator; all other letters are counted as ambiguous.
}
