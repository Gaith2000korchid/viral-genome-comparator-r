# Changelog

## 0.3.0 — biological interpretation

- expanded the analysis from a six-genome pilot to a 29-genome RefSeq panel;
- added a dated v0.3 accession/provenance manifest;
- moved production whole-genome alignment to MAFFT;
- added pairwise identity, coverage and p-distance outputs;
- added within-genus versus between-genus distance analysis;
- added nearest-neighbor genus concordance;
- added Adjusted Rand Index and a 999-permutation label test;
- added site-level conservation, occupancy and Shannon-entropy summaries;
- added windowed conservation profiles;
- added exploratory Neighbor-Joining tree output;
- added explicit diagnostics for the divergent three-genome Teetrevirus branch;
- added a dedicated end-to-end v0.3 GitHub Actions workflow;
- kept R CMD check and unit tests green.

## 0.2.0 — intergenomic pilot

- added a frozen six-genome bacteriophage panel;
- added automated NCBI FASTA retrieval;
- added Biostrings FASTA QC;
- added complete-genome multiple alignment;
- added pairwise identity and coverage metrics;
- added uncorrected p-distance matrix, clustering, heatmap and dendrogram.

## 0.1.0 — sequence fundamentals

- added DNA normalization and validation;
- added reverse complement and sequence QC;
- added transparent Needleman-Wunsch global alignment;
- added testthat coverage and R CMD check CI.
