# viral-genome-comparator-r

Comparative viral genomics in R: sequence alignment, intergenomic distances, conservation and phylogenetic exploration of bacteriophage genomes.

## Scientific question

**Can whole-genome sequence similarity and alignment recover biologically meaningful relationships among related bacteriophages?**

This repository is a reproducible R case study in comparative viral genomics. The workflow connects sequence quality control, multiple sequence alignment, pairwise identity and coverage, whole-genome distance, clustering and later phylogenetic interpretation.

## Current status

### v0.1 — sequence and alignment fundamentals

- DNA normalization and validation;
- GC and ambiguous-base QC;
- reverse complement;
- educational Needleman-Wunsch implementation;
- unit tests and GitHub Actions R CMD check.

The custom Needleman-Wunsch implementation is intentionally limited to short teaching examples because its O(n*m) time and memory costs make it unsuitable for complete viral genomes.

### v0.2 — pilot intergenomic comparison

The current pilot panel contains six complete bacteriophage genomes from two current genera:

- **Teseptimavirus**: T7 and phiA1122;
- **Teetrevirus**: T3, phiYeO3-12, phiSG-JL2 and vB_YenP_AP5.

The accession versions and taxonomy metadata are frozen in `inst/extdata/genome_manifest.csv`.

The v0.2 workflow:

1. downloads the frozen accession versions from NCBI;
2. computes genome-level QC;
3. aligns the six complete genomes with DECIPHER;
4. computes pairwise nucleotide identity and alignment coverage separately;
5. computes an uncorrected p-distance matrix;
6. performs average-linkage hierarchical clustering;
7. renders a distance heatmap and dendrogram.

Identity and coverage are intentionally kept separate so that high identity over a short shared region is not misinterpreted as high whole-genome similarity.

## Reproduce the analyses

Install the required packages:

```r
install.packages(c("testthat", "devtools"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Biostrings", "DECIPHER"))

devtools::load_all()
```

Run the v0.2 workflow in order:

```r
source("analysis/03_download_genomes.R")
source("analysis/04_genome_qc.R")
source("analysis/05_align_genomes.R")
source("analysis/06_pairwise_metrics.R")
source("analysis/07_cluster_heatmap.R")
```

Expected outputs include:

- `results/genome_qc.csv`
- `results/pilot_alignment.fasta`
- `results/pairwise_alignment_metrics.csv`
- `results/p_distance_matrix.csv`
- `results/cluster_order.csv`
- `results/p_distance_heatmap.png`
- `results/p_distance_dendrogram.png`

## Data policy

The exact accession versions are frozen in the manifest. Raw FASTA files are downloaded programmatically from NCBI and ignored by Git. Derived results can be versioned so that the biological analysis remains auditable without silently replacing source genomes.

## Roadmap

### v0.3 — biological interpretation

The next version will add site-by-site variation, conservation, comparison with curated taxonomy and a more explicit phylogenetic interpretation.

## Inspiration

The project is conceptually inspired by established viral comparative-genomics workflows such as VIRIDIC, ViralMSA, Nextclade and Vclust, while remaining an independent educational implementation. It does not aim to replace those production tools.

## License

MIT.
