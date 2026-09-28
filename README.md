# viral-genome-comparator-r

Comparative viral genomics in R: sequence alignment, intergenomic distances, conservation and phylogenetic exploration of bacteriophage genomes.

## Scientific question

**Can whole-genome sequence similarity and alignment recover biologically meaningful relationships among related bacteriophages?**

This repository is being developed as a reproducible R case study in comparative viral genomics. The long-term workflow will connect sequence quality control, alignment, intergenomic distances, site-by-site variation, conservation, clustering and phylogenetic exploration.

## v0.1 scope

The first version deliberately focuses on fundamentals before whole-genome analysis:

- normalize and validate DNA sequences;
- compute length, GC percentage and ambiguous-base statistics;
- compute reverse complements;
- implement global pairwise alignment with Needleman-Wunsch;
- expose the dynamic-programming score matrix for inspection;
- test expected behavior with `testthat`;
- run automated R package checks with GitHub Actions.

The custom Needleman-Wunsch implementation is educational. Its O(n*m) time and memory requirements make it unsuitable as the production aligner for complete viral genomes.

## Roadmap

### v0.1 — Sequence and alignment fundamentals
Short DNA sequences, QC, reverse complement, Needleman-Wunsch, tests and CI.

### v0.2 — Intergenomic comparison
Validated bacteriophage genome panel, production alignment workflow, pairwise similarity/coverage metrics, distance matrix and heatmap.

### v0.3 — Biological interpretation
Site-by-site variation, conservation, clustering, phylogenetic exploration and comparison with curated taxonomy/metadata.

## Reproduce v0.1

```r
install.packages("testthat")
devtools::load_all()
testthat::test_dir("tests/testthat")

sequence_qc("ACGTNN")
reverse_complement("ATGC")
needleman_wunsch("GATTACA", "GCATGCT")
```

## Data policy

No biological genome dataset is frozen in v0.1. Accessions, taxonomy and provenance for the bacteriophage pilot panel will be validated before they become part of the reproducible dataset.

## Inspiration

The project is conceptually inspired by established viral comparative-genomics workflows such as VIRIDIC, ViralMSA, Nextclade and Vclust, while remaining an independent educational implementation. It does not aim to replace those production tools.

## License

MIT.
