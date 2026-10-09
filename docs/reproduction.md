# Frozen v0.3 analysis environment

Use Linux or Ubuntu under WSL2. The direct dependency specification is `environment.yml`; the exact Linux package/build resolution is `environment-linux-64.explicit.txt`.

```bash
micromamba create -y -n viral-genome-comparator -f environment-linux-64.explicit.txt
micromamba run -n viral-genome-comparator R CMD INSTALL .
micromamba run -n viral-genome-comparator Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'
micromamba run -n viral-genome-comparator Rscript analysis/run_v03.R
```

R 4.4.3, Biostrings 2.74.0, DECIPHER 3.2.0 and MAFFT 7.525 are specified. The explicit file also fixes their resolved dependencies. It does not provide a Windows-native/macOS lock or a bit-identical OS environment. Linux package downloads, including the Bioconductor annotation-data post-install, still need network access.

The runner loads `viralGenomeComparator` before sourcing analysis functions. Run it from the repository root, not by passing its path while working in an unrelated directory. It downloads exact accession versions, aligns the panel, computes metrics, tests label concordance, summarizes conservation and generates an exploratory tree. It records `sessionInfo()` and software versions alongside the generated results.

GitHub Actions uses the same explicit environment for package checks and the v0.3 analysis. Full FASTA/alignment/distance outputs are workflow artifacts; the small committed result tables are the dated interpretation baseline. New executions can differ if an algorithm version differs from the historical run; inspect recorded versions before replacing a numerical baseline.

## Metric definitions

`identity_percent` compares canonical A/C/G/T bases at columns where both sequences provide one. Coverage measures shared non-gap columns relative to each sequence's ungapped length. DECIPHER's `method="overlap"`, `includeTerminalGaps=FALSE`, `penalizeGapLetterMatches=TRUE`, `correction=NA` defines the separately reported distance. It is not automatically `1 - identity`, ANI or VIRIDIC similarity.

Nearest-neighbor label agreement is descriptive within the selected 29-genome panel. The 100% value is not external test-set classification accuracy. The permutation test permutes genome labels; it does not treat all pairwise distances as independent observations. The Neighbor-Joining plot is exploratory and is not a supported maximum-likelihood phylogeny.
