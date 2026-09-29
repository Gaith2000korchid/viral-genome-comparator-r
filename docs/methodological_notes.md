# Methodological notes — v0.3

## Scope

v0.3 is a focused comparative-genomics case study using 29 complete RefSeq bacteriophage genomes assigned in the frozen manifest to two genera: Teseptimavirus and Teetrevirus.

The analysis asks whether whole-genome sequence structure is concordant with these metadata labels. It does **not** attempt to redefine viral taxonomy or establish classification thresholds.

## Frozen inputs

The file `inst/extdata/genome_manifest_v0.3.csv` records exact accession versions and the snapshot date. Raw FASTA files are downloaded programmatically and are not committed to Git.

This prevents a later database update from silently changing the biological input while preserving a reproducible route back to the source records.

## Alignment

v0.3 uses MAFFT for the complete-genome multiple alignment. The repository retains the custom Needleman-Wunsch implementation only as an educational demonstration on short sequences.

The two should not be conflated: the custom implementation is intentionally transparent but quadratic in time and memory, whereas the v0.3 analysis requires a production-scale aligner.

## Identity, coverage and distance

The workflow reports pairwise nucleotide identity and alignment coverage separately.

The DECIPHER distance used downstream is not interpreted as `1 - identity`. Gap-containing alignment columns can affect the distance calculation, so identity, coverage and distance are kept as distinct outputs.

The distance should also not be called ANI or VIRIDIC similarity. Those terms refer to specific methods and conventions that are not implemented here.

## Taxonomy concordance

Three complementary summaries are used:

1. **nearest-neighbor genus accuracy** — whether the nearest genome in the distance matrix has the same genus label;
2. **Adjusted Rand Index (ARI)** after cutting the hierarchical tree into the same number of groups as the manifest;
3. **permutation test** comparing observed within-genus versus between-genus distance separation with randomized labels.

These measures answer different questions. In v0.3, nearest-neighbor accuracy is 100% and the permutation test is significant, while a forced two-cluster ARI is slightly negative. That combination is possible because three divergent Teetrevirus form a separate global branch.

The ARI is therefore reported, not discarded.

## Conservation

Site-level conservation is summarized from the multiple alignment using canonical A/C/G/T observations, occupancy thresholds, major-base frequency and Shannon entropy.

Low-occupancy alignment columns are kept distinct from genuinely variable high-occupancy sites. This prevents gap-rich regions from being interpreted as ordinary nucleotide diversity.

## Tree interpretation

The Neighbor-Joining tree is an exploratory visualization of the distance matrix.

It is **not** a maximum-likelihood or Bayesian phylogeny and should not be used to make strong evolutionary claims about branching order, ancestral states or divergence times.

## Generalization limits

The current panel contains only two genera from a related bacteriophage context. Results should not be generalized to all bacteriophages or all viral families.

A future expansion would require broader taxon sampling, explicit outgroups or rooting strategy, sensitivity analyses across alignment/distance methods, and comparison against dedicated viral-genomics tools.
