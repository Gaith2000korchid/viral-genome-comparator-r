# Pilot genome panel

This pilot panel is intentionally small. It is designed to test the v0.2 workflow before scaling.

## Question

Do genome-wide sequence distances separate historically "T7-like" bacteriophages in a way that is consistent with their current genus-level taxonomy?

## Design

The panel contains two current genera in the subfamily Studiervirinae:

- **Teseptimavirus**: T7 and phiA1122.
- **Teetrevirus**: T3, phiYeO3-12, phiSG-JL2 and vB_YenP_AP5.

The historical label "T7-like" is not used as a current taxonomic assignment. Current taxonomy and sequence identifiers are kept as separate metadata fields because database organism names and ICTV taxon names can differ.

## Reproducibility rule

The accession version (for example, `NC_001604.1`) is frozen in `data/genome_manifest.csv`. Downstream analyses must report the accession versions actually used.

Raw FASTA files should be obtained programmatically from NCBI and must not be silently replaced by newer accession versions.

## Scope

Six genomes are enough for a pilot but not for a taxonomic benchmark. The panel can be expanded only after the alignment and distance definitions are validated.
