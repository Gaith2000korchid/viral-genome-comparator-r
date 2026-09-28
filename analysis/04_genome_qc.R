# Build a reproducible QC table from downloaded complete-genome FASTA files.

manifest <- read.csv("inst/extdata/genome_manifest.csv", stringsAsFactors = FALSE)

qc_rows <- lapply(manifest$accession, function(accession) {
  path <- file.path("data/raw", paste0(accession, ".fasta"))
  if (!file.exists(path)) {
    stop("Missing FASTA file: ", path, ". Run analysis/03_download_genomes.R first.")
  }

  qc <- fasta_qc(path)
  if (nrow(qc) != 1) {
    stop("Expected exactly one FASTA record for ", accession)
  }

  data.frame(accession = accession, qc, stringsAsFactors = FALSE)
})

qc <- do.call(rbind, qc_rows)
qc <- merge(manifest, qc, by = "accession", all.x = TRUE, sort = FALSE)

dir.create("results", showWarnings = FALSE)
write.csv(qc, "results/genome_qc.csv", row.names = FALSE)

print(qc[, c("accession", "virus_name", "length", "gc_percent",
             "ambiguous_count", "ambiguous_percent")])
