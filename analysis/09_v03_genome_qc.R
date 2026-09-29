# QC the expanded v0.3 RefSeq panel.

manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)

qc_rows <- lapply(seq_len(nrow(manifest)), function(i) {
  accession <- manifest$accession[i]
  path <- file.path("data/raw_v03", paste0(accession, ".fasta"))

  if (!file.exists(path)) {
    stop("Missing FASTA file: ", path, ". Run analysis/08_download_v03_panel.R first.")
  }

  qc <- fasta_qc(path)
  if (nrow(qc) != 1) {
    stop("Expected exactly one FASTA record for ", accession)
  }

  observed_length <- qc$length[1]
  if (observed_length != manifest$expected_length_bp[i]) {
    stop(
      "Length mismatch for ", accession,
      ": expected ", manifest$expected_length_bp[i],
      " bp but downloaded ", observed_length, " bp."
    )
  }

  data.frame(accession = accession, qc, stringsAsFactors = FALSE)
})

qc <- do.call(rbind, qc_rows)
qc <- merge(manifest, qc, by = "accession", all.x = TRUE, sort = FALSE)

dir.create("results", showWarnings = FALSE)
write.csv(qc, "results/v0.3_genome_qc.csv", row.names = FALSE)

print(qc[, c(
  "accession", "virus_name", "current_genus", "length",
  "gc_percent", "ambiguous_count"
)])
