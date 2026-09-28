# Download the accession versions frozen in data/genome_manifest.csv.
# Network access is deliberately kept outside package tests.

manifest <- read.csv("data/genome_manifest.csv", stringsAsFactors = FALSE)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

for (accession in manifest$accession) {
  destination <- file.path("data/raw", paste0(accession, ".fasta"))

  query <- paste0(
    "?db=nuccore&id=", utils::URLencode(accession, reserved = TRUE),
    "&rettype=fasta&retmode=text"
  )

  download.file(
    paste0(base_url, query),
    destination,
    mode = "wb",
    quiet = FALSE
  )

  if (!file.exists(destination) || file.info(destination)$size == 0) {
    stop("Download failed for accession: ", accession)
  }

  Sys.sleep(0.4)
}
