# Download the frozen v0.3 RefSeq census panel from NCBI.

manifest <- read.csv("inst/extdata/genome_manifest_v0.3.csv", stringsAsFactors = FALSE)
dir.create("data/raw_v03", recursive = TRUE, showWarnings = FALSE)

base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

for (accession in manifest$accession) {
  destination <- file.path("data/raw_v03", paste0(accession, ".fasta"))

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

  header <- readLines(destination, n = 1, warn = FALSE)
  if (length(header) != 1 || !startsWith(header, paste0(">", accession))) {
    stop("Unexpected FASTA header for accession: ", accession)
  }

  Sys.sleep(0.4)
}
