load("R/sysdata.rda")
suppressMessages(library(dplyr))
no.pcode <- data.frame(taxon = c("Acila castrenis","Spiophanes bombx"), stringsAsFactors=FALSE)
cat("fuzzyjoin available:", requireNamespace("fuzzyjoin", quietly = TRUE), "\n")
if (requireNamespace("fuzzyjoin", quietly = TRUE)) {
  res <- fuzzyjoin::stringdist_left_join(no.pcode, ed.14.ptaxa, by = "taxon", max_dist = 2)
  res <- dplyr::select(res, submitted_taxon = taxon.x, did_you_mean_taxon = taxon.y, p_code)
  print(res)
  cat("ncol:", ncol(res), "names:", paste(names(res),collapse=","), "\n")
}
