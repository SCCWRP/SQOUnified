# This test verifies the taxa names in the reference files exactly match the
# taxa names in the template.

library(SccwrpRivpacs)

# Read in taxa reference file.

taxa.reference <- read.csv(paste(system.file(package = "SccwrpRivpacs"), 
                                 "/extdata/TaxaList.csv", sep = ""), 
                           stringsAsFactors = FALSE)$TaxonName

taxa.reference <- gsub(" ", "_", taxa.reference, fixed = TRUE)
taxa.reference <- gsub("(", "_", taxa.reference, fixed = TRUE)
taxa.reference <- gsub(")", "_", taxa.reference, fixed = TRUE)

# SoCal test.

socal <- colnames(socal.reference.taxa)

socal.matches <- socal %in% taxa.reference

# Non-matching taxa.
socal[!(socal %in% taxa.reference)]

# Do all taxa in the habitat reference match entries in the taxa reference?
stopifnot(length(socal) == length(socal.matches[socal.matches == TRUE]))

sfbay <- colnames(sfbay.reference.taxa)

sfbay.matches <- sfbay %in% taxa.reference

# Non-matching taxa.
sfbay[!(sfbay %in% taxa.reference)]

# Do all taxa in the habitat reference match entries in the taxa reference?
stopifnot(length(sfbay) == length(sfbay.matches[sfbay.matches == TRUE]))
