ConvertSoCalFiles <- function() {
  
  # ----- Reference Files -----

  # Observed presence/absence (1/0) at reference sites for all taxa. Remove
  # Group and Site colums, and convert to presence/absence (0/1).
  socal.reference.taxa <- read.table("SoCalReferenceTaxa.txt", header = TRUE, stringsAsFactors = FALSE) 
  socal.reference.taxa <- socal.reference.taxa[, !names(socal.reference.taxa) %in% c("Group", "Site")]
  socal.reference.taxa[socal.reference.taxa > 1] <- 1 
  #save(socal.reference.taxa, file = "SoCalReferenceTaxa.Rdata")
  
  # Group IDs of reference sites in socal.reference.taxa. Assign names from socal.reference.taxa.
  socal.reference.groups <- as.vector(read.table("SoCalReferenceGroups.txt", stringsAsFactors = FALSE)[1, ], mode = "integer")
  names(socal.reference.groups) <- row.names(socal.reference.taxa)
  #save(socal.reference.groups, file = "SoCalReferenceGroups.Rdata")
  
  # Group (cluster) means at reference sites.
  socal.reference.group.means <- as.matrix(read.table("SoCalReferenceGroupMeans.txt", header = TRUE, stringsAsFactors = FALSE))
  #save(group.means, file = "SoCalReferenceGroupMeans.Rdata")
  
  # Inverse of pooled covariance matrix among final predictor variables at
  # reference sites.
  socal.reference.covariance <- as.matrix(read.table("SoCalReferenceCovarianceMatrix.txt", header = TRUE, stringsAsFactors = FALSE))
  #save(reference.cov, file = "SoCalReferenceCovarianceMatrix.Rdata")
  
  
  save(socal.reference.taxa, socal.reference.groups, socal.reference.group.means, socal.reference.covariance, file = "SoCalReference.RData")
  
  
  
  
  # ----- Example User Files -----
  
  # Observed species matrix of sites (rows) by taxa (columns).
  socal.example.taxa <- read.table("SoCalExampleTaxa.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  #save(observed.taxa, file = "SoCalExampleTaxa.Rdata")
  
  # Predictor variables (columns) at all new sites/samples (rows).
  socal.example.habitat <- as.matrix(read.table("SoCalExampleHabitat.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1))
  #save(observed.predictors, file = "SoCalExampleHabitat.Rdata")
  
  
  save(socal.example.taxa, socal.example.habitat, file = "SoCalExample.RData")
    
}
