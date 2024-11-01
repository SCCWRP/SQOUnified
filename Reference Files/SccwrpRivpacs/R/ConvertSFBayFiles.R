ConvertSFBayFiles <- function() {
  
  # ----- Reference Files -----

  # Observed presence/absence (1/0) at reference sites for all taxa. Remove
  # Group and Site colums, and convert to presence/absence (0/1).
  sfbay.reference.taxa <- read.table("SFBayReferenceTaxa.txt", header = TRUE, stringsAsFactors = FALSE) 
  sfbay.reference.taxa <- sfbay.reference.taxa[, !names(sfbay.reference.taxa) %in% c("Group", "Site")]
  sfbay.reference.taxa[sfbay.reference.taxa > 1] <- 1 
  
  # Group IDs of reference sites in sfbay.reference.taxa. Assign names from sfbay.reference.taxa.
  sfbay.reference.groups <- as.vector(read.table("SFBayReferenceGroups.txt", stringsAsFactors = FALSE)[1, ], mode = "integer")
  names(sfbay.reference.groups) <- row.names(sfbay.reference.taxa)
  
  # Group (cluster) means at reference sites.
  sfbay.reference.group.means <- as.matrix(read.table("SFBayReferenceGroupMeans.txt", header = TRUE, stringsAsFactors = FALSE))
  
  # Inverse of pooled covariance matrix among final predictor variables at
  # reference sites.
  sfbay.reference.covariance <- as.matrix(read.table("SFBayReferenceCovarianceMatrix.txt", header = TRUE, stringsAsFactors = FALSE))

  
  save(sfbay.reference.taxa, sfbay.reference.groups, sfbay.reference.group.means, sfbay.reference.covariance, file = "SFBayReference.RData")
  
 
  
  # ----- Example User Files -----
  
  # Observed species matrix of sites (rows) by taxa (columns).
  sfbay.example.taxa<- read.table("SFBayExampleTaxa.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
  # Predictor variables (columns) at all new sites/samples (rows).
  sfbay.example.habitat <- as.matrix(read.table("SFBayExampleHabitat.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1))
  
  
  save(sfbay.example.taxa, sfbay.example.habitat, file = "SFBayExample.RData")
    
}
