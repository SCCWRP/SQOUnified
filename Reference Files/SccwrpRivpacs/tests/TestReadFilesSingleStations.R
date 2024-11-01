# This test verifies the ability to read in and format .csv files with the
# server script.

library(SccwrpRivpacs)

# Run using input csv files.

ts <- "TestSingleStations"
uf <- "R CMD check Read CSV Test Single Stations"

# Read in user files.

station <- read.csv(paste(system.file(package = "SccwrpRivpacs"), "/extdata/TestSingleStationInfo.csv", sep = ""), stringsAsFactors = FALSE) 

benthic <- read.csv(paste(system.file(package = "SccwrpRivpacs"), "/extdata/TestSingleBenthicInfo.csv", sep = ""), stringsAsFactors = FALSE) 

# Split to SoCal and SFBay.

scb.station <- station[station$HabitatCode == "C", ]

sfb.station <- station[station$HabitatCode == "D", ]

# If data exists for habitat, format data.  

if(nrow(scb.station) > 0) {
  
  scb.predictors <- data.frame(Latitude = scb.station$Latitude,
                               Longitude = scb.station$Longitude,
                               SampleDepth = scb.station$SampleDepth)
  
  row.names(scb.predictors) <- scb.station$StationID
  
  scb.predictors <- as.matrix(scb.predictors)
  
  scb.taxa <- benthic[benthic$StationID %in% scb.station$StationID, ]
  
  scb.taxa$Taxa <- gsub(" ", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub("(", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub(")", "_", scb.taxa$Taxa, fixed = TRUE)
  
  scb.taxa <- reshape(data = scb.taxa, v.names = "Abundance", timevar = "Taxa", 
                      idvar = "StationID", direction = "wide")
  
  row.names(scb.taxa) <- scb.taxa$StationID
  
  scb.taxa <- scb.taxa[, -1]
  
  colnames(scb.taxa) <- gsub("Abundance.", "", colnames(scb.taxa))
  
  # Replace NAs with zero.
  scb.taxa[is.na(scb.taxa)] <- 0
  
  # RIVPACS calculations. By default the functions use the example user data.
  socal <- SoCalRivpacs(observed.predictors = scb.predictors, observed.taxa = scb.taxa)
  
  # Create the HTML files displaying the output.
  HtmlOutput(rivpacs = socal, timestamp = ts, user.filename = uf, path = "")
  
}  

if(nrow(sfb.station) > 0) {
  
  sfb.predictors <- data.frame(SampleDepth = sfb.station$SampleDepth,
                               Hab_G = rep(0, times = nrow(sfb.station)),
                               Longitude = sfb.station$Longitude)
  
  row.names(sfb.predictors) <- sfb.station$StationID
  
  sfb.predictors <- as.matrix(sfb.predictors)
  
  sfb.taxa <- benthic[benthic$StationID %in% sfb.station$StationID, ]
  
  sfb.taxa$Taxa <- gsub(" ", "_", sfb.taxa$Taxa, fixed = TRUE)
  sfb.taxa$Taxa <- gsub("(", "_", sfb.taxa$Taxa, fixed = TRUE)
  sfb.taxa$Taxa <- gsub(")", "_", sfb.taxa$Taxa, fixed = TRUE)
  
  sfb.taxa <- reshape(data = sfb.taxa, v.names = "Abundance", timevar = "Taxa", 
                      idvar = "StationID", direction = "wide")
  
  row.names(sfb.taxa) <- sfb.taxa$StationID
  
  sfb.taxa <- sfb.taxa[, -1]
  
  colnames(sfb.taxa) <- gsub("Abundance.", "", colnames(sfb.taxa))
  
  # Replace NAs with zero.
  sfb.taxa[is.na(sfb.taxa)] <- 0
  
  # RIVPACS calculations. By default the functions use the example user data.
  sfbay <- SFBayRivpacs(observed.predictors = sfb.predictors, observed.taxa = sfb.taxa)
  
  # Create the HTML files displaying the output.
  HtmlOutput(rivpacs = sfbay, timestamp = ts, user.filename = uf, path = "")
  
}

# Checks of O/E output.

# From dput(socal$oe.table$O.over.E)
socal.check <- 1.4286

stopifnot(socal$oe.table$O.over.E == socal.check)

# From dput(sfbay$oe.table$O.over.E)
sfbay.check <- 0.3944

stopifnot(sfbay$oe.table$O.over.E == sfbay.check)
