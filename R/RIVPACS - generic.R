#' River Invertebrate Prediction and Classification System (RIVPACS) Index calculation
#' 
#'  The RIVPACS index is an Observed to Expected index comparing the infaunal community observed in a sample to that which would be expected to be in the
#'  location under reference conditions based on a predictive model based on latitude, longitude, and water depth
#'
#' #   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality Assessment Techincal Support Manual. SCCWRP
#      Technical Report 777
#    Details on validation of the index can be found in Ranasinghe et al. 2009 Calibration and evaluation of five indicators of benthic community condition
#       in two California bay and estuary habitats. Marine Pollution Bulletin 59:5-13
#   Background concepts of the index can be found in Van Sickle et al. 2006. Selecting discriminant function models for predicting the expected richness
#       of aquatic macroinvertebrates. Freshwater Biology 51:359-372

# this function, and O:E style indices in general, require an expected taxa model to be run as part of the index calculation. We have not changed the 
# framework of that model since the original conception of the index in 2008. Consequently this function loads and runs that modelling function (SoCalRivpacs v2.R)
# and its associated data. The associated data file (socal.reference.taxa) can be updated as SQO-related taxonomy is updated.

#       NOTE: This code is designed to be run in the Bight Program's Generic SQO BLOE calculator wrapper function. However, the function can
#             be run independently but the user must load the following into their R environment:
#                 1. file_id - this is a quoted string for you to identify the data the RIVPACS scores are associated with
#                   e.g., "Bight 23" or "2024 San Diego Bay" - this will be used to name all of the output files
#                 2. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
#                 3. BenthicData - a data frame containing the benthic data and the station information, detailed above
#
# This function will produce a csv file of final RIVPACS scores for each sample,


RIVPACS.generic <- function(BenthicData, output_path, file_id)
  {
  # R packages needed to run the function
  require(tidyverse)
  # loading in the discriminant function used to assign samples to their reference clusters
  source("Reference Files/SccwrpRivpacs/R/SoCalRivpacs v2.R")

  # RIVPACS model information used to calibrate the discriminant function model, based off of Southern California data
  attach("Reference Files/SccwrpRivpacs/data/SoCalReference.RData")
  attach("Reference Files/SccwrpRivpacs/data/SoCalExample.RData")

  
  #rectifying field and dataframe naming conventions to match pre-existing RIVPACS code
  benthic_data <- BenthicData %>%
    rename(Taxa = taxon) %>%
    mutate(sample_id=paste(stationid, sampledate, replicate, sep="_")) %>%
    filter(Taxa!="NoOrganismsPresent") #removing defaunated samples from analysis

  #selecting sample id information to associate w/ OE scores later

  sampleids<-benthic_data %>%
    select(sample_id, stationid, sampledate, replicate) %>%
    distinct()


  #in case a sample had no animals (e.g., taxon=NoOrganismsPresent), we force it into the High Disturbance category.
  #the calculator would not be able to process that sample and would drop it, so we deal with it apriori
  defaunated<-BenthicData %>%
    filter(taxon=="NoOrganismsPresent") %>%
    mutate(index="Rivpacs",score=NaN, condition_category="High Disturbance", condition_category_score=4, note="Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)

  
  # selecting location and depth information for a sample, which are used to assign it to a reference cluster(s) and thereby predict reference community composition
  scb.predictors<-benthic_data %>%
    select(sample_id, Latitude=latitude, Longitude=longitude, SampleDepth=depth) %>%
    distinct() %>%
    column_to_rownames("sample_id") %>%
    as.matrix()

    # selecting the observed taxa in each sample to compare to the model's expected reference community
   scb.taxa<-benthic_data %>%
    select(sample_id, Taxa, Abundance=abundance) %>%
    mutate(Taxa=str_replace_all(Taxa, "[ ()]", "_")) # changing naming convention to matcht that in the discrimanant function


  
  #Putting the taxa data into a taxon-abundance matrix
   scb.taxa.2 <- scb.taxa %>%
    pivot_wider(id_cols = sample_id, names_from = Taxa,
                       values_from = Abundance,  values_fill = 0) %>%
    column_to_rownames("sample_id")


 
   #submitting abiotic predictors and biotic response data to the RIVPACS discriminant function

   rivpacs.mod<-SoCalRivpacs.2(observed.predictors = scb.predictors, observed.taxa = scb.taxa.2)
 
    # extract the observed taxa and expected taxa information from the RIVPACS function output
   oe.tab<-rivpacs.mod$oe.table %>%
    as_tibble() %>%
    left_join(., sampleids, by="sample_id")


   # gathering the station information associated with each sample to associate theme later
   oe.stations<-BenthicData %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()




  # Assign condition categories based upon O:E scores and thresholds

############

  rivpacs.scores <- oe.tab %>%
    rename(score=O.over.E) %>%
    select(-sample_id) %>%
    mutate(condition_category = case_when((score <= 0.32) ~ "High Disturbance",
                                       ((score > 0.32 & score <= 0.74) | (score >= 1.26)) ~ "Moderate Disturbance",
                                       ((score > 0.74 & score <= 0.90) | score >= 1.10 & score < 1.26) ~ "Low Disturbance",
                                       (score > 0.90 | score < 1.10) ~ "Reference"),
           condition_category_score = case_when(condition_category == "Reference" ~ 1,
                                               condition_category == "Low Disturbance" ~ 2,
                                               condition_category == "Moderate Disturbance" ~ 3,
                                               condition_category == "High Disturbance" ~ 4),
           #David Gillett added this flag to call the user's attention that the sample may be outside the calibration dataset for the RIVPACS model
           #David Gillett thinks this is not necessarily a reason to dismiss the RIVPACS index score, but to be skeptical and examine it further
           note=case_when(outlier.01=='FAIL'~"Caution-Sample Outside RIVPACS habitat model",
                          TRUE~ NA),
           index="Rivpacs") %>%
    bind_rows(defaunated,.) #add the defaunated samples back in

   #output RIVPACS O:E details for the use to review
   write.csv(rivpacs.scores, file=paste(output_path,"/", file_id, " SQO RIVPACS interim file - OE details.csv", sep="" ), row.names = FALSE)

  #clean up scores and add associated station information to match other index outputs
   rivpacs.scores.2<-rivpacs.scores %>%
    select(-O, -E, -outlier.05, -outlier.01) %>%
    left_join(oe.stations, ., by=c("stationid", "sampledate", "replicate"))

   write.csv(rivpacs.scores.2, file=paste(output_path,"/", file_id, " SQO RIVPACS scores.csv", sep="" ), row.names = FALSE)


 

  return(rivpacs.scores.2)
}






