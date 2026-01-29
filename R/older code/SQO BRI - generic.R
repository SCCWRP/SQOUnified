# Compute the benthic response index (BRI) score and BRI condition category.
#

#   The BRI is the 4th root abundance weighted pollution tolerance score of the organisms present in a sample relative to the 4th root abundance
#   of all taxa in the sample with a tolerance score (aka p-code, p-value). The higher the BRI score, the more degraded the sample. 
#

#   The BRI is the 4th root relative abundance weighted pollution tolerance score of the organisms present in a benthic sample. The higher
#   the BRI score, the more degraded the benthic community represented by the sample.
#
#   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality Assessment Techincal Support Manual. SCCWRP
#      Techincal Report 777
#    Details on validation of the index can be found in Ranasinghe et al. 2009 Calibration and evaluation of five indicators of benthic community condition
#       in two California bay and estuary habitats. Marine Pollution Bulletin 59:5-13
#     Background concepts of the index can be found in Smith et al  2001. Benthic Resonspe Index for Assessing Infaunal Communities on the Southern
#         California Mainland Shelf. Ecological Applications 11: 1073-1087
#
#
#       NOTE: This code is designed to be run in the Bight Program's Generic SQO BLOE calculator wrapper function. However, the function can
#             be run independently but the user must load the following into their R environment:
#                 1. file_id - this is a quoted string for you to identify the data the BRI scores are associated with
#                   e.g., "Bight 23" or "2024 San Diego Bay" - this will be used to name all of the output files
#                 2. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
#                 3. BenthicData - a data frame containing the benthic data and the station information, detailed above
#
# This function will produce a csv file of final SQO BRI scores for each sample, as well as csv files of interim tables detailing the taxa in the
#     submitted data that have tolerance values assigned and one for taxa in the submitted data that do not have an tolerance value assigned.


SQO.BRI.generic <- function(BenthicData, output_path=output_path, file_id=file_id) 
  
{
  #loading in packages needed to run function
  require(tidyverse)
  #loading in SQO species list that contains p codes, amongst other things

  load("Reference Files/SoCal SQO LU 4_7_20.RData")
  
  #create empty dataframe to populate w/ bri scores
  bri.out.null<-tibble(stationid="dummy",
                       replicate=NaN,
                       sampledate=ymd("2000-01-1"),
                       index="BRI",
                       score=NaN,
                       condition_category=NA,
                       condition_category_score=NA,
                       note=NA)

  #in case a sample had no animals (e.g., taxon=NoOrganismsPresent), we force it into the High Disturbance category.
  #the calculator would not be able to process that sample and would drop it, so we deal with it apriori
  defaunated<-BenthicData %>%
    filter(taxon=="NoOrganismsPresent") %>%
    mutate(index="BRI",score=NaN, condition_category="High Disturbance", condition_category_score=4, note="Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)

  #matching p codes to taxa in the submitted data
  all.for.bri <- BenthicData %>%
    filter(taxon!="NoOrganismsPresent") %>% #removing samples without any animals so the calculator doesn't get confused
  left_join(., select(sqo.list.4_7_20, TaxonName, ToleranceScore), by = c('taxon' = 'TaxonName'))

  #identify the taxa in the submitted data with a tolerance score and how many samples they occur in

  taxa_w_pvalue<-all.for.bri %>%
    group_by(taxon, ToleranceScore) %>%
    summarise(Freq_of_Occ=length(stationid)) %>%
    ungroup() %>%
    drop_na(ToleranceScore)
  #export as an interim file of taxa with a tolerance value that the user can review
  
  write.csv(taxa_w_pvalue, paste(output_path,"/", file_id, " SQO BRI interim 1 - taxa with a tolerance score.csv", sep=""), row.names = FALSE)

  #identify those taxa in the submitted data without a tolerance score and how many samples they occur in
  taxa_wo_pvalue<-all.for.bri%>%
    group_by(taxon, ToleranceScore) %>%
    summarise(Freq_of_Occ=length(stationid)) %>%
    ungroup() %>%
    filter(is.na(ToleranceScore)) %>%
    select(-ToleranceScore)

  #export as an interim file of taxa missing a tolerance value that the user can review
  
  write.csv(taxa_wo_pvalue, paste(output_path, "/", file_id, " SQO BRI interim 2 - taxa without a tolerance score.csv", sep=""), row.names = FALSE)

  
  bri.out<-all.for.bri %>%
  drop_na(ToleranceScore) %>%
  mutate(fourthroot_abun = abundance ** 0.25,
         tolerance_value = fourthroot_abun * ToleranceScore) %>%
  group_by(stationid, sampledate, replicate) %>%
  summarize(numerator = sum(tolerance_value, na.rm = T), denomenator= sum(fourthroot_abun, na.rm = T), score=numerator/denomenator) %>%
    select(stationid, sampledate, replicate, score) %>%
    # Output the BRI category given the BRI score and the thresholds for Southern California Marine Bays
    mutate(
      condition_category = case_when( (score < 39.96) ~ "Reference",
                                (score >= 39.96 & score < 49.15) ~ "Low Disturbance",
                                (score >= 49.15 & score < 73.27) ~ "Moderate Disturbance",
                                (score >= 73.27) ~ "High Disturbance"
    )) %>%
    # Output the BRI category score given the category for thresholds for Southern CA Marine Bays
    mutate(
      condition_category_score = case_when( (condition_category == "Reference") ~ 1,
                                      (condition_category == "Low Disturbance") ~ 2,
                                      (condition_category == "Moderate Disturbance") ~ 3,
                                      (condition_category == "High Disturbance") ~ 4 ),
      #index="BRI", 
      note=NA)

  bri.stations<-BenthicData %>% 
    select(-taxon, -abundance, -exclude) %>% 
    distinct()
  
  
  bri.out.2<-bri.out.null %>%
    bind_rows(bri.out, defaunated) %>% 
    left_join(bri.stations, ., by=c("stationid", "sampledate", "replicate")) %>% 
    mutate(note=case_when(is.na(condition_category)~"Unable to calculate BRI, likely no p-code taxa", 
                          TRUE~note),
           index="BRI")
    

  if(length(bri.out.2$stationid)>1)#if BRI scores are calculated, function will drop the dummy data placeholders and report the data for the submitted samples
  {
    bri.out.3<-bri.out.2 %>%
      filter(stationid!="dummy") 
  }

else
  {
  bri.out.3<-bri.out.2 %>% # if BRI scores are not calculated, the function will only report the dummy data placeholders
    mutate(note="BRI scores not caculated")
}
 
  write.csv(bri.out.3, paste(output_path, "/", file_id, " SQO BRI scores.csv", sep=""), row.names = FALSE) 
  return(bri.out.3)
}


