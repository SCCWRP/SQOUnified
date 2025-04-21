# Compute the Index of Biotic Integrity (IBI) and IBI condition category.
#

#   The IBI is a a multi-metric index that compares the values of four different metrics to the ranges expected under 
#   reference conditions. The score increases by one for each metric that is outside of the reference range. The four
#   metrics are:
#
#   (1) the total number of taxa - measure of biodiversity 
#   (2) the total number of mollusc taxa - measure of sensitivity to eutrophication and potentially invasive taxa
#   (3) the abundance of Notomastus sp. - measure of the presence of organic matter indicative taxa
#   (4) the number of sensitive taxa. - measure of the presence of pollution sensitive taxa designated by Thompson and Lowe, 2004
#
#   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality Assessment Techincal Support Manual. SCCWRP
#      Technical Report 777
#    Details on validation of the index can be found in Ranasinghe et al. 2009 Calibration and evaluation of five indicators of benthic community condition
#       in two California bay and estuary habitats. Marine Pollution Bulletin 59:5-13
#   Background concepts of the index can be found in Thompson and Lowe 2004. Assessment of macrobenthos response to sediment contamination in the
#       San Francisco Estuary, California USA. Environmental Toxicology and Chemistry 23:2178-2187
#   
#       NOTE: This code is designed to be run in the Bight Program's Generic SQO BLOE calculator wrapper function. However, the function can
#             be run independently but the user must load the following into their R environment:
#                 1. file_id - this is a quoted string for you to identify the data the IBI scores are associated with
#                   e.g., "Bight 23" or "2024 San Diego Bay" - this will be used to name all of the output files
#                 2. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
#                 3. BenthicData - a data frame containing the benthic data and the station information, detailed above
#
# This function will produce a csv file of final IBI scores for each sample, as well as csv files of interim tables detailing the taxa in the
#     submitted data and how they are classified by the index, as well as, the benthic data classified into IBI metrics.

IBI.generic <- function(BenthicData, output_path, file_id)
  
{
  require(tidyverse)
  require(naniar)
  load("Reference Files/SoCal SQO LU 4_7_20.RData")

  #create an empty dataframe to populate with IBI scores
  ibi.out.null<-tibble(stationid="dummy",
                       replicate=NaN,
                       sampledate=ymd("2000/01/1"),
                       index="IBI",
                       score=NaN,
                       condition_category=NA,
                       condition_category_score=NA,
                       note=NA)
  #in case a sample had no animals (e.g., taxon=NoOrganismsPresent), we force it into the High Disturbance category.
  #the calculator would not be able to process that sample and would drop it, so we deal with it apriori
  defaunated<-BenthicData %>%
    filter(taxon=="NoOrganismsPresent") %>%
    mutate(index="IBI",score=NaN, condition_category="High Disturbance", condition_category_score=4, note="Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)

  # Prepare the given data frame so that we can compute the IBI score and categories
  ibi_data <- BenthicData %>%
    filter(exclude!="Yes") %>%
    left_join(sqo.list.4_7_20, by = c("taxon"="TaxonName")) %>%
    replace_na_with(., "") %>%
    filter(taxon!="NoOrganismsPresent")

  #Export data so the user knows what is going to used in subsequent calculations
  
  ibi_data.review<-ibi_data %>% 
    mutate(Notomastus_flag=if_else(str_detect(taxon,"Notomastus"), 1, 0)) %>% 
    select(stationid, sampledate, replicate, taxon, abundance, exclude, Mollusc, IBISensitive, Notomastus_flag)
    
  write.csv(ibi_data.review, file = paste(output_path, "/", file_id, " SQO IBI interim 1 - data to be analyzed with SQO designations assigned.csv", sep=""),
            row.names = FALSE)



  # Calculate taxa richness
  ibi1 <- ibi_data %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfTaxa =length(taxon))

 
  # Calculate mollusc taxa richness
  ibi2 <- ibi_data %>%
    mutate(flag=(if_else(Mollusc=="Mollusc", 1, 0))) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfMolluscTaxa=sum(flag)) %>%
    ungroup()


  # calculate Notomastus spp. abundance
  ibi3 <- ibi_data %>%
    mutate(flag=if_else(str_detect(taxon,"Notomastus"), abundance, 0)) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NotomastusAbun = sum(flag)) %>%
    ungroup()


  # Calculate % Sensitive Taxa
  ibi4 <- ibi_data %>%
    mutate(flag=if_else(IBISensitive=="S", 1, 0) )%>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(sensitive_S=sum(flag)) %>%
    ungroup() %>%
    left_join(ibi1, by=c("stationid", "sampledate", "replicate")) %>%
    mutate(PctSensTaxa=(sensitive_S/NumOfTaxa)*100) %>%
    select(stationid, sampledate, replicate, PctSensTaxa)

 
  ### Reference ranges for IBI metrics in Southern California Marine Bays against which sample values are compared
  ### [ Table 4.19 CASQO Technical Manual 3rd edition 2021 - page 68 ]
  ibi_ref_ranges_table <- data.frame(ref_low = c(13, 2, 0, 19),
                                     ref_high = c(99, 25, 59, 47.1))
  row.names(ibi_ref_ranges_table) <- c("NumOfTaxa", "NumOfMolluscTaxa", "NotomastusAbun", "PctSensTaxa")


  ### IBI Metrics:
  # We stitch together all the necessary IBI metrics to determine the IBI index.
  # Each of the metrics is then compared to the tables listed above (Table 5.4 and Table 5.5) to determine the IBI score,
  # the IBI Category, and IBI Category Score
  ibi_metrics <- ibi1 %>%
    full_join(ibi2, by = c( "sampledate", "stationid", "replicate")) %>%
    full_join(ibi3, by = c( "sampledate", "stationid", "replicate")) %>%
    full_join(ibi4, by = c( "sampledate", "stationid", "replicate"))

  # Export a table of the data classified into IBI metrics for the user to reveiw before they are integrated
  write.csv(ibi_metrics, file = paste(output_path, "/", file_id, " SQO IBI interim 2 - IBI metric values.csv", sep=""),
            row.names = FALSE)
  
  
  
  ### Reference ranges for IBI metrics in Southern California Marine Bays against which sample values are compared
  ### [ Table 4.19 CASQO Technical Manual 3rd edition 2021 - page 68 ]
  ibi_ref_ranges_table <- data.frame(metric = c("NumOfTaxa", "NumOfMolluscTaxa", "NotomastusAbun", "PctSensTaxa"),
                                     ref_low = c(13, 2, 0, 19),
                                     ref_high = c(99, 25, 59, 47.1))
 
  
  # Calculate IBI scores
  ibi.scores<-ibi_metrics %>% 
    #rearrange data to make scoring cleaner
    pivot_longer(., cols=c(-stationid, -sampledate, -replicate), names_to = "metric", values_to = "value") %>% 
    #add in the reference site ranges for comparison
    left_join(.,ibi_ref_ranges_table, by="metric") %>% 
    #IBI score starts at 0 and for each metric out of reference rance (< or >), the score increases by 1
    mutate(out_of_range=case_when(value<ref_low | value>ref_high ~1,
                                  TRUE~0)) %>% 
    group_by(stationid, sampledate, replicate) %>% 
    summarise(score=sum(out_of_range)) %>% 
    ungroup() %>% 
    #convert IBI scores to SQO categories and category scores
    mutate(index="IBI", .before=score) %>% 
    mutate(condition_category = case_when(score == 0 ~ "Reference", 
                                          score == 1 ~ "Low Disturbance", 
                                          score == 2 ~ "Moderate Disturbance", 
                                          score %in%c(3,4) ~ "High Disturbance"),
           condition_category_score = case_when(score == 0 ~ 1, 
                                                score == 1 ~ 2, 
                                                score == 2 ~ 3, 
                                                score %in%c(3,4) ~ 4))
  
 
  #gathering the station information (i.e., non-taxonomic data) for each site

  ibi.stations<-BenthicData %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()

  ibi.out<-ibi.scores %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score) %>%
    bind_rows(ibi.out.null, ., defaunated) %>%
    full_join(ibi.stations, ., by=c("stationid", "sampledate", "replicate"))

  if(length(ibi.out$stationid)>1)#if IBI scores are calculated, function will drop the dummy data placeholders and report the data for the submitted samples
  {
    ibi.out.2<-ibi.out %>%
      filter(stationid!="dummy") %>%
      mutate(note=if_else(is.na(note), "none", note))
  }

  else
  {
    ibi.out.2<-ibi.out %>% # if IBI scores are not calculated, the function will only report the dummy data placeholders
      mutate(note="IBI scores not caculated")
  }




 
# Export the final scores to the location designated by the output path   
  write.csv(ibi.out.2, paste(output_path, "/", file_id, " SQO IBI scores.csv", sep=""), row.names = FALSE)

# Export the final scores to the general R environment
   return(ibi.out.2)


}

