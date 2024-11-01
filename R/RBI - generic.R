# Compute the relative benthic index (RBI) score.
#


#   The RBI is a multi-metric index calculated from the weighted sum of: (a) four community metrics related to biodiversity (total number of taxa, number of crustacean taxa, abundance
#   of crustacean individuals, and number of mollusc taxa), (b) abundances of three positive indicator taxa, and (c) the presence of two negative
#   indicator species.
#
#   The metrics used in the  RBI are:
#   (1) Total number of taxa - a measure of biodiversity,
#   (2) Number of mollusc taxa - measure of taxa sensitive to eutrophication and potentially invasive taxa,
#   (3) Number of crustacean individuals - measure of taxa sensitive to pesticides and low DO, as well potentially invasive taxa,
#   (4) Number of individuals of \emph{Monocorophium insidiosum} - abundance of a pollution sensitive amphipod
#   (5) Number of individuals of \emph{Asthenothaerus diegensis} - abundance of a pollution sensitive bivalve,
#   (6) Number of individuals of \emph{Goniada littorea} - abundance of a pollution sensitive polychaete
#   (7) Whether the data has the presence of \emph{Capitella capitata} complex - presence of pollution indicative polychaete
#   (8) Whether the data has the presence of Oligochaeta - presence of pollution indicative annelids
#
#   In brief, the observed values of each metric are scaled relative to values observed at reference sites in the southern California calibration data set.
#   The scaled metric values are then combined into three meta-metrics: 1-3 as Taxa Richness Weighted Value (TWV), 4-6 as Positive Indicator Taxa (PIT), and 7-8
#   as Negative Indicator Taxa (NIT). These three meta-metrics are combined and scaled to values observed at reference sites in the southern California 
#   calibration data set. 
#   
#   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality Assessment Techincal Support Manual. SCCWRP
#      Technical Report 777
#    Details on validation of the index can be found in Ranasinghe et al. 2009 Calibration and evaluation of five indicators of benthic community condition
#       in two California bay and estuary habitats. Marine Pollution Bulletin 59:5-13
#     Background concepts of the index can be found in Hunt et al  2001. A large-scale categorization of sites in San Francisco Bay, USA,
#         based on the sediment quality triad, toxicity identification evaluations, and gradient studies. Environmental Toxicology and Chemistry 20:1252-1265



RBI.generic <- function(BenthicData, output_path, file_id)
{

  
require(tidyverse)
require(naniar)
  load("Reference Files/SoCal SQO LU 4_7_20.RData")

  #create an empty dataframe to populate with RBI scores
  rbi.out.null<-tibble(stationid="dummy",
                       replicate=NaN,
                       sampledate=ymd("2000/01/1"),
                       index="RBI",
                       score=NaN,
                       condition_category=NA,
                       condition_category_score=NA,
                       note=NA)
  #incase a sample had no animals (e.g., taxon=NoOrganismsPresent), we force it into the High Disturbance category.
  #the calculator would not be able to process that sample and would drop it, so we deal with it apriori
  defaunated<-BenthicData %>%
    filter(taxon=="NoOrganismsPresent") %>%
    mutate(index="RBI",score=NaN, condition_category="High Disturbance", condition_category_score=4, note="Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)


  # Prepare the given data frame so that we can compute the RBI score and categories
  rbi_data <- BenthicData %>%
    filter(exclude!="Yes") %>%
    left_join(sqo.list.4_7_20, by = c("taxon"="TaxonName")) %>%
    replace_na_with(., "") %>%
    filter(taxon!="NoOrganismsPresent")

 #Export data so the user knows what is going to used in subsequent calculations
   write.csv(rbi_data, file = paste(output_path, "/", file_id, " SQO RBI interim 1 - data to be analyzed with SQO designations assigned.csv", sep=""),
            row.names = FALSE)

  
  ####Calculate the different metrics for the RBI
  # calculate taxa richness
  rbi1<- rbi_data %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfTaxa = length(taxon)) %>%
    ungroup()



  # calculate mollusc taxa richness
  rbi2 <- rbi_data %>%
    mutate(flag=(if_else(Mollusc=="Mollusc", 1, 0))) %>%
    group_by(stationid, sampledate, replicate) %>%
      summarise(NumOfMolluscTaxa=sum(flag)) %>%
    ungroup()

  
  # calculate crustacean richness
  rbi3 <- rbi_data %>%
    mutate(flag=if_else(Crustacean=="Crustacean", 1,0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(NumOfCrustaceanTaxa = sum(flag)) %>%
      ungroup()

 
  # calculate crustacean abundance
  rbi4 <- rbi_data %>%
    mutate(flag=if_else(Crustacean=="Crustacean", abundance,0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(CrustaceanAbun = sum(flag)) %>%
    ungroup()

  
  # calculate abundance of M. insidiosum
  rbi5 <- rbi_data %>%
    mutate(flag=if_else(taxon == "Monocorophium insidiosum",abundance, 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(M_insidiosumAbun = sum(flag)) %>%
    ungroup()

  
  # calculate abundance of A. diegensis
  rbi6 <- rbi_data %>%
    mutate(flag=if_else(taxon == "Asthenothaerus diegensis",abundance, 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(A_diegensisAbun = sum(flag)) %>%
    ungroup()

  
  # calculate abundance of G. littorea
  rbi7 <- rbi_data %>%
    mutate(flag=if_else(taxon == "Goniada littorea",abundance, 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(G_littoreaAbun = sum(flag)) %>%
    ungroup()

  
  # calculate negative indicator taxa score (NIT)
  rbi8 <- rbi_data %>%
   mutate(badness=if_else(taxon %in% c( "Capitella capitata Cmplx","Oligochaeta"), -0.1, 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(NIT = sum(badness)) %>%
    ungroup()

 # combine all the metrics into one data frame
  rbi_metrics <- rbi1 %>%
    dplyr::full_join(rbi2, by = c( "stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi3, by = c( "stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi4, by = c( "stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi5, by = c( "stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi6, by = c( "stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi7, by = c( "stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi8, by = c( "stationid", "replicate", "sampledate")) 
   
  #Export an interim file with all rbi metrics for each sample
  write.csv(rbi_metrics, file=paste(output_path, "/", file_id, " SQO RBI interim 2 - raw RBI metrics.csv", sep=""), row.names = FALSE)
  
  ### RBI Category Thresholds for Southern California Marine Bays
  RBI_category_thresholds <- data.frame(ref_low = c(0.27, 0.16, 0.08, 0.08),
                                        ref_high = c(0.27, 0.27, 0.16, 0.08),
                                        condition_category = as.factor(c("Reference",
                                                               "Low Disturbance",
                                                               "Moderate Disturbance",
                                                               "High Disturbance")),
                                        condition_category_score = c(1, 2, 3, 4))

  
  # Scale observed values to the maxima observed in index calibration dataset
  
  rbi_scaled <- rbi_metrics %>%
    #scale the richness metrics
    mutate(scaled_NumTaxa = (NumOfTaxa/99),
           scaled_NumMolluscTaxa = (NumOfMolluscTaxa/28),
           scaled_NumCrustaceanTaxa = (NumOfCrustaceanTaxa/29),
           scaled_CrustaceanAbun = (CrustaceanAbun/1693),
    # calculate Taxa Richness Weighted Value (TWV) from the four richness/abundance type metrics
           TWV = (scaled_NumTaxa + scaled_NumMolluscTaxa + scaled_NumCrustaceanTaxa + (0.25 * scaled_CrustaceanAbun)),
    # scale the indicator taxa
            scaled_M_insid=(((M_insidiosumAbun)^(1/4)) / ((473)^(1/4))),
            scaled_A_dieg=(((A_diegensisAbun)^(1/4) )/ ((27)^(1/4) )),
            scaled_G_littor=(((G_littoreaAbun)^(1/4) )/ ((15)^(1/4) )),
    # calculate Positive Indicator Taxa (PIT)
            PIT=scaled_M_insid+scaled_A_dieg+scaled_G_littor,
    # integrate TWV, NIT, and PIT
            Raw_RBI = TWV + NIT + (2 * PIT),
    # scaling the RBI score
            score = ((Raw_RBI - 0.03)/ 4.69),
    # RBI Categories based on RBI scores
            condition_category = case_when( (score > 0.27) ~ "Reference",
                                  (score > 0.16 & score <= 0.27) ~ "Low Disturbance",
                                  (score >= 0.09 & score <= 0.16) ~ "Moderate Disturbance",
                                  (score <0.09)  ~ "High Disturbance" ),
    # RBI Category Scores based on RBI scores
            condition_category_score = case_when( (condition_category == "Reference") ~ 1,
                                          (condition_category == "Low Disturbance") ~ 2,
                                          (condition_category == "Moderate Disturbance") ~ 3,
                                          (condition_category == "High Disturbance") ~ 4),
            index = "RBI")

  #Export an interim file with the raw and scaled rbi metrics for user review
  write.csv(rbi_scaled, file = paste(output_path, "/", file_id, " SQO RBI interim 3 - scaled RBI metrics.csv ", sep=""), row.names = FALSE)


  #gathering the station information (i.e., non-taxonomic data) for each site
  rbi.stations<-BenthicData %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()

  rbi.out<-rbi_scaled %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score) %>%
    bind_rows(rbi.out.null, ., defaunated) %>%
    full_join(rbi.stations, ., by=c("stationid", "sampledate", "replicate"))

  if(length(rbi.out$stationid)>1)#if RBI scores are calculated, function will drop the dummy data placeholders and report the data for the submitted samples
  {
   rbi.out.2<-rbi.out %>%
      filter(stationid!="dummy") %>%
      mutate(note=if_else(is.na(note), "none", note))
  }

  else
  {
    rbi.out.2<-rbi.out %>% # if RBI scores are not calculated, the function will only report the dummy data placeholders
      mutate(note="RBI scores not caculated")
  }
  return(rbi.out.2)
  write.csv(rbi.out.2, paste(output_path, "/", file_id, "  SQO RBI score.csv", sep=""), row.names = FALSE)


   



}

