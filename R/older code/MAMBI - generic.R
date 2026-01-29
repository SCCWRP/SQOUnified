# Compute the multivariate AMBI (M-AMBI) index score.
#
#     This is a function to calculate the United States version of the multivariate AMBI (M-AMBI)
#     index scores following Pelletier et al. 2018, which is in turn built upon the work of
#     Sigovini et al. 2013 and Muxica et al. 2007.
#
#     The function is designed for use in US estuarine waters and requires three arguments:
#     BenthicData, EG_Ref_values, and EG_Scheme. More details are given below.
#
#     BenthicData is a data frame with that must have the following fields (additional fields related to
#     station information are okay):
#            1. stationid - an alpha-numeric identifier of the sample site
#            2. replicate - serial number identifying the number of replicate samples during a given event
#            3. sampledate - the date of sample collection in yyyy-mm-dd date format
#            4. latitude - latitude of sample location in decimal degrees
#            5. longitude - longitude of sample location in decimal degrees (remember negative sign for western hemisphere)
#            6. salinity - bottom water salinity at time of sampling in PSU
#            7. taxon - name of the fauna, ideally in SCAMIT ed14 format, do not use sp. or spp.,
#              use sp only or just the Genus. If no animals were present in the sample use
#              NoOrganismsPresent with 0 abundance
#             8. abundance - a numeric value indicating the number of individuals counted for the specified taxon from the sample
#             9. exclude - Yes or No value, indicating if the taxon name is ambiguous relative to other taxa in the sample
#       EG_Ref_Values - is a data frame with Ecological Group (i.e., tolerance values) assignements for the infauna (I, II, III, IV, or V) and
#             Oligochaeta - a Yes/No designation if the organism is an oliogichaete (important for tidal freshwater samples). The
#             default value is to NULL, which will use the accompanying EG list file upon Gillett et al. 2015 and updated with each Bight cycle
#       EG_Scheme is a quoted string specifying the name of the field in the EG_Ref_Values data frame that contains the EG values. The default
#             is "Hybrid", after Gillett et al. 2015. Other standard options include "US_East", "Us_West", "US_Gulf" coasts, "US" and "Standard"
#
#     The function requires two additional dataframes that are supplied: Saline and Tidal freshwater
#     good-bad standards for the M-AMBI. These files will be loaded from TidalFresh_Standards.RData and Saline_Standards.RData
#
#     For the function to run, the following packages NEED to be installed:  tidyverse and vegan.
#        Additionally the EQR.R function is required and will be loaded as part of the function.
#
#       NOTE: This code is designed to be run in the Bight Program's Generic SQO BLOE calculator wrapper function. However, the function can
#             be run independently but the user must load the following into their R environment:
#                 1. file_id - this is a quoted string for you to identify the data the M-AMBI scores are associated with
#                   e.g., "Bight 23" or "2024 San Diego Bay" - this will be used to name all of the output files
#                 2. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
#                 3. BenthicData - a data frame containing the benthic data and the station information, detailed above
#
# This function will produce a csv file of final M-AMBI scores for each sample, as well as csv files of interim tables detailing the taxa in the
#     submitted data that have EG values assigned and one for taxa in the submitted data that do not have an EG value assigned.
MAMBI.generic<-function(BenthicData, EG_Ref_values = NULL, EG_Scheme="Hybrid", output_path, file_id)
{

  # loading in required reference files
  load("Reference Files/Saline_Standards.RData") # Good-Bad Benchmarks following Pelletier et al. 2018
  load("Reference Files/TidalFresh_Standards.RData") #Good-Bad Benchmarks following Pelletier et al. 2018
  load("Reference Files/us mambi egvalues 04-23-24.RData") #Standard EG value schemes following and Gillett et al. 2014 and USEPA
  source("Reference Files/EQR.R")


  #loading in required r packages
  require(tidyverse)
  require(vegan)



#remnant option that allows user to support their own EG scheme, if they have one. If not, it defaults to the standard list
  if (is.null(EG_Ref_values)) {
    EG_Ref_values <- us.mambi.eg.values.04_23_24
  }




  Input_File.0 <- BenthicData %>%
    mutate(
      Species_ended_in_sp = (str_detect(taxon," sp$")),
      taxon=(str_replace(taxon, " sp$",""))) %>%
    #splitting out west coast sites to facilitate use of different criteria for west coast high salinities
    mutate(
      Coast = (ifelse(longitude<=-115,"West","Gulf-East")))%>%
    # associating a sample with one of the different salinity zones and their standards
    mutate(
      SalZone = case_when(
        salinity > 30 & salinity <= 40 & Coast == "Gulf-East" ~ "EH",
        salinity > 18 & salinity <= 30 & Coast == "Gulf-East" ~ "PH",
        salinity > 5 & salinity <= 18 ~ "MH",
        salinity > 0.2 & salinity <= 5 ~ "OH",
        salinity >= 0 & salinity <= 0.2 ~ "TF",
        salinity > 40 ~ "HH",
        salinity > 30 & salinity <= 40 & Coast == "West" ~ "WEH",
        salinity > 18 & salinity <= 30 & Coast == "West" ~ "WPH"
      )
    )



  EG_to_use <- EG_Ref_values %>%
    #something in the update of Tidyverse/dplyr prevents the usage of a vector (i.e., EG_Scheme) in a select statement
    #the appropriate syntax is pass the vector into a all_of() or any_of() statement. FWIW, rename is wrapper around select
    rename(EG=any_of(EG_Scheme)) %>%
    mutate(EG = ifelse(Taxon=="Oligochaeta", "V", EG)) %>%
    select(Taxon, Exclude, EG)



  #in case a sample had no animals (e.g., taxon=NoOrganismsPresent), we force it into the High Disturbance category.
  #the calculator would not be able to process that sample and would drop it, so we deal with it apriori
  defaunated<-Input_File.0 %>%
    filter(taxon=="NoOrganismsPresent") %>%
    mutate(mambi_score=0, Orig_mambi_condition= "Bad", SQO_mambi_condition="High Disturbance", ambi_score=7, H=NaN, S=NaN, oligo_pct=NaN,
           use_ambi="Yes", use_mambi="Yes", note="Defaunated Sample") %>%
    select(stationid, sampledate, replicate,latitude, longitude, mambi_score, Orig_mambi_condition, SQO_mambi_condition,
           ambi_score, H, S, oligo_pct, use_mambi, use_ambi, note)



# dropping azoic samples from analysis
  Input_File<-Input_File.0 %>% filter(taxon != "NoOrganismsPresent")



  Sample.info<-Input_File %>%
    select(stationid, replicate, sampledate, latitude, longitude, salinity, Coast, SalZone) %>%
    distinct()

  Input_File2<-Input_File %>%
    filter(!is.na(SalZone))

  EG.Assignment<-Input_File %>%
    left_join(., EG_to_use, by=c("taxon"="Taxon")) %>% #filter(Exclude!="Yes") #%>%
    group_by(stationid, replicate, sampledate) %>%
    mutate(tot_abun=sum(abundance)) %>%
    ungroup() %>%
    #left_join(.,total.abundance, by=c("StationID", "Replicate", "SampleDate")) %>%
    mutate(rel_abun=((abundance/tot_abun)*100))

  # Export interim file with all taxa in submitted data with an assigned EG value for user to review
  taxa_w_EG<-EG.Assignment %>%
    filter(EG%in%c("I", "II", "III", "IV", "V")) %>%
    group_by(taxon, EG) %>%
    summarise(total_abundance=sum(abundance))

  write.csv(taxa_w_EG, file=paste(output_path, "/", file_id, " M-AMBI interim 1 - taxa with EG values.csv", sep=""), row.names = FALSE)

  # Export interim file with all taxa in submitted data without an assigned EG value for user to review
  taxa_wo_EG<-EG.Assignment %>%
    filter(is.na(EG)|EG=="") %>%
    group_by(taxon) %>%
    summarise(total_abundance=sum(abundance))

  write.csv(taxa_wo_EG, file=paste(output_path, "/", file_id, " M-AMBI interim 2 - taxa without EG values.csv", sep=""), row.names = FALSE)


    # EG.Assignment.cast<-data.frame(NoEG=numeric(),
    #                              YesEG=numeric())

  # Calculating aplicability of the index to each sample based upon the % of abundance that can be assigned and EG value
  # Guidelines follow Pelleteier et al 2018, which follow Borja and Muxica 2005
  AMBI.applicability<-EG.Assignment %>%
    mutate(EG_Test=ifelse(is.na(EG),"NoEG", "YesEG")) %>%
    pivot_wider(id_cols = c(stationid, replicate, sampledate), names_from = EG_Test, values_from = rel_abun, values_fill = 0, values_fn = sum) %>%
    mutate( use_ambi = case_when(
        NoEG <= 20 ~ "Yes",
        NoEG > 20 & NoEG <= 50 ~ "With Care",
        NoEG > 50 ~ "Not Recommended",
        is.na(NoEG) ~ "Yes"))

  # Recommendations for applicability of the M-AMBI based upon AMBI suitability and presence of salinity data
  # MAMBI.applicability <- Sample.info %>%
  #   mutate(use_mambi = ifelse(is.na(SalZone),"No - No Salinity Value","Yes")) %>%
  #   select(stationid, replicate, sampledate, use_mambi)

  MAMBI.applicability<-AMBI.applicability %>%
    left_join(., Sample.info, by=c("stationid", "replicate", "sampledate")) %>%
    mutate(use_mambi=case_when(is.na(SalZone)~"No - No Salinity Value",
                               use_ambi=="With Care"~"Caution - Sparse AMBI Coverage",
                               use_ambi=="Not Reccommended"~ "Not Recommended - Poor AMBI Coverage",
                               TRUE~"Yes")) %>%
    select(stationid, replicate, sampledate, use_ambi, use_mambi)



# establish the salinity zones in the data set
  # used later to identify which sub-routine and metric strandards to use
  Sal_range.dataset<-unique(Input_File2$SalZone)


  ######Saline calcs ################
  # if a sample is from oligohaline through euhaline salinities, M-AMBI is calculated using AMBI, Species Richness, and Species Diversity

  # Calculating AMBI for each sample
  AMBI.Scores<-EG.Assignment %>%
    group_by(stationid, replicate, sampledate,tot_abun,EG) %>%
    summarise(sum_rel=sum(rel_abun)) %>%
    ungroup() %>%
    replace_na(list(EG="NoEG")) %>%
    mutate(
      EG_Score = case_when(
        EG == "I" ~ sum_rel*0,
        EG == "II" ~ sum_rel*1.5,
        EG == "III" ~ sum_rel*3,
        EG == "IV" ~ sum_rel*4.5,
        EG == "V" ~ sum_rel*6,
        EG == "NoEG" ~ 0)) %>%
    mutate(EG_Score=ifelse(tot_abun==0,7,EG_Score)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(ambi_score=(sum(EG_Score, na.rm=TRUE)/100)) %>%
    ungroup()

 #Calculating taxa richness for each sample
   Rich<-Input_File %>%
    filter(exclude!="Yes") %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(S=length(taxon)) %>%
    ungroup()

  #Rich$S<-as.numeric(Rich$S)


   # Claculating base 2 shannon wiener diversity for each sample using vegan package
   Divy<-Input_File %>%
     filter(exclude!="Yes") %>%
     pivot_wider(id_cols=c(stationid, replicate, sampledate), names_from = taxon, values_from=abundance, values_fill = 0) %>%
     mutate(H=vegan::diversity(select(., 4:ncol(.)), index="shannon", base=2)) %>%
     select(stationid, replicate, sampledate, H)


  # Combining the three metrics
   metrics<-AMBI.Scores %>%
    left_join(.,Rich, by=c("stationid", "replicate", "sampledate")) %>%
    left_join(.,Divy, by=c("stationid", "replicate", "sampledate"))

  # Combining the metrics with the station information for each sample
   metrics.1 <- Sample.info %>%
    left_join(., metrics, by=c("stationid", "replicate", "sampledate")) %>%
    select(stationid, replicate, sampledate,ambi_score, S, H, SalZone)

  # Isolating samples for which salinity has not been submitted
   no.SalZone.data<- metrics.1 %>%
     filter(is.na(SalZone)) %>%
    left_join(.,Sample.info, by=c("stationid", "replicate", "sampledate", "SalZone")) %>%
    select(stationid, replicate, sampledate, ambi_score, S, H, latitude, longitude, SalZone)

  # Selecting those samples for which salinity was submitted and adding on the standard good/bad endpoints for each salinity zone
   metrics.2<- metrics.1 %>%
    filter(!is.na(SalZone)) %>%
    bind_rows(.,Saline_Standards)

  # sumbit metrics into the saline MAMBI calculator factor analysis function for all non-tidal freshwater samples
  #  follows from Pelletier et al. 2018

   # written as a purrr::map statement calling and creating a new function "sal" to run the factor analysis across the different salinity zones
   # note the use of map at the beginning and list_rbind at the end to mimic behavior of map_df (superceded) while maintaining practices of purrr v1.0.0
  saline.mambi<-purrr::map(Sal_range.dataset, function(sal)
  {
    sal.df<- filter(metrics.2, SalZone == sal) #iterating through each salinity zone present in the samples
    METRICS.tot<-select(sal.df, ambi_score, S, H) #islolating the metrics to be fed into the factor analysis

    # doing the factor analysis among the data submitted and the Pelletier et al. 2018 standards for best and worse conditions in a given salinity zone
    options(warn = -1)
    METRICS.fa2 <- princomp(METRICS.tot, cor = T, covmat = cov(METRICS.tot))
    options(warn = 0)
    METRICS.fa2.load <- loadings(METRICS.fa2) %*% diag(METRICS.fa2$sdev)
    METRICS.fa2.load.varimax <- loadings(varimax(METRICS.fa2.load))
    METRICS.scores2 <- scale(METRICS.tot) %*% METRICS.fa2.load.varimax
    colnames(METRICS.scores2) <- c("x", "y", "z")
    METRICS.tr <- METRICS.scores2


    # Transforming factor analysis 3-axis loadings into M-AMBI scores using the EQR function
    eqr <-EQR(METRICS.tr)
    colnames(eqr)<-c("mambi_score")
    eqr<-data.frame(eqr)

    # Adding back in station information and classifying scores into traditional M-AMBI (5) and SQO (4) condition categories
    results<-sal.df %>%
      bind_cols(.,eqr) %>%
      left_join(.,Sample.info, by=c("stationid", "replicate", "sampledate", "SalZone")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, SalZone, ambi_score, S, H, mambi_score) %>%
      filter(!stationid%in%Saline_Standards$stationid, SalZone!="TF") %>% #dropping the the standards, leaving only submitted data
      mutate(
        Orig_mambi_condition = case_when(
          mambi_score < 0.2 ~ "Bad",
          mambi_score >= 0.2 & mambi_score < 0.39 ~ "Poor",
          mambi_score >= 0.39 & mambi_score < 0.53 ~ "Moderate",
          mambi_score >= 0.53 & mambi_score < 0.77 ~ "Good",
          mambi_score >= 0.77 ~ "High"
        ),
        SQO_mambi_condition = case_when(
          mambi_score <= 0.387 ~ "High Disturbance",
          mambi_score > 0.387 & mambi_score < 0.483 ~ "Moderate Disturbance",
          mambi_score >= 0.483 & mambi_score < 0.578 ~ "Low Disturbance",
          mambi_score >= 0.578 ~ "Reference")
      )

  }
  )%>% list_rbind() # Combine the results from each SalZone/iteration of the function into one data frame






  ###################
  # if there are tidal freshwater samples in the submitted data, they need to be seperated out and run thorugh the TF subroutine and then joined to
  # the saline results

  if(any(Sal_range.dataset=="TF"))
  {

    TF.EG.Assignment <- EG.Assignment %>% filter(SalZone=="TF")
    TF.EG_Ref_values <- us.mambi.eg.values.04_23_24 %>%
      select(.,Taxon, Exclude, EG=any_of(EG_Scheme), Oligochaeta)

   # calculate AMBI scores for each sample
     TF.AMBI.Scores <- TF.EG.Assignment %>%
      group_by(stationid, replicate, sampledate, tot_abun, EG, rel_abun) %>%
      summarise(sum_rel=sum(rel_abun)) %>%
      ungroup() %>%
      replace_na(list(EG="NoEG")) %>%
      mutate(
        EG_Score = case_when(
          EG == "I" ~ sum_rel*0,
          EG == "II" ~ sum_rel*1.5,
          EG == "III" ~ sum_rel*3,
          EG == "IV" ~ sum_rel*4.5,
          EG == "V" ~ sum_rel*6,
          EG == "NoEG" ~ 0)) %>%
      group_by(stationid, replicate, sampledate) %>%
      summarise(ambi_score = sum(EG_Score)/100) %>%
      ungroup()

    # calculate % oligochaetes in each sample
     TF.Oligos <- Input_File %>%
      group_by(stationid, replicate, sampledate) %>%
      mutate(tot_abun=sum(abundance)) %>%
      ungroup() %>%
      left_join(., TF.EG_Ref_values, by=c("taxon"="Taxon") ) %>%
      filter(Oligochaeta=="Yes", SalZone=="TF") %>%
      group_by(stationid, replicate, sampledate) %>%
      summarise(oligo_pct = sum(abundance/tot_abun) * 100) %>%
      ungroup()


    # calculate shannon wiener diversity (base 2) for each sample
     TF.Divy <- Input_File %>%
      filter(SalZone=="TF") %>%
      pivot_wider(id_cols = c(stationid, replicate, sampledate), names_from=taxon, values_from=abundance, values_fill = 0) %>%
      mutate(H = diversity((select(.,4:(ncol(.)))), index = "shannon", base = 2)) %>%
      select(.,stationid, replicate, sampledate, H)


    # join the three metrics together
     TF.metrics <- TF.AMBI.Scores %>%
      left_join(.,TF.Divy, by=c("stationid", "replicate", "sampledate")) %>%
      left_join(.,TF.Oligos, by=c("stationid", "replicate", "sampledate"))

    # add on sample information
     TF.metrics.1<-Sample.info %>%
      filter(SalZone=="TF") %>%
      left_join(., TF.metrics, by=c("stationid", "replicate", "sampledate")) %>%
      select(stationid, replicate, sampledate, ambi_score, H, oligo_pct, SalZone)

    # add the best and worst standard values for the tidal freshwater salinity zone
     TF.metrics.2<-bind_rows(TF.metrics.1, TidalFresh_Standards)

    TF.METRICS.tot<-TF.metrics.2 %>%
      select(ambi_score, H, oligo_pct)


    # sumbit metrics into the tidal freshwater MAMBI calculator factor analysis function for all non-tidal freshwater samples
    #  follows from Pelletier et al. 2018

    # doing the factor analysis among the data submitted and the Pelletier et al. 2018 standards for best and worse conditions in a given salinity zone
    options(warn = -1)
    TF.METRICS.fa2 <- princomp (TF.METRICS.tot, cor = T, covmat = cov(TF.METRICS.tot))
    options(warn = 0)
    TF.METRICS.fa2.load <- loadings(TF.METRICS.fa2) %*% diag(TF.METRICS.fa2$sdev)
    TF.METRICS.fa2.load.varimax <- loadings(varimax(TF.METRICS.fa2.load))
    TF.METRICS.scores2 <- scale(TF.METRICS.tot) %*% TF.METRICS.fa2.load.varimax
    colnames(TF.METRICS.scores2) <- c("x", "y", "z")
    TF.METRICS.tr <- TF.METRICS.scores2

    # Transforming factor analysis 3-axis loadings into M-AMBI scores using the EQR function
    TF.eqr <-EQR(TF.METRICS.tr)
    colnames(TF.eqr) <- c("mambi_score")
    TF.eqr <- data.frame(TF.eqr)

    # Adding back in station information and classifying scores into traditional M-AMBI (5) and SQO (4) condition categories
     TF.mambi <- TF.metrics.2 %>%
      bind_cols(.,TF.eqr) %>%
      left_join(.,Sample.info, by=c("stationid", "replicate", "SalZone", "sampledate")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, ambi_score, H, oligo_pct, mambi_score) %>%
      filter(!stationid%in%TidalFresh_Standards$stationid) %>%
      mutate(
        Orig_mambi_condition = case_when(
          mambi_score < 0.2 ~ "Bad",
          mambi_score >= 0.2 & mambi_score < 0.39 ~ "Poor",
          mambi_score >= 0.39 & mambi_score < 0.53 ~ "Moderate",
          mambi_score >= 0.53 & mambi_score < 0.77 ~ "Good",
          mambi_score >= 0.77 ~ "High"),
        SQO_mambi_condition = case_when(
          mambi_score <= 0.387 ~ "High Disturbance",
          mambi_score > 0.387 & mambi_score < 0.483 ~ "Moderate Disturbance",
          mambi_score >= 0.483 & mambi_score < 0.578 ~ "Low Disturbance",
          mambi_score >= 0.578 ~ "Reference"))

    # adding empty columns for %oligochaetes and taxa richness to saline and TF mambi output tables respectively so they line up when joined
    saline.mambi.2<- saline.mambi %>%
      mutate(oligo_pct=NA) %>%
      select(stationid, replicate, sampledate, latitude, longitude, ambi_score, S, H, oligo_pct, mambi_score,
             Orig_mambi_condition, SQO_mambi_condition, SalZone)

    TF.mambi.2<-TF.mambi %>%
      mutate(S=NA, SalZone="TF") %>%
      select(stationid, replicate, sampledate, latitude, longitude, ambi_score, S, H, oligo_pct, mambi_score,
             Orig_mambi_condition, SQO_mambi_condition, SalZone)


    Overall.Results <- bind_rows(saline.mambi.2, TF.mambi.2, no.SalZone.data) %>%
      left_join(., MAMBI.applicability, by=c("stationid", "replicate", "sampledate")) %>%
      mutate(note=case_when(is.na(SalZone)~"No salnity data - cannot calculate M-AMBI",
                            use_mambi=="Yes"~"None",
                            use_mambi=="Caution - Sparse AMBI Coverage"~ "Check sample - M-AMBI ok, but limited taxa coverage",
                            use_mambi=="Not Recommended - Poor AMBI Coverage"~ "M-AMBI not recommended - poor taxa coverage")) %>%

      select(stationid, replicate, sampledate, latitude, longitude, mambi_score, Orig_mambi_condition, SQO_mambi_condition,
             ambi_score, H, S, oligo_pct, use_mambi, use_ambi, note) %>%

      bind_rows(.,defaunated) %>%
      mutate(index = "M-AMBI", .before=latitude)
  }


  else
  {
    #if there are no tidal freshwater samples, only the saline calculations are needed
    Overall.Results<-saline.mambi %>%
      bind_rows(.,no.SalZone.data) %>%
      left_join(., MAMBI.applicability, by=c("stationid", "replicate", "sampledate")) %>%
      mutate(oligo_pct=NaN,
             note=case_when(is.na(SalZone)~"No salnity data - cannot calculate M-AMBI",
                            use_mambi=="Yes"~"None",
                            use_mambi=="Caution - Sparse AMBI Coverage"~ "Check sample - M-AMBI ok, but limited taxa coverage",
                            use_mambi=="Not Recommended - Poor AMBI Coverage"~ "M-AMBI not recommended - poor taxa coverage")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, mambi_score, Orig_mambi_condition, SQO_mambi_condition,
             ambi_score, H, S, oligo_pct, use_mambi, use_ambi, note) %>%

      bind_rows(.,defaunated) %>%
      mutate(index = "M-AMBI", .before=latitude)
  }

  #gathering all of the site/sample information that was initially submitted with the benthic data to attach to the m-ambi scores
  station.info<-BenthicData %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()


  Overall.Results.2<-Overall.Results %>%
    left_join(., station.info, by=c("stationid", "sampledate", "replicate", "latitude", "longitude"))

  return(Overall.Results.2)
  write.csv(Overall.Results.2, file=paste(output_path, "/", file_id, " M-AMBI Scores.csv", sep=""), row.names = FALSE)
}



