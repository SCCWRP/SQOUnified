#Function designed to input benthic data for SQO calculation, calculate SQO Benthic Line of Evidence (BLOE), calculate US M-AMBI,
#and produce both interim and final output files

#This function requires that you have the following packages installed on your machine:
#       1. tidyverse
#       2. naniar
#       3. openxlsx
#       4. vegan
#
# in your r environments you can install packages with:  install.packages(".......")

##################################################################################################################################
# Instructions for use
# The function has 3 required inputs:
#      1. file_id - this is a quoted string for you to identify the data the BRI scores are associated with
#           e.g., "Bight 23" or "2024 San Diego Bay" - this will be used to name all of the output files
#      2. infauna_path - a quoted string detailing the name and location of the .csv file with infauna abundance
#         data. Remember with path names, that R use forward slashes "/" not the normal backslash "\" that windows typically uses.
#         All column names should be lower case. The expectation is that the file will have, at a minimum, a column for each of:
#           stationid - unique identifier for that station, preferably formatted as character/text value
#           sampledate - the date on which the sample was collected, must be in a mm/dd/yyyy format
#           replicate - a number identifying the replicate infauna sample collected from the specified station on the specified date
#           taxon - character string identifying the organism. naming conventions should follow SCAMIT edition 12
#           abundance - a numeric value indicating the number of individuals counted for the specified taxon from the sample
#           exclude - Yes or No value, indicating if the taxon name is ambiguous relative to other taxa in the sample
#           depth - station depth in meters
#           latitude - station latitude in decimal degrees
#           longitude - station longitude in decimal degrees (negative values for west longitudes)
#           salinity - bottom salinity at time of sample collection in PSU
#       3. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
###########################################################  ########################################################################


SQO_BLOE.generic<-function(file_id, infauna_path, output_path)
{
  #loading required packages and the individual functions to calculate each of the individual SQO benthic indices
  require(tidyverse)
  require(naniar)
  require(openxlsx)
  require(vegan)
  source("R/SQO BRI - generic.R")
  source("R/IBI - generic.R")
  source("R/RBI - generic.R")
  source("R/RIVPACS - generic.R")
  source("R/MAMBI - generic.R")

  print(paste("saving files to ", output_path, sep=""))

  #cleaning up the infauna data to be analyzed
  #standardizing the nomenclature of ommitted salinity or depth data to NA, so that it can be dealt with later
  BenthicData<-read.csv(infauna_path) %>%
  replace_with_na(., list(salinity=c(-88, -99), depth=c(-88, -99))) %>%
    mutate(sampledate=mdy(sampledate))

  #running each of the individual benthic indices

  bri.scores.x<-SQO.BRI.generic(BenthicData, output_path, file_id)

  ibi.scores.x<-IBI.generic(BenthicData, output_path, file_id)

  rbi.scores.x<-RBI.generic(BenthicData, output_path, file_id)

  rivpacs.scores.x<-RIVPACS.generic(BenthicData, output_path, file_id)

  mambi.scores.x<-MAMBI.generic(BenthicData, EG_Ref_values = NULL, EG_Scheme = "Hybrid", output_path, file_id)



  # aggregating all of the individual index scores so they can be reviewed by the user and used to calculate the BLOE score
  all.sqo.scores.x<-bind_rows(bri.scores.x, ibi.scores.x, rbi.scores.x, rivpacs.scores.x)


  write.csv(all.sqo.scores.x, file=paste(output_path, "/", file_id, "SQO BLOE interim - all benthic indices scores.csv", sep=""), row.names = FALSE)


  BLOE.scores.x<-all.sqo.scores.x %>%

    group_by(stationid, sampledate, replicate) %>%
    #Integrating the four seperate indices by calculating the median of the scores and rounding up to the next integer
    mutate(BLOE_score=ceiling(median(condition_category_score, na.remove=TRUE)),
              Note=str_flatten(note, collapse=",")) %>%
    ungroup() %>%
    pivot_wider(id_cols=c(-score, -condition_category,  -note),
                names_from = index, values_from = condition_category_score) %>%
    # The notes column provides the user some context about the scores. Particularly about the appropriateness of the BLOE indices relative to salinity at time of sampling
    # It will also carry over any notes generated during calculation of the individual indices
    mutate(notes=case_when(is.na(salinity)&is.na(Note)~"Salinity value unknown-confirm salinity >=27 PSU for BLOE to be accurate.",
                           is.na(salinity) & !is.na(Note)~paste("Salinity value unknown-confirm salinity >=27 PSU for BLOE to be accurate.", Note, sep="; "),
                           salinity<27 & is.na(Note)~"Caution, salinity value <27 PSU - BLOE may not be accurate. Consider M-AMBI",
                           salinity<27 & !is.na(Note)~paste("Caution, salinity value <27 PSU - BLOE may not be accurate. Consider M-AMBI", Note, sep="; "),
                           salinity>=27 & !is.na(Note)~Note,
                           salinity>=27 & is.na(Note)~"None",
                           TRUE~"bob"),
           BLOE_category=case_when(BLOE_score==1~"Reference",
                                    BLOE_score==2~"Low Disturbance",
                                    BLOE_score==3~"Moderate Disturbance",
                                    BLOE_score==4~"High Disturbance")
           ) %>%
    relocate(stationid, sampledate, replicate, BLOE_score, BLOE_category, BRI_cond=BRI, IBI_cond=IBI, RBI_cond=RBI, Rivpacs_cond=Rivpacs) %>%
    select(-Note)

  mambi.scores.sqoformat<-mambi.scores.x %>% #adjusting MAMBI output table to match the SQO BLOE indices format
      mutate(MAMBI_cond=case_when(SQO_mambi_condition=="Reference"~1,
                                SQO_mambi_condition=="Low Disturbance"~2,
                                SQO_mambi_condition=="Moderate Disturbance"~3,
                                SQO_mambi_condition=="High Disturbance"~4,
                                TRUE~NA), .after=SQO_mambi_condition)
  BLOE.scores.w.MAMBI<-BLOE.scores.x %>%
    left_join(., select(mambi.scores.sqoformat, stationid, sampledate, replicate, MAMBI_cond, note ), by=c("stationid", "sampledate", "replicate")) %>%
    mutate(notes=case_when(note=="None"~notes,
                             is.na(note)~notes,
                             TRUE~paste(notes, note, sep="; "))) %>%
    select(-note) %>%
    relocate(MAMBI_cond, .after=Rivpacs_cond)

  write.csv(BLOE.scores.x, file=paste(output_path, "/", file_id, " SQO integrated BLOE scores.csv", sep=""), row.names = FALSE)

  write.csv(BLOE.scores.w.MAMBI, file=paste(output_path, "/", file_id, " SQO integrated BLOE scores plus M-AMBI scores.csv", sep=""), row.names = FALSE)

  #exporting everything into a single excel workbook with the BLOE scores, but also a tab for each individual index

  #define sheet names for each data frame
  # as M-AMBI is not included in the calculation of the traditional BLOE it is provided on a separate tab to avoid confusion.

  dataset_names <- list("sqo_bloe"=BLOE.scores.x, "sqo_bri"=bri.scores.x, "ibi"=ibi.scores.x, "rbi"=rbi.scores.x, "rivpacs"=rivpacs.scores.x,
                        "sqo_bloe_w_mambi" =BLOE.scores.w.MAMBI, "mambi"=mambi.scores.sqoformat)

  #export each data frame to separate sheets in same Excel file
  write.xlsx(dataset_names, file = paste(output_path, "/", file_id, " SQO output summary.xlsx", sep=""))

return(BLOE.scores.w.MAMBI)

}
