#Function designed to input benthic data for SQO calculation, calculate SQO Benthic Line of Evidence (BLOE), calculate US M-AMBI,
#and produce both interim and final output files

#This function requires that you have the following packages installed on your machine:
#       1. tidyverse
#       2. naniar
#       3. openxlsx
#       4. vegan
#
# in your r environments you can install packages with:  install.packages(".......")

####################################################################################################################################
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
### IMPORTANT: this "alt" version of the code transforms modern taxonomy to ed5 taxonomy to sidestep DJG's hypothesized diversity/richness problem arising from
###             the improving taxonomy of modern times. the indices are scaled to a certain level of diversity and fully adopting new, more specios species
###             lists will alter the index scores to a degree. This version can be compared to a modern version to test this effect


alt.SQO_BLOE.generic<-function(file_id, infauna_path, output_path)
{
  #loading required packages and the individual functions to calculate each of the individual SQO benthic indices
  require(tidyverse)
  require(naniar)
  require(openxlsx)
  require(vegan)
  source("alternate species list approach/alt SQO BRI - generic.R")
  source("alternate species list approach/alt IBI - generic.R")
  source("alternate species list approach/alt RBI - generic.R")
  source("alternate species list approach/alt RIVPACS - generic.R")
  source("alternate species list approach/alt MAMBI - generic.R")
  load("Reference Files/ed14_to_SQO.RData")#loading in files to convert from Ed14 taxonomy to SQO compatible taxonomy

  print(paste("saving files to ", output_path, sep=""))

  #cleaning up the infauna data to be analyzed
  #standardizing the nomenclature of ommitted salinity or depth data to NA, so that it can be dealt with later
  BenthicData<-read.csv(infauna_path) %>%
  replace_with_na(., list(salinity=c(-88, -99), depth=c(-88, -99))) %>%
    mutate(sampledate=mdy(sampledate))

#### retrofitting modern taxonomy back and taking daughter taxa and combining them up where appropriate ####


  #Step 1 is to take one to one changes and swap the new names back to their old version
  one.to.one<-SoCal.SQO.ed14.link %>%
    select(Original.SQO.Taxon, Ed.14.Taxon, change, type) %>%
    filter(type%in%c("one-to-one", "spelling", "worms", "removed"))

  BenthicData.0<-BenthicData %>%
    left_join(., one.to.one, by=c("taxon"="Ed.14.Taxon")) %>%
    relocate(Original.SQO.Taxon, .before = taxon) %>%
    mutate(taxon.2=if_else(is.na(Original.SQO.Taxon), taxon, Original.SQO.Taxon),
           change_type=if_else(is.na(type), "",type),
           taxa_changed=if_else(is.na(Original.SQO.Taxon),"",taxon)) %>%

    select(-taxon, -type, -change, -Original.SQO.Taxon)




  #Step two is to combine daughter taxa to the higher taxonomic level the indices recognize
  BenthicData.1<-BenthicData.0 %>%
    relocate(abundance,taxon.2, .before = stationid) %>%
    left_join(., ed14.rollups, by=c("taxon.2"="ed.14_daughters")) %>%
    mutate(taxon.3=if_else(is.na(sqo.name), taxon.2, sqo.name),
           rolled_up=if_else(is.na(sqo.name), "no", "yes"), .before = taxon.2) %>%
    #mutate(flag=if_else(taxon==sqo.name, 1,0), .before=abundance) %>%
    select(-sqo.name) %>%
    group_by(stationid, replicate, sampledate, taxon.3) %>%
    mutate(abundance.2=sum(abundance),
           #when rolling up taxa we may be combining taxa w/ different exclude codes into one. To keep the most the taxa richness/exclude coding
           # if any of the combined taxa has a No (i.e., it is unique), we give that prescedance to maintain the richness
           all.excludes=str_flatten_comma(exclude),
           exclude.prob.flag=case_when(str_detect(all.excludes, "No, Yes")~1,
                                       str_detect(all.excludes, "Yes, No")~1,
                                       TRUE~0),
           exclude.2=case_when(rolled_up=="yes" & exclude.prob.flag==1 ~ "No",
                               TRUE~exclude),
           .before=abundance) %>%
    ungroup() %>%
    mutate(change_type.2=if_else(rolled_up=="yes", paste("rolled up to ",join_level, sep="" ),change_type),
           taxa_changed.2=if_else(rolled_up=="yes", taxon.2, taxa_changed),
           .before = stationid) %>%
    select(-all.excludes, -exclude.prob.flag, -exclude) %>%
    rename(exclude=exclude.2)


## Step 3 is to deal with complex changes where a single ed 14 taxon used to be multiple taxa on the SQO lu list
## The taxa to choose among the old SQO names was prioritized to select a taxon that had a tolerance score, sensitivity designation, etc vesus one that didn't
  ed.14.complex.2<-ed.14.complex %>%
    filter(Priority=="yes")

  BenthicData.2<-BenthicData.1 %>%
    left_join(., ed.14.complex.2, by=c("taxon.3"="Ed.14.Taxon")) %>%
    mutate(taxon.4=case_when(Priority=="yes"~ Original.SQO.Taxon,
                             TRUE~taxon.3), .before = taxon.3) %>%
    mutate(change_type.3=case_when(Priority=="yes"~"complex",
                                   TRUE~change_type.2), .before = change_type.2)

#Step 4 is to create an interim output file detailing all of the changes made to the submitted data

    taxa.changes<-BenthicData.2 %>%
    group_by(stationid, replicate, sampledate, taxon.4, change_type.3) %>%
  mutate(
         taxa_changed.3= str_flatten_comma(taxa_changed.2),
         .before=abundance) %>%
    select(stationid, sampledate, replicate, taxon_used=taxon.4, taxon_submitted=taxon.2, changed_taxa=taxa_changed.3,
           abundance_used=abundance.2,
           abundance_submitted=abundance, type_of_change=change_type.3)


  write.csv(taxa.changes, paste(output_path, "/", file_id," interim data prep file - ed14 to SQO LU list taxa changes.csv", sep=""), row.names = FALSE)

# Step 5 is to create final file to run through the different indices.
#Modern names have been changed to older, SQO-valid names and daughter taxa have been collapsed to the appropriate higher level
#Where taxa were collapsed, their abundances were summed

  BenthicData.3<-BenthicData.2 %>%
   select(-abundance, -taxon.2, -taxon.3, -rolled_up, -change, -change_type, -change_type.2, -change_type.3,
          -taxa_changed, -taxa_changed.2, -join_level, -Original.SQO.Taxon, -change, -type, -Priority) %>%
   distinct() %>%
   rename(taxon=taxon.4, abundance=abundance.2) %>%
   relocate(taxon, abundance, .after=sampledate)



  #### running each of the individual benthic indices on the newly created "retrofit" data ####
  #Each of these individual functions are available in the repository to inspected by the user, as they may like

  bri.scores.x<-alt.SQO.BRI.generic(BenthicData.3, output_path, file_id)

  ibi.scores.x<-alt.IBI.generic(BenthicData.3, output_path, file_id)

  rbi.scores.x<-alt.RBI.generic(BenthicData.3, output_path, file_id)

  rivpacs.scores.x<-alt.RIVPACS.generic(BenthicData.3, output_path, file_id)

  mambi.scores.x<-alt.MAMBI.generic(BenthicData, EG_Ref_values = NULL, EG_Scheme = "Hybrid", output_path, file_id)
  #Note that the MAMBI script uses modern taxonomy, as it can be updated unlike the SQO calculations



  # aggregating all of the individual index scores so they can be reviewed by the user and used to calculate the BLOE score
  all.sqo.scores.x<-bind_rows(bri.scores.x, ibi.scores.x, rbi.scores.x, rivpacs.scores.x) %>%
    arrange(., stationid, sampledate, replicate, index)


  write.csv(all.sqo.scores.x, file=paste(output_path, "/", file_id, " SQO BLOE interim - all benthic indices scores.csv", sep=""), row.names = FALSE)


  BLOE.scores.x<-all.sqo.scores.x %>%

    group_by(stationid, sampledate, replicate) %>%
    #Integrating the four separate indices by calculating the median of the scores and rounding up to the next integer
    mutate(BLOE_score=ceiling(median(condition_category_score, na.rm=TRUE)),
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
    left_join(., select(mambi.scores.sqoformat, stationid, sampledate, replicate, SQO_mambi_condition, MAMBI_cond, note ), by=c("stationid", "sampledate", "replicate")) %>%
    mutate(notes=case_when(note=="None"~notes,
                             is.na(note)~notes,
                             TRUE~paste(notes, note, sep="; "))) %>%
    select(-note) %>%
    rename(MAMBI_SQO_Cat=SQO_mambi_condition) %>%
    relocate(MAMBI_SQO_Cat, MAMBI_cond, .after=Rivpacs_cond)

  write.csv(BLOE.scores.x, file=paste(output_path, "/", file_id, " SQO integrated BLOE scores.csv", sep=""), row.names = FALSE)

  write.csv(BLOE.scores.w.MAMBI, file=paste(output_path, "/", file_id, " SQO integrated BLOE scores plus M-AMBI scores.csv", sep=""), row.names = FALSE)

  #exporting everything into a single excel workbook with the BLOE scores, but also a tab for each individual index

  #define sheet names for each data frame
  # as M-AMBI is not included in the calculation of the traditional BLOE it is provided on a separate tab to avoid confusion.

  dataset_names <- list("sqo_bloe"=BLOE.scores.x, "sqo_bri"=bri.scores.x, "ibi"=ibi.scores.x, "rbi"=rbi.scores.x, "rivpacs"=rivpacs.scores.x,
                        "sqo_bloe_w_mambi" =BLOE.scores.w.MAMBI, "mambi"=mambi.scores.sqoformat)

  #export each data frame to separate sheets in same Excel file
  write.xlsx(dataset_names, file = paste(output_path, "/", file_id, " SQO BLOE output summary.xlsx", sep=""))

  #putting the BLOE+MAMBI final scores into the R environment for the user to view
return(BLOE.scores.w.MAMBI)

}
