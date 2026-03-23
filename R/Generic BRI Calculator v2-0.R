#Generic BRI calculator v2.0

#Function designed to input benthic data and calculate the offshore BRI of Smith et al. (2001)
# it will produce an excel workbook and individual csv worksheets containing the BRI scores, as well as
# a number of interim files so that the user can review how the index works and how it processes the data the
# user submitted

#V2.0 calculates BRI scores for samples from the overlapping depth zones by calculating a score from each of the two depth
#zones and averaging them, as opposed to the creation and usage of overlap p-values like previous versions. This decision was based
#upon the Bight 2023 index code sub-committee's discussions. In response to the index code committee's feedback, this version also
#produces two copies of the calculator outputs:
# #1 is an excel workbook containing the final scores and all of the interim files compiled together; #2 is a list of r objects/dataframes
#in your local r environment so that you can call them (i.e., objectname$taxa_w_o_pcode) and work on them without having to
#load in the csv/excel file.
#Additionally, and also based on committee feedback, this version will
# 1. produce a warning in the notes column if the sample is located north of Point Conception or south of the US-Mexico border since the
#BRI is nominally only applicable to the Southern California Bight.
# 2. the interim table 3 - taxa that were submitted but couldn't be matched with code now provides some suggestions that account for
#potential mis-spelling of the specimen's name in the data file. If it is a typo, then it is up to the user to correct their input
#file and re-run the calculator. the "example infauna.csv" file now has some intentional misspellings to highlight/test this functionality

#This function requires that you have the following packages installed on your machine:
#       1. tidyverse
#       2. openxlsx
#       3. janitor
#       4. messydates
#
# in your r environments you can install packages with:  install.packages(".......")




####################################################################_#############################################################
# Instructions for use
# The function has 4 required inputs:
#      1. file_id - this is a quoted string for you to identify the data the BRI scores are associated with
#           e.g., "Bight 23" or "2024 southern shelf" - this will be used to name all of the output files
#      2. infauna_path - a quoted string detailing the name and location of the .csv file with infauna abundance
#         data. Remember with path names, that R use forward slashes "/" not the normal backslash "\" that windows typically uses.
#         All column names should be lower case. The expectation is that the file will have, at a minimum, a column for each of:
#           stationid - unique identifier for that station, preferably formatted as character/text value
#           sampledate - the date on which the sample was collected, must be in a mm/dd/yyyy format
#           replicate - a number identifying the replicate infauna sample collected from the specified station on the specified date
#           taxon - character string identifying the organism. naming conventions should follow SCAMIT edition 14
#           abundance - a numeric value indicating the number of individuals counted for the specified taxon from the sample
#      3. station_path - a quoted string detailing the name and location of the .csv file with station information for each infauna
#           sample. All column names should be lower case. The expectation is that the file  will have, at a minimum,a column for:
#           stationid - unique identifier for that station, preferably as character/text value
#           depth - station depth in meters
#           latitude - station latitude in decimal degrees
#           longitude - station longitude in decimal degrees (negative values for west longitudes)
#       4. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
 ###################################################################################################################################


offshore_bri_calc_ed14_v2<-function(file_id, infauna_path, station_path, output_path)
  {
require(tidyverse)
require(openxlsx)
require(janitor)
require(messydates)

####Input files
infauna <- read.csv(infauna_path) #the user's infauna to be submitted
station_info<-read.csv(station_path) #the user's station information to be submitted

print(paste("Saving files to ", output_path, sep=""))

load("Reference Files/pcode_14.RData") #pcode values by depth and taxa associated with pcodes

###            Prep the data ####

#ensuring data are in the correct formats

#think about adding a latitude gaurdrails


station_info.2<-station_info %>%
  clean_names(., case="snake" , sep_out="") %>% #trying to be a little flexible incase odd field names are used
  mutate(stationid=as.character(stationid),
         depth=as.numeric(depth))


infauna.2<-infauna %>%
  clean_names(., case="snake" , sep_out="") %>% #trying to be a little flexible incase odd field names are used
  mutate(stationid=as.character(stationid),
         abundance=as.numeric(abundance),
         sampledate=as_messydate(sampledate),#trying to be flexible with date format given the wierdness of excel vs .CSV vs. R
         sampledate=as_date(sampledate),
         sample_id=paste(stationid, sampledate, replicate, sep="_"))

#joining taxa names, and counts to depth and geographic information by the stationid field
taxa_to_calc<-infauna.2 %>% select(sample_id, stationid, replicate,sampledate,taxon, abundance) %>%
  left_join(., select(station_info.2,stationid, depth, latitude, longitude), by="stationid") %>%
  arrange(sample_id, desc(abundance))

#output the joined taxa-station info for review
write.csv(taxa_to_calc, paste(output_path, "/", file_id, " interim file 1 - taxa to be analyzed.csv", sep=""), row.names = FALSE)


all.4.bri<-taxa_to_calc %>%
  left_join(., ed.14.ptaxa, by="taxon") %>%
  left_join(., pcodes, by="p_code") %>%
  select(-shallow_mid, -mid_deep)


# output taxa with their p-codes for review
with.pcode<-all.4.bri %>% distinct(taxon, p_code) %>% drop_na(p_code) %>% arrange(taxon)

write.csv(with.pcode, paste(output_path, "/", file_id, " interim file 2 - taxa w assigned pcodes.csv", sep=""), row.names = FALSE)
#


# output taxa without p-codes for review
no.pcode<-all.4.bri %>% distinct(taxon, p_code) %>%
  filter(is.na(p_code)) %>%
  select(-p_code) %>%
  arrange(taxon) %>%
  fuzzyjoin::stringdist_left_join(., ed.14.ptaxa,  by=c("taxon")) %>%
  select(submitted_taxon=taxon.x, did_you_mean_taxon=taxon.y, p_code)


write.csv(no.pcode, paste(output_path, "/", file_id, " interim file 3 - taxa w-o pcodes.csv", sep=""), row.names = FALSE)

write.csv(all.4.bri, paste(output_path, "/", file_id, " interim file 4 - taxa and pcodes by sample.csv",sep=""), row.names=FALSE)





####          Calculate Scores

#pulling out defaunated samples and classifying them as Defaunation and assigning a score of 5 (i.e., the worst)
defaunated<-all.4.bri %>%
  filter(taxon=="NoOrganismsPresent") %>%
  select(sample_id, stationid, sampledate, replicate) %>%
  mutate(tot_bri_abun_dz=0,

         bri_score=NaN,
         bri_cond="Defaunation",
         bri_class=5)

defaunated.samps<-unique(defaunated$sample_id)


#pulling out the samples outside the range of the index. They get no scores nor condition classification
too.deep<-all.4.bri %>%
  filter(depth>324) %>%
  distinct(sample_id, stationid, sampledate, replicate) %>%
  mutate(tot_bri_abun_dz=NaN,

         bri_score=NaN,
         bri_cond="BRI not Applicable",
         bri_class=NaN)

too.deep.samps<-unique(too.deep$sample_id)

bri_scores.shallow<-all.4.bri %>%
  mutate(drop_flag=case_when(sample_id%in%too.deep.samps~1,
                             sample_id%in%defaunated.samps~1,
                             TRUE~0)) %>%
  filter(drop_flag==0) %>%
  drop_na(shallow) %>% #removing taxa without a shallow pcode
  mutate(cube_abun=(abundance)^(1/3)) %>%
  group_by(sample_id, stationid, sampledate, replicate) %>% #calculating scores by sample_id (station, date, and replicate)
  mutate(tot_bri_abun=sum(abundance), #summing abundance of all taxa w/ a pcode
         tot_cube_abun=sum(cube_abun), #summing abundance of cube root abundances for all all taxa w/ a pcode
         tol_score=shallow*cube_abun) %>% #multiplying cube root abundance by pcode tolerance values
  ungroup() %>%
  group_by(sample_id, stationid, depth, sampledate,replicate, tot_bri_abun, tot_cube_abun) %>%
  summarise(numerator=sum(tol_score), .groups = "drop_last") %>% #summing the cube root abundance weighted tolerance scores by station (numerator in BRI calculation)
  ungroup() %>%
  mutate(shallow_bri_score=numerator/tot_cube_abun)

bri_scores.mid<-all.4.bri %>%
  mutate(drop_flag=case_when(sample_id%in%too.deep.samps~1,
                             sample_id%in%defaunated.samps~1,
                             TRUE~0)) %>%
  filter(drop_flag==0) %>%
  drop_na(mid) %>% #removing taxa without a mid pcode
  mutate(cube_abun=(abundance)^(1/3)) %>%
  group_by(sample_id, stationid, sampledate, replicate) %>% #calculating scores by sample_id (station, date, and replicate)
  mutate(tot_bri_abun=sum(abundance), #summing abundance of all taxa w/ a pcode
         tot_cube_abun=sum(cube_abun), #summing abundance of cube root abundances for all all taxa w/ a pcode
         tol_score=mid*cube_abun) %>% #multiplying cube root abundance by pcode tolerance values
  ungroup()%>%
  group_by(sample_id, stationid, depth, sampledate,replicate, tot_bri_abun, tot_cube_abun) %>%
  summarise(numerator=sum(tol_score), .groups = "drop_last") %>% #summing the cube root abundance weighted tolerance scores by station (numerator in BRI calculation)
  ungroup() %>%
  mutate(mid_bri_score=numerator/tot_cube_abun)

bri_scores.deep<-all.4.bri %>%
  mutate(drop_flag=case_when(sample_id%in%too.deep.samps~1,
                             sample_id%in%defaunated.samps~1,
                             TRUE~0)) %>%
  filter(drop_flag==0) %>%
  drop_na(deep) %>% #removing taxa without a deep pcode
  mutate(cube_abun=(abundance)^(1/3)) %>%
  group_by(sample_id, stationid, sampledate, replicate) %>% #calculating scores by sample_id (station, date, and replicate)
  mutate(tot_bri_abun=sum(abundance), #summing abundance of all taxa w/ a pcode
         tot_cube_abun=sum(cube_abun), #summing abundance of cube root abundances for all all taxa w/ a pcode
         tol_score=deep*cube_abun) %>% #multiplying cube root abundance by pcode tolerance values
  ungroup() %>%

  group_by(sample_id, stationid, depth, sampledate,replicate, tot_bri_abun, tot_cube_abun) %>%
  summarise(numerator=sum(tol_score), .groups = "drop_last") %>% #summing the cube root abundance weighted tolerance scores by station (numerator in BRI calculation)
  ungroup() %>%
  mutate(deep_bri_score=numerator/tot_cube_abun)

bri_scores.1<-bri_scores.shallow %>%
  left_join(., bri_scores.mid, by=c("sample_id", "stationid", "sampledate", "replicate", "depth" ), suffix = c("_s","_m")) %>%
  left_join(., bri_scores.deep, by=c("sample_id", "stationid", "sampledate", "replicate", "depth" ), suffix = c("","_d")) %>%
  mutate(bri_score=case_when(depth<25~shallow_bri_score,
                             depth>=25&depth<=35~(shallow_bri_score + mid_bri_score)/2,
                             depth>35&depth<110~mid_bri_score,
                             depth>=110&depth<=130~(mid_bri_score + deep_bri_score)/2,
                             depth>130&depth<=324~deep_bri_score),
         depth_zone=case_when(depth<25~"shallow",
                             depth>=25&depth<=35~"shallow_mid",
                             depth>35&depth<110~"mid",
                             depth>=110&depth<=130~"mid_deep",
                             depth>130&depth<=324~"deep",
                             TRUE~"out_of_range"),
         tot_bri_abun_dz=case_when(depth<25~tot_bri_abun_s,
                                depth>=25&depth<=35~ceiling((tot_bri_abun_s + tot_bri_abun_m)/2),#keeping total BRI abundance reports as an integer for ease of communication
                                depth>35&depth<110~tot_bri_abun_m,
                                depth>=110&depth<=130~ceiling((tot_bri_abun_m + tot_bri_abun)/2),
                                depth>130&depth<=324~tot_bri_abun) #note, this doesn't have a "_d" suffix due to the way the join was structured
           )



bri_scores.2<-bri_scores.1 %>%
  select(sample_id, stationid, sampledate,replicate, depth, depth_zone, shallow_bri_score, mid_bri_score,
         deep_bri_score, tot_bri_abun_dz,bri_score,) %>%

  mutate(
         bri_cond=case_when(bri_score<25 ~"Reference", #assigning traditional condition classes, based on score
                            bri_score>=25&bri_score<34~"Marginal Deviation",
                            bri_score>=34&bri_score<44~"Biodiversity Loss",
                            bri_score>=44&bri_score<72~"Function Loss",
                            bri_score>=72~"Defaunation"),
         bri_class=case_when(bri_score<25 ~1, #assigning numeric condition class
                             bri_score>=25&bri_score<34~2,
                             bri_score>=34&bri_score<44~3,
                             bri_score>=44&bri_score<72~4,
                             bri_score>=72~5)) %>%
  ungroup() %>%
  bind_rows(., defaunated, too.deep)

bri_w_station_info<-bri_scores.2 %>% #attaching station information to the BRI scores
  left_join(., select(station_info.2, stationid, latitude, longitude ), by=c("stationid")) %>%
  #add in index useage notes for each sample so the user can be aware
  mutate(note=case_when(sample_id%in%c(defaunated.samps)~"No fauna in sample",
                        tot_bri_abun_dz==0~"Caution - No organisms with P-code in Sample",
                        latitude>34.45 ~"Caution - Sample outside of the geographic range of the BRI",
                        latitude<32.52 ~"Caution - Sample outside of the geographic range of the BRI",
                        depth>200 & depth<=324 ~"Caution - Sample deeper than Bight reccommendations but within BRI calibration",
                        depth>324~"Do Not Use - Sample deeper than BRI calibration",
                        TRUE~"None")) %>%
  #remove columns used in calculation but probably not relevant to user
  select(., -tot_bri_abun_dz) %>%
  relocate(depth, latitude, longitude, .after = replicate)



#Saving the final table of BRI scores to the specified output path
write.csv(bri_w_station_info, paste(output_path, "/", file_id, " final file - BRI scores by station and replicate.csv", sep=""), row.names = FALSE)

#exporting everything into a single excel workbook with the BRI Scores and all of the interim files  for review

#define sheet names for each data frame
dataset_names <- list("final - BRI Scores"=bri_w_station_info, "int_1 - taxa submitted"=taxa_to_calc, "int_2 - taxa w pcode"=with.pcode,
                      "int_3 - taxa wo pcode"=no.pcode, "int_4 - all taxa by sample"=all.4.bri)

#export each data frame to separate sheets in same Excel file
write.xlsx(dataset_names, file = paste(output_path, "/", file_id, " Offshore BRI output summary.xlsx", sep=""))

out<-list(taxa_to_calc, with.pcode, no.pcode, all.4.bri, bri_w_station_info )

bri.outputs.2<-list(final_bri_scores=out[[5]],
                    int1_taxa_sumbitted=out[[1]],
                    int2_taxa_with_pcode=out[[2]],
                    int3_taxa_wo_pcode=out[[3]],
                    int4_all_taxa_by_sample=out[[4]])

return(bri.outputs.2) #so the user can see the scores and interim files in the in R Studio environment

}





