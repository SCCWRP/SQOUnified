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
#           taxon - character string identifying the organism. naming conventions should follow SCAMIT edition 12
#           abundance - a numeric value indicating the number of individuals counted for the specified taxon from the sample
#      3. station_path - a quoted string detailing the name and location of the .csv file with station information for each infauna
#           sample. All column names should be lower case. The expectation is that the file  will have, at a minimum,a column for:
#           stationid - unique identifier for that station, preferably as character/text value
#           depth - station depth in meters
#           latitude - station latitude in decimal degrees
#           longitude - station longitude in decimal degrees (negative values for west longitudes)
#       4. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
 ###########################################################  ######################################################################## 


offshore_bri_calc_ed12<-function(file_id, infauna_path, station_path, output_path)
  {
require(tidyverse)
  
####Input files
infauna <- read.csv(infauna_path) #the user's infauna to be submitted
station_info<-read.csv(station_path) #the user's station information to be submitted

print(paste("Saving files to ", getwd(),"/", output_path, sep=""))

load("Reference Files/pcode.RData") #pcode values by depth and taxa associated with pcodes

###            Prep the data

#ensuring data are in the correct formats
station_info.2<-station_info %>% 
  mutate(stationid=as.character(stationid),
         depth=as.numeric(depth))
infauna.2<-infauna %>% 
  mutate(stationid=as.character(stationid),
         abundance=as.numeric(abundance),
         sampledate=mdy(sampledate),
         sample_id=paste(stationid, sampledate, replicate, sep="_"))

#joining taxa names, and counts to depth and geographic information by the station field
taxa_to_calc<-infauna.2 %>% select(sample_id, stationid, replicate,sampledate,taxon, abundance) %>% 
  left_join(., select(station_info.2,stationid, depth, latitude, longitude), by="stationid") %>% 
  arrange(sample_id, desc(abundance))

#output the joined taxa-station info for review
write.csv(taxa_to_calc, paste(output_path, "/", file_id, " interim file 1 - taxa to be analyzed.csv", sep=""), row.names = FALSE) 

#rearranging the pcode values to make them easier to join to the taxa
pcodes.2<-pcodes %>% 
  pivot_longer(cols=-p_code, names_to="bri_dz", values_to="tol_val") %>%
  drop_na(tol_val)

#joining pcode values to the taxa based upon the depth zone (bri_dz)
all.4.bri<-taxa_to_calc %>% mutate(bri_dz=case_when(depth<25~"shallow", depth>=25&depth<=35~"shallow_mid",depth>35&depth<110~"mid",
                                                    depth>=110&depth<=130~"mid_deep",depth>130&depth<=324~"deep",
                                                    TRUE~"out_of_range")) %>% 
  #filter(abundance>0) %>% 
  left_join(.,ptaxa,by=c("taxon"="taxon_name")) %>% 
  left_join(.,pcodes.2, by=c("p_code", "bri_dz"))

# output taxa with their p-codes for review
with.pcode<-all.4.bri %>% distinct(taxon, bri_dz,p_code) %>% drop_na(p_code) %>% arrange(taxon)
write.csv(with.pcode, paste(output_path, "/", file_id, " interim file 2 - taxa w assigned pcodes.csv", sep=""), row.names = FALSE) 



# output taxa without p-codes for review
no.pcode<-all.4.bri %>% distinct(taxon, bri_dz,p_code) %>% filter(is.na(p_code)) %>% arrange(taxon)
write.csv(no.pcode, paste(output_path, "/", file_id, " interim file 3 - taxa w-o pcodes.csv", sep=""), row.names = FALSE) 

write.csv(all.4.bri, paste(output_path, "/", file_id, " interim file 4 - taxa and pcodes by sample.csv",sep=""), row.names=FALSE)

####          Calculate Scores

#pulling out defaunated samples and classifying them as Defaunation with a score of 5
defaunated<-all.4.bri %>% 
  filter(taxon=="NoOrganismsPresent") %>% 
  select(sample_id, stationid, sampledate, replicate) %>% 
  mutate(tot_bri_abun=0,
         tot_cube_abun=0,
         numerator=NaN,
         bri_score=NaN,
         bri_cond="Defaunation",
         bri_class=5)

defaunated.samps<-unique(defaunated$sample_id)


#pulling out the samples outside the range of the index. They get no scores nor condition classification
too.deep<-all.4.bri %>% 
  filter(bri_dz=="out_of_range") %>% 
  distinct(sample_id, stationid, sampledate, replicate) %>% 
  mutate(tot_bri_abun=NaN,
         tot_cube_abun=NaN,
         numerator=NaN,
         bri_score=NaN,
         bri_cond="BRI not Applicable",
         bri_class=NaN)

bri_scores<-all.4.bri %>%  drop_na(tol_val) %>% #drop taxa w/o a pcode
  mutate(cube_abun=(abundance)^(1/3), ) %>% #calculate cube root abundance
  group_by(sample_id, stationid, sampledate, replicate) %>% #calculating scores by sample_id (station, date, and replicate)
  mutate(tot_bri_abun=sum(abundance), #summing abundance of all taxa w/ a pcode
         tot_cube_abun=sum(cube_abun), #summing abundance of cube root abundances for all all taxa w/ a pcode
         tol_score=tol_val*cube_abun) %>% #multiplying cube root abundance by pcode tolerance values
  ungroup()

bri_scores.2<-bri_scores %>% 
  group_by(sample_id, stationid, sampledate,replicate, tot_bri_abun, tot_cube_abun) %>%
  summarise(numerator=sum(tol_score)) %>% #summing the cube root abundance weighted tolerance scores by station (numerator in BRI calculation)
  ungroup() %>%
  mutate(bri_score=numerator/tot_cube_abun, #calculation of BRI score
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

bri_station_info<-bri_scores.2 %>% #attaching station information to the BRI scores
  left_join(., station_info.2, by=c("stationid")) %>% 
  #add in index useage notes for each sample so the user can be aware
  mutate(note=case_when(sample_id%in%c(defaunated.samps)~"No fauna in sample",
                        tot_bri_abun==0~"Caution - No organisms with P-code in Sample",
                        depth>200 & depth<=324 ~"Caution - Sample deeper than Bight reccommendations but within BRI calibration",
                        depth>324~"Do Not Use - Sample deeper than BRI calibration",
                        TRUE~"None")) %>% 
  #remove cloumns used in calculation but probably not relevant to user
  select(., -tot_bri_abun, -tot_cube_abun, -numerator) %>% 
  relocate(depth, latitude, longitude, .after = replicate)



#Saving the final table of BRI scores to the specified output path
write.csv(bri_station_info, paste(output_path, "/", file_id, " final file - BRI scores by station and replicate.csv", sep=""), row.names = FALSE)

return(bri_station_info) #so you can see the scores in R Studio environment
}





