#Function designed to input benthic data and station information for offshore BRI calculation
#producing both interim and final output files

#This function requires that you have the following packages installed on your machine:
#       1. tidyverse
#       
#
# in your r environments you can install packages with:  install.packages(".......")

# The BRI is an abundance-weighted pollution tolerance index that assesses the condition of a given sample
# based upon the tolerance/sensitivity to pollution of the fauna observed in a sample. Specifics on calibration and validateion
# of the index can be found in Smith et al. 2001. Benthic Response Index for Assessing Infaunal Communities on the Southern California 
# Mainland Shelf. Ecological applications 11:1073-1087


####################################################################_#############################################################
# Instructions for use
# The function has 4 required inputs:
#      1. file_id - this is a quoted string for you to identify the data the BRI scores are associated with
#           e.g., "Bight 23" or "2024 southern shelf" - this will be used to name all of the output files
#      2. infauna_path - a quoted string detailing the name and location of the .csv file with infauna abundance
#         data. Remember with path names, that R use forward slashes "/" not the normal backslash "\" that windows typically uses.
#         All column names should be lower case. The expectation is that the file will have, at a minimum, a column for each of:
#           station_id - unique identifier for that station, preferably formatted as character/text value
#           sample_date - the date on which the sample was collected, must be in a mm/dd/yyyy format
#           replicate - a number identifying the replicate infauna sample collected from the specified station on the specified date
#           taxon - character string identifying the organism. naming conventions should follow SCAMIT edition 12
#           abundance - a numeric value indicating the number of individuals counted for the specified taxon from the sample
#      3. station_path - a quoted string detailing the name and location of the .csv file with station information for each infauna
#           sample. All column names should be lower case. The expectation is that the file  will have, at a minimum,a column for:
#           station_id - unique identifier for that station, preferably as character/text value
#           depth - station depth in meters
#           latitude - station latitude in decimal degrees
#           longitude - station longitude in decimal degrees (negative values for west longitudes)
#       4. output_path - a quoted string detailing the location where you want the output files to be saved. Remember to use "/" not "\"
 ###########################################################  ######################################################################## 


offshore_bri.generic<-function(file_id, infauna_path, station_path, output_path)
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
  mutate(station_id=as.character(station_id),
         depth=as.numeric(depth))
infauna.2<-infauna %>% 
  mutate(station_id=as.character(station_id),
         abundance=as.numeric(abundance),
         sample_date=ymd(sample_date),
         sample_id=paste(station_id, sample_date, replicate, sep="_"))

#joining taxa names, and counts to depth and geographic information by the station field
taxa_to_calc<-infauna.2 %>% select(sample_id, station_id, replicate,sample_date,taxon, abundance) %>% 
  left_join(., select(station_info.2,station_id, depth, latitude, longitude), by="station_id") %>% 
  arrange(sample_id, desc(abundance))

#output the joined taxa-station info for review
write.csv(taxa_to_calc, paste(output_path, "/", file_id, " interim file 1 - taxa to be analyzed.csv", sep=""), row.names = FALSE) 

#rearranging the pcode values to make them easier to join to the taxa
pcodes.2<-pcodes %>% 
  pivot_longer(cols=-p_code, names_to="bri_dz", values_to="tol_val") %>%
  drop_na(tol_val)

#joining pcode values to the taxa based upon the depth zone (bri_dz)
all.4.bri<-taxa_to_calc %>% mutate(bri_dz=case_when(depth<25~"shallow", depth>=25&depth<=35~"shallow_mid",depth>35&depth<110~"mid",
                                                    depth>=110&depth<=130~"mid_deep",depth>130&depth<=324~"deep")) %>% filter(abundance>0) %>% 
  left_join(.,ptaxa,by=c("taxon"="taxon_name")) %>% 
  left_join(.,pcodes.2, by=c("p_code", "bri_dz"))

# output taxa with their p-codes for review
with.pcode<-all.4.bri %>% distinct(taxon, bri_dz,p_code) %>% drop_na(p_code) %>% arrange(taxon)
write.csv(with.pcode, paste(output_path, "/", file_id, " interim file 2 - taxa w assigned pcodes.csv", sep=""), row.names = FALSE) 



# output taxa without p-codes for review
no.pcode<-all.4.bri %>% distinct(taxon, bri_dz,p_code) %>% filter(is.na(p_code)) %>% arrange(taxon)
write.csv(no.pcode, paste(output_path, "/", file_id, " interim file 3 - taxa w-o pcodes.csv", sep=""), row.names = FALSE) 

write.csv(all.4.bri, paste(output_path, "/", file_id, " interim file 4 - taxa and pcodes by sample.csv",sep=""), row.names=FALSE)

# provide some guidance on index application based on water depth

samp.depths<-station_info.2 %>% 
  mutate(note=case_when(depth<=200~"none", #Bight practice is to only use the BRI to 200m depth
                        depth>200 & depth <=324 ~ "Index calibrated to 324m, but Bight practice is only to apply to 200m",
                        #The original Smith etl 2001 paper calibrated the BRI to 324m depth
                        depth>324 ~ "BRI not reccommended at sample depth"))
                        #There is some evidence (Gillett, D.J., Gilbane, L., and Schiff, K.C. 2019. Benthic Infauna of the Southern 
                          #California Bight Continental Slope: Characterizing Community Structure for the Development of an Index of Disturbance. 
                          #Camarillo (CA): US Department of the Interior, Bureau of Ocean Energy Management. OCS Study BOEM 2019-050. 157 p.) that
                          #the BRI could be applied to 400m depth but that has not been fully vetted

# output depth range applicability information for the user to review
write.csv(samp.depths, paste(output_path, "/", file_id, " interim file 5 - BRI applicability for each sample.csv", sep=""), row.names=FALSE)



####          Calculate Scores

bri_scores<-all.4.bri %>%  drop_na(tol_val) %>% #drop taxa w/o a pcode
  mutate(cube_abun=(abundance)^(1/3), ) %>% #calculate cube root abundance
  group_by(sample_id, station_id, sample_date, replicate) %>% #calculating scores by sample_id (station, date, and replicate)
  mutate(tot_bri_abun=sum(abundance), #summing abundance of all taxa w/ a pcode
         tot_cube_abun=sum(cube_abun), #summing abundance of cube root abundances for all all taxa w/ a pcode
         tol_score=tol_val*cube_abun) %>% #multiplying cube root abundance by pcode tolerance values
  ungroup()

bri_scores.2<-bri_scores %>% 
  group_by(sample_id, station_id, sample_date,replicate, tot_bri_abun, tot_cube_abun) %>%
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
  ungroup()

bri_station_info<-bri_scores.2 %>% #attaching station information to the BRI scores
  left_join(., samp.depths, by=c("station_id"))



#Saving the final table of BRI scores to the specified output path
write.csv(bri_station_info, paste(output_path, "/", file_id, " final file - BRI scores by station and replicate.csv", sep=""), row.names = FALSE)

return(bri_station_info) #so you can see the scores in R Studio environment
}





