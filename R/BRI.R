# ORIGINAL BRI ------------------------------------------------------------------------------
#' Compute the benthic response index (BRI) score and BRI condition category.
#'
#' @description
#'   The BRI is the abundance weighted pollution tolerance score of the organisms present in a benthic sample. The higher
#'   the BRI score, the more degraded the benthic community represented by the sample.
#'
#' @details
#'   The BRI is the 4th root relative abundance weighted pollution tolerance score of the organisms present in a benthic sample. The higher
#'   the BRI score, the more degraded the benthic community represented by the sample.
#'
#'   Two types of data are needed to calculate the BRI:
#'
#'   (1) the abundance of each species
#'   (2) species-specific pollution tolerance score (aka, P Value)
#'
#'   Tolerance Values are stored in the Southern California SQO Species List provided with this coded. Species names are periodically
#'   updated by benthic experts.
#'
#'   The BRI is only calculated from those taxa with a tolerance score. The first step in the BRI calculation is to compute the 4th root
#'   of the abundance of each taxon in the sample that have an associated tolerance score
#'   The next step is to multiply the 4th root abundance value by the tolerance score for each taxon.
#'   The next step is to sum all of the 4th root abundance values in a given sample.
#'   The actual BRI score is calculated as:
#'
#'   \deqn{ \frac{\sum \left(\sqrt[p]{\textrm{Abundance}} \right) \times P}{\sum \sqrt[p]{\textrm{Abundance}}} }
#'
#'   The last step is to convert the BRI score to condition category using the category thresholds listed in Table 5.
#'
#'   <Table 5. To be included in R markdown file>
#'
#'
#'
#' @param BenthicData a data frame with the following headings
#'
#'    \strong{\code{StationID}} - an alpha-numeric identifier of the location;
#'
#'    \strong{\code{Replicate}} - a numeric identifying the replicate number of samples taken at the location;
#'
#'    \strong{\code{SampleDate}} - the date of sample collection;
#'
#'    \strong{\code{Latitude}} - latitude in decimal degrees;
#'
#'    \strong{\code{Longitude}} - longitude in decimal degrees.
#'    Make sure there is a negative sign for the Western coordinates;
#'
#'    \strong{\code{Species}} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#' @usage
#' BRI(benthic_data)
#'
#' @examples
#' data(benthic_sampledata) # load sample data
#' BRI(benthic_sampledata) # see the output
#'
#' @import vegan
#' @import reshape2
#' @importFrom dplyr left_join filter rename select mutate group_by summarize summarise case_when
#'
#' @export
BRI <- function(BenthicData, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = F)
{
  # This (as of now) is the main function used by the R package

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\n## BEGIN: BRI function.\n', logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'benthic.bri.input.RData'
  if (verbose) {
    save(BenthicData, file = file.path( dirname(logfile), rawinput.filename ))
  }

  # Create code block and download link to BRI input
  writelog(
    'Input to BRI - BRI-step0.csv',
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') ### This will load a dataframe called 'BenthicData' into your environment"),
    data = BenthicData,
    verbose = verbose
  )
  create_download_link(data = BenthicData, logfile = logfile, filename = 'BRI-step0.csv', linktext = 'Download BRI initial input', verbose = verbose)


  # SQO List New (Gets joined to the initial input for BRI)
  writelog(
    'SQO List New (Gets joined to the initial input for BRI)',
    logfile = logfile,
    data = sqo.list.new,
    verbose = verbose
  )
  create_download_link(data = sqo.list.new, logfile = logfile, filename = 'BRI-sqo.list.new.csv', linktext = 'Download BRI "sqo.list.new" dataframe', verbose = verbose)

  # BRI Step 1
  # Join with sqo.list.new
  bri1 <- BenthicData %>%
    left_join(sqo.list.new, by = c('Taxon' = 'TaxonName'))

  # Write to the logs for BRI Step 1
  writelog(
    '\nBRI Step 1 - Join with sqo.list.new',
    logfile = logfile,
    code = "
      bri1 <- BenthicData %>%
        left_join(sqo.list.new, by = c('Taxon' = 'TaxonName'))
    ",
    data = bri1,
    verbose = verbose
  )
  create_download_link(data = bri1, logfile = logfile, filename = 'BRI-step1.csv', linktext = 'Download BRI step 1', verbose = verbose)



  # BRI Step 2
  # Remove missing tolerance scores
  bri2 <- bri1 %>%
    # I assume that the next line is something they had in there as a method of removing duplicates
    # for this reason, this next line will likely be eliminated.
    # They grouped by all the columns that were selected (In query BRI - 1)
    # Instead, if need be we can use something from dplyr that deals with duplicates
    # I actually found that it didn't appear to make a difference
    filter(!is.na(ToleranceScore)) %>%
    #rename(Stratum) %>%
    select(Stratum, StationID, SampleDate, Replicate, Taxon, Abundance, ToleranceScore)

  # Write to the logs for BRI Step 2
  writelog(
    '\nBRI Step 2 - Remove missing tolerance scores',
    logfile = logfile,
    code = "
      bri2 <- bri1 %>%
        filter(!is.na(ToleranceScore)) %>%
        select(Stratum, StationID, SampleDate, Replicate, Taxon, Abundance, ToleranceScore)
    ",
    data = bri2,
    verbose = verbose
  )
  create_download_link(data = bri2, logfile = logfile, filename = 'BRI-step2.csv', linktext = 'Download BRI step 2', verbose = verbose)



  # BRI Step 3
  # Take the fourth root of the abundance
  bri3 <- bri2 %>%
    mutate(
      fourthroot_abun = Abundance ** 0.25,
      tolerance_score = fourthroot_abun * ToleranceScore
    )

  # Write to the logs for BRI Step 3
  writelog(
    '\nBRI Step 3 - Take the fourth root of the abundance',
    logfile = logfile,
    code = "
      bri3 <- bri2 %>%
        mutate(
          fourthroot_abun = Abundance ** 0.25,
          tolerance_score = fourthroot_abun * ToleranceScore
        )
    ",
    data = bri3,
    verbose = verbose
  )
  create_download_link(data = bri3, logfile = logfile, filename = 'BRI-step3.csv', linktext = 'Download BRI step 3', verbose = verbose)



  writelog('\nNext get the Score - group by Stratum, StationID, SampleDate, Replicate and do: (sum of the tolerance scores)/(sum of fourthroot abundances)', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('\nThen Get Categories (CASQO Technical Manual 3rd Edition Page 72 - Table 4.24)', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- < 39.96 is Reference', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- >=39.96 and <49.15 is Low', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- >=49.15 and <73.27 is Reference', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- >=73.27 is High', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)




  # BRI Step 4
  # Sum of tolerance scores divided by fourthroot of abundance
  bri4 <- bri3 %>%
    group_by(
      Stratum, StationID, SampleDate, Replicate
    ) %>%
    summarize(
      Score = sum(tolerance_score, na.rm = T) / sum(fourthroot_abun, na.rm = T)
    )

  # Write to the logs for BRI Step 4
  writelog(
    '\nBRI Step 4 - Sum of tolerance scores divided by fourthroot of abundance',
    logfile = logfile,
    code = "
      bri4 <- bri3 %>%
        group_by(
          Stratum, StationID, SampleDate, Replicate
        ) %>%
        summarize(
          Score = sum(tolerance_score, na.rm = TRUE) / sum(fourthroot_abun, na.rm = TRUE)
        )
    ",
    data = bri4,
    verbose = verbose
  )
  create_download_link(data = bri4, logfile = logfile, filename = 'BRI-step4.csv', linktext = 'Download BRI step 4', verbose = verbose)




  # BRI Step 5
  # Output the BRI category given the BRI score and the thresholds for Southern California Marine Bays - [CASQO Technical Manual 3rd Edition Page 72 - Table 4.24]
  bri5 <- bri4 %>%
    mutate(
      Category = case_when( (Score < 39.96) ~ "Reference",
                            (Score >= 39.96 & Score < 49.15) ~ "Low Disturbance",
                            (Score >= 49.15 & Score < 73.27) ~ "Moderate Disturbance",
                            (Score >= 73.27) ~ "High Disturbance"
      ))

  # Write to the logs for BRI Step 5
  writelog(
    '\nBRI Step 5 - Output the BRI category given the BRI score and the thresholds for Southern California Marine Bays - [CASQO Technical Manual 3rd Edition Page 72 - Table 4.24]',
    logfile = logfile,
    code = "
      bri5 <- bri4 %>%
        mutate(
          Category = case_when(
            Score < 39.96 ~ 'Reference',
            Score >= 39.96 & Score < 49.15 ~ 'Low Disturbance',
            Score >= 49.15 & Score < 73.27 ~ 'Moderate Disturbance',
            Score >= 73.27 ~ 'High Disturbance'
          )
        )
    ",
    data = bri5,
    verbose = verbose
  )
  create_download_link(data = bri5, logfile = logfile, filename = 'BRI-step5.csv', linktext = 'Download BRI step 5', verbose = verbose)




  # BRI Final Step
  # Output the BRI category score given the category for thresholds for Southern CA Marine Bays
  bri_final <- bri5 %>%
    mutate(
      `Category Score` = case_when( (Category == "Reference") ~ 1,
                                    (Category == "Low Disturbance") ~ 2,
                                    (Category == "Moderate Disturbance") ~ 3,
                                    (Category == "High Disturbance") ~ 4 )
    ) %>%
    dplyr::mutate(Index = "BRI")

  # Write to the logs for BRI Final Step
  writelog(
    '\nBRI Final Step - Output the BRI category score given the category for thresholds for Southern CA Marine Bays',
    logfile = logfile,
    code = "
      bri_final <- bri5 %>%
        mutate(
          `Category Score` = case_when(
            Category == 'Reference' ~ 1,
            Category == 'Low Disturbance' ~ 2,
            Category == 'Moderate Disturbance' ~ 3,
            Category == 'High Disturbance' ~ 4
          )
        ) %>%
        dplyr::mutate(Index = 'BRI')
    ",
    data = bri_final,
    verbose = verbose
  )
  create_download_link(data = bri_final, logfile = logfile, filename = 'BRI-final.csv', linktext = 'Download BRI Final Step', verbose = verbose)



  writelog('\n## END: BRI function.\n', logfile = logfile, verbose = verbose)

  return(bri_final)
}







# ------------------------------------------------------------------ GENERIC OFFSHORE BRI --------------------------------------------------------------------------

# Nick Haring and David Gillett both pushed their own versions of the offshore BRI function - I am putting those functions in here
# tampering with the above BRI function may break certain things in terms of adjusting the package for the SQO audit of our calculations used for the '18 synth report
# So until that is overwith, I would like to be the one making any changes to the BRI function,
# (Although TBH it should be perfectly fine as long as the input args remain the same and the output dataframe has the same column name/datatype structure)
# But the below functions may be modified in any way without reservation
# So long as the have the roxygen comments typed correctly, they will be available when the package is loaded
#  - Robert (6/27/2024)

# eventually, I believe we will want to replace the above "original" BRI function with the ones below


# Nick's Generic BRI ------------------------------------------------------------------------------------------------------------
#' Compute the benthic response index (BRI) score and BRI condition category. This version is the generic offshore one
#'
#' @description
#'   The BRI is the abundance weighted pollution tolerance score of the organisms present in a benthic sample. The higher
#'   the BRI score, the more degraded the benthic community represented by the sample.
#'   This function should also work with offshore data
#'
#' @import vegan
#' @import reshape2
#' @importFrom dplyr left_join filter rename select mutate group_by summarize summarise case_when
#' @importFrom lubridate ymd
#'
#' @export

BRI.GenericOffshore.NH <- function(BenthicData) #BenthicData will need to be the species abundances for each sample in the correct format noted above and in support material
{


  # I put these in the import statements in the roxygen comments
  # I think the only function that may have came from tidyverse that this function needed was ymd from lubridate
  # If others come up we can add them as needed, or if it becomes too much we can import the whole tidyverse package (@import tidyverse)
  # - Robert 06.10.2024

  #loading in packages needed to run function
  # require(tidyverse)

  # It appears that this version of the function works with all lowercase column names - Robert 06.10.2024
  names(BenthicData) <- names(BenthicData) %>% tolower()

  #loading in SQO species list that contains p codes, amongst other things

  load("data/SoCal SQO LU 4_7_20.RData")
  #I've created an issue, but we will need to periodically update the support files for the different indices,
  # e.g., the SQO look up list or BRI ptaxa list.
  # Do we want to date stamp the names of the dataframes and the RData files as they are updated? e.g., sqo.list.4_7_20 vs. a more generic name like sqo.list.
  # if the former, we will need to update the internal call of the index functions to make sure it is pulling the correct version. However it is
  # explicit as to what version is being used. Alternatively: if we do not date stamp the files, then the code would not need to be upadated each time the
  # support file is updated. This requires less maintenance. However, it then becomes less clear which exact version of the the support file is being used.

  #create empty dataframe to populate w/ bri scores
  bri.out.null<-tibble(stationid="dummy",
                       replicate=NaN,
                       sampledate=ymd("2000-01-1"),
                       index="BRI",
                       score=NaN,
                       condition.category=NA,
                       condition.category.score=NA,
                       note=NA)

  #incase a sample had no animals (e.g., taxon=NoOrganismsPresent), we force it into the High Disturbance category.
  #the calculator would not be able to process that sample and would drop it, so we deal with it apriori
  defaunated<-BenthicData %>%
    filter(taxon=="NoOrganismsPresent") %>%
    mutate(index="BRI",score=NaN, condition.category="High Disturbance", condition.category.score=4, note="Defaunated Sample") %>%
    select(-taxon, -salinity,-exclude, -abundance, -latitude, -longitude )

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
  #export as an interim file that the user can review
  return(taxa_w_pvalue)
  write.csv(taxa_w_pvalue, paste(output.path,"/", file.name, " interim - taxa with a tolerance score.csv", sep=""), row.names = FALSE)

  #identify those taxa in the submitted data without a tolerance score and how many samples they occur in
  taxa_wo_pvalue<-all.for.bri%>%
    group_by(taxon, ToleranceScore) %>%
    summarise(Freq_of_Occ=length(stationid)) %>%
    ungroup() %>%
    filter(is.na(ToleranceScore)) %>%
    select(-ToleranceScore)

  #export as an interim file that the user can review
  return(taxa_wo_pvalue)
  write.csv(taxa_wo_pvalue, paste(output.path, "/", file.name, " interim taxa without a tolerance score.csv", sep=""), row.names = FALSE)

  #calculate some summary values for context of index utility evaluation
  # tol.inventory<-all.for.bri %>%
  #   mutate(tol.flag=if_else(is.na(ToleranceScore), "wo_tol_value", "w_tol_value")) %>% #group taxa by those without and with a tolerance score
  #   group_by(stationid, sampledate, replicate) %>%
  #   mutate(tot_abun=sum(abundance), S=length(taxon)) %>% #calcualte total abundance and taxa richness for each sample
  #   ungroup() %>%
  #   group_by(stationid, sampledate, replicate, tot_abun, S,tol.flag) %>%
  #   summarise(tol.abun=sum(abundance), tol.s=length(taxon)) %>% #calculating percent of abundance or richness with a tolerance value per sample
  #   ungroup() %>%
  #   mutate(pct_abun=round((tol.abun/tot_abun)*100, digits = 1), pct_taxa=round((tol.s/S)*100, digits=1)) %>%
  #   pivot_wider(id_cols =c(stationid, sampledate, replicate), names_from = tol.flag, values_from = c(pct_abun, pct_taxa) )#manipulating the shape of the data



  bri.out<-all.for.bri %>%
    drop_na(ToleranceScore) %>%
    mutate(fourthroot_abun = abundance ** 0.25,
           tolerance_value = fourthroot_abun * ToleranceScore) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarize(numerator = sum(tolerance_value, na.rm = T), denomenator= sum(fourthroot_abun, na.rm = T), score=numerator/denomenator) %>%
    select(stationid, sampledate, replicate, score) %>%
    # Output the BRI category given the BRI score and the thresholds for Southern California Marine Bays
    mutate(
      condition.category = case_when( (score < 39.96) ~ "Reference",
                                      (score >= 39.96 & score < 49.15) ~ "Low Disturbance",
                                      (score >= 49.15 & score < 73.27) ~ "Moderate Disturbance",
                                      (score >= 73.27) ~ "High Disturbance"
      )) %>%
    # Output the BRI category score given the category for thresholds for Southern CA Marine Bays
    mutate(
      condition.category.score = case_when( (condition.category == "Reference") ~ 1,
                                            (condition.category == "Low Disturbance") ~ 2,
                                            (condition.category == "Moderate Disturbance") ~ 3,
                                            (condition.category == "High Disturbance") ~ 4 ),
      index="BRI")

  bri.out.2<-bri.out.null %>%
    bind_rows(bri.out, defaunated)

  if(length(bri.out.2$stationid)>1)#if BRI scores are calculated, functin will drop the dummy data placeholders and report the data for the submitted samples
  {
    bri.out.3<-bri.out.2 %>%
      filter(stationid!="dummy") %>%
      mutate(note=if_else(is.na(note), "none", note))
  }

  else
  {
    bri.out.3<-bri.out.2 %>% # if BRI scores are not calculated, the functin will only report the dummy data placeholders
      mutate(note="BRI scores not caculated")
  }
  return(bri.out.3)
}


# David's Generic BRI ------------------------------------------------------------------------------------------------------------
#################################################################################################################################-
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
###################################################################################################################################-


offshore_bri_calc_ed12<-function(file_id, infauna_path, station_path, output_path)
{
  require(tidyverse)

  ####Input files
  infauna <- read.csv(infauna_path) #the user's infauna to be submitted
  station_info<-read.csv(station_path) #the user's station information to be submitted

  load("data/pcode.RData") #pcode values by depth and taxa associated with pcodes

  ###            Prep the data

  #ensuring data are in the correct formats
  station_info.2<-station_info %>%
    mutate(station_id=as.character(station_id),
           depth=as.numeric(depth))
  infauna.2<-infauna %>%
    mutate(station_id=as.character(station_id),
           abundance=as.numeric(abundance),
           sample_date=mdy(sample_date),
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
    left_join(.,ptaxa,by=c("taxon"="TaxonName")) %>% left_join(.,pcodes.2, by=c("PCode"="p_code", "bri_dz"))

  # output taxa with their p-codes for review
  with.pcode<-all.4.bri %>% distinct(taxon, bri_dz,PCode) %>% drop_na(PCode) %>% arrange(taxon)
  write.csv(with.pcode, paste(output_path, "/", file_id, " interim file 2 - taxa w assigned pcodes.csv", sep=""), row.names = FALSE)



  # output taxa without p-codes for review
  no.pcode<-all.4.bri %>% distinct(taxon, bri_dz,PCode) %>% filter(is.na(PCode)) %>% arrange(taxon)
  write.csv(no.pcode, paste(output_path, "/", file_id, " interim file 3 - taxa w-o pcodes.csv", sep=""), row.names = FALSE)

  write.csv(all.4.bri, paste(output_path, "/", file_id, " interim file 4 - taxa and pcodes by sample.csv",sep=""), row.names=FALSE)

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
    left_join(., station_info.2, by=c("station_id"))

  return(bri_station_info) #so you can see the scores in R Studio environment

  #Saving the final table of BRI scores to the specified output path
  write.csv(bri_station_info, paste(output_path, "/", file_id, " final file - BRI scores by station and replicate.csv", sep=""), row.names = FALSE)


}








