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
#' BRI_(benthic_data)
#'
#' @examples
#' data(benthic_sampledata) # load sample data
#' BRI_(benthic_sampledata) # see the output
#'
#' @import vegan
#' @import reshape2
#' @importFrom dplyr left_join filter rename select mutate group_by summarize summarise case_when
#' @importFrom lubridate ymd
#'
#' @export

BRI_ <- function(BenthicData) #BenthicData will need to be the species abundances for each sample in the correct format noted above and in support material
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


