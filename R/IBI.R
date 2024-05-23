#' Compute the Index of Biotic Integrity (IBI) and IBI condition category.
#'
#' @description
#'   The IBI compares the values of four different metrics to the ranges expected under reference conditions. Each metric
#'   that is outside of the reference range increases the IBI score by one. Therefore, if all four metrics were inside
#'   the reference range, the score would be 0. Conversely, if all four metrics were outside the reference range, the
#'   value would be 4.
#'
#' @details
#'   The IBI compares the values of four different metrics to the ranges expected under reference conditions. Each metric
#'   that is outside of the reference range increases the IBI score by one. Therefore, if all four metrics were inside
#'   the reference range, the score would be 0. Conversely, if all four metrics were outside the reference range, the
#'   value would be 4.
#'
#'   The data needed to calculate the IBI are:
#'   (1) the total number of taxa,
#'   (2) the total number of mollusc taxa,
#'   (3) the abundance of \emph{Notomastus} sp., and
#'   (4) the number of sensitive taxa.
#'
#'   The total numnber of taxa, number of mollusc taxa, and abundance of \emph{Notomastus} sp. can be obtained directly
#'   from the data. The list of sensitive species should be based on the species list for Southern California Marine Bays
#'   and the percentage of sensitive taxa present is calulated as:
#'
#'   \deqn{\% \textrm{sensitive taxa} = (\textrm{number of sensistive taxa} / \textrm{total number of taxa}) \times 100}
#'
#'   The value for each metric is then compared to a reference range for that metric (Table 2).
#'   The IBI score is set to zero before comparison to the reference range. For each metric that is out of the reference
#'   range (above or below), the IBI score goes up by one.
#'
#'   <Include Table 2>
#'
#'   The IBI score is then compared to condition category thresholds (Table 3) in order to determine the IBI category and
#'   score.
#'
#'   <Include Table 3>
#'
#'
#' @param BenthicData a data frame with AT LEAST the following information with these headings:
#'
#'    \code{StationID} - an alpha-numeric identifier of the location;
#'
#'    \code{Replicate} - a numeric identifying the replicate number of samples taken at the location;
#'
#'    \code{SampleDate} - the date of sample collection;
#'
#'    \code{Latitude} - latitude in decimal degrees;
#'
#'    \code{Longitude} - longitude in decimal degrees. Make sure there is a negative sign for the Western coordinates;
#'
#'    \code{Species} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#'    \code{Abundance} - the number of each Species observed in a sample;
#'
#'    \code{Salinity} - the salinity observed at the location in PSU, ideally at time of sampling.
#'
#' @usage
#' IBI(benthic_data)
#'
#' @examples
#' data(benthic_sampledata)
#' IBI(BenthicData)
#'
#' @import dplyr



##########################################################################################################################
#
##########################################################################################################################
#' @export
IBI <- function(BenthicData, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ), verbose = T)
{
  # load("data/SoCal_SQO_Infauna_LU_updated_4.7.20.RData")

  # Prepare the given data frame so that we can compute the IBI score and categories
  ibi_data <- BenthicData %>%
    filter(Exclude!="Yes") %>%
    left_join(sqo.list.new, by = c('Taxon'='TaxonName')) %>%
    mutate_if(is.numeric, list(~na_if(., -88))) %>%
    select('StationID','SampleDate', 'Replicate','Taxon','Abundance','Stratum', 'Phylum', 'IBISensitive', "Mollusc", "Crustacean") %>%
    #rename(B13_Stratum = Stratum) %>%
    mutate(n=if_else(Taxon=="NoOrganismsPresent", 0,1))

  ### SQO IBI - 1
  # columns needed in RBI: B13_Stratum, StationID, Replicate, Phylum, NumofTaxa
  ibi1 <- ibi_data %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NumOfTaxa =sum(n))


  ### SQO IBI - 2
  ibi2 <- ibi_data %>%
    filter(Mollusc == "Mollusc") %>%
    group_by(Stratum, StationID, SampleDate,Replicate) %>%
    summarise(NumOfMolluscTaxa = length(Taxon))




  ### SQO RBI - 3 - 2
  ibi3_2 <- ibi_data %>%
    filter(str_detect(Taxon,"Notomastus")) %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NotomastusAbun = sum(Abundance))



  ### SQO IBI - 4 - 2
  ibi4_2 <- ibi_data %>%
    filter(IBISensitive=="S") %>%
    left_join(ibi1, by=c("Stratum", "StationID", "SampleDate", "Replicate")) %>%
    group_by(Stratum, StationID, SampleDate, Replicate, NumOfTaxa) %>%
    summarise(SensTaxa = length(Taxon)) %>%
    mutate(PctSensTaxa=(SensTaxa/NumOfTaxa)*100) %>%
    select(Stratum, StationID, SampleDate, Replicate, PctSensTaxa)



  ### Reference ranges for IBI metrics in Southern California Marine Bays
  ### [ Table 5.4 (p. 77, Technical Manual, 2014) ]
  ibi_ref_ranges_table <- data.frame(ref_low = c(13, 2, 0, 19),
                                     ref_high = c(99, 25, 59, 47.1))
  row.names(ibi_ref_ranges_table) <- c("NumOfTaxa", "NumOfMolluscTaxa", "NotomastusAbun", "PctSensTaxa")


  ### IBI category response ranges for Southern California Marine Bays
  ### [ Table 5.5 (p. 77, Technical Manual, 2014) ]
  ibi_category_response_table <- data.frame(ibi_score = as.factor(c(0, 1, 2, 3, 4)),
                                            category = as.factor(c("Reference",
                                                                   "Low Disturbance",
                                                                   "Moderate Disturbance",
                                                                   "High Disturbance",
                                                                   "High Disturbance")),
                                            category_score = as.factor(c(1, 2, 3, 4, 4)))

  ### B13 IBI Metrics:
  # We stitch together all the necessary IBI metrics to determine the IBI index.
  # Each of the metrics is then compared to the tables listed above (Table 5.4 and Table 5.5) to determine the IBI score,
  # the IBI Category, and IBI Category Score
  ibi_metrics <- ibi1 %>%
    full_join(ibi2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    full_join(ibi3_2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    full_join(ibi4_2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    replace(.,is.na(.),0) %>%
    #select(B13_Stratum, StationID, SampleDate, Replicate, NumOfTaxa, NumOfMolluscTaxa, NotomastusAbun, PctSensTaxa) %>%
    # We replace any NAs with 0 so that we can compare the values to the tables listed above
    # The IBI score is set to zero before comparison the reference range.
    mutate(Score = 0) %>%
    # For each metric that is out of the reference range (above or below), the IBI score goes up by one.
    mutate(Score = if_else((NumOfTaxa < ibi_ref_ranges_table["NumOfTaxa",]$ref_low  | NumOfTaxa > ibi_ref_ranges_table["NumOfTaxa",]$ref_high),
                                      Score + 1, Score)) %>%
    mutate(Score = if_else((NumOfMolluscTaxa < ibi_ref_ranges_table["NumOfMolluscTaxa",]$ref_low  | NumOfMolluscTaxa > ibi_ref_ranges_table["NumOfMolluscTaxa",]$ref_high),
                                      Score + 1, Score)) %>%
    mutate(Score = if_else((NotomastusAbun < ibi_ref_ranges_table["NotomastusAbun",]$ref_low  | NotomastusAbun > ibi_ref_ranges_table["NotomastusAbun",]$ref_high),
                                      Score + 1, Score)) %>%
    mutate(Score = if_else((PctSensTaxa < ibi_ref_ranges_table["PctSensTaxa",]$ref_low  | PctSensTaxa > ibi_ref_ranges_table["PctSensTaxa",]$ref_high),
                                      Score + 1, Score)) %>%
    # The IBI score is then compared to condition category response ranges (Table 5.5) to determine the IBI category and category score.
    mutate(Category = case_when(Score == 0 ~ "Reference", Score == 1 ~ "Low Disturbance", Score == 2 ~ "Moderate Disturbance", (Score == 3 | Score == 4) ~ "High Disturbance")) %>%
    mutate(`Category Score` = case_when(Score == 0 ~ 1, Score == 1 ~ 2, Score == 2 ~ 3, (Score == 3 | Score == 4) ~ 4)) %>%
    mutate(Index = "IBI") %>%
    distinct()

    return(ibi_metrics)

  #write.csv(ibi_metrics, file = "data/ibi-metrics.csv", row.names = FALSE)

}
