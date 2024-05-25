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
IBI <- function(BenthicData, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt'), verbose = TRUE) {

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\nBEGIN: IBI function.\n', logfile = logfile, verbose = verbose)

  # Log initial input data
  writelog('*** DATA *** Input to IBI function - IBI-step0.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(BenthicData, logfile = file.path(dirname(logfile), 'IBI-step0.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  # Prepare the given data frame
  ibi_data <- BenthicData %>%
    filter(Exclude != "Yes") %>%
    left_join(sqo.list.new, by = c('Taxon' = 'TaxonName')) %>%
    mutate_if(is.numeric, list(~na_if(., -88))) %>%
    select('StationID', 'SampleDate', 'Replicate', 'Taxon', 'Abundance', 'Stratum', 'Phylum', 'IBISensitive', "Mollusc", "Crustacean") %>%
    mutate(n = if_else(Taxon == "NoOrganismsPresent", 0, 1))

  writelog('*** DATA *** Prepared IBI data - IBI-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi_data, logfile = file.path(dirname(logfile), 'IBI-step1.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### SQO IBI - 2
  ibi1 <- ibi_data %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NumOfTaxa = sum(n))

  writelog('number of taxa', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('*** DATA *** IBI Step 2 - IBI-step2.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi1, logfile = file.path(dirname(logfile), 'IBI-step2.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### SQO IBI - 3
  ibi2 <- ibi_data %>%
    filter(Mollusc == "Mollusc") %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NumOfMolluscTaxa = length(Taxon))

  writelog('number of molluscs', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('*** DATA *** IBI Step 3 - IBI-step3.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi2, logfile = file.path(dirname(logfile), 'IBI-step3.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### SQO IBI - 4
  ibi3_2 <- ibi_data %>%
    filter(str_detect(Taxon, "Notomastus")) %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NotomastusAbun = sum(Abundance))

  writelog('sum abundance of Notomastus', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('*** DATA *** IBI Step 4 - IBI-step.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi3_2, logfile = file.path(dirname(logfile), 'IBI-step4.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### SQO IBI - 5
  ibi4_2 <- ibi_data %>%
    filter(IBISensitive == "S") %>%
    left_join(ibi1, by = c("Stratum", "StationID", "SampleDate", "Replicate")) %>%
    group_by(Stratum, StationID, SampleDate, Replicate, NumOfTaxa) %>%
    summarise(SensTaxa = length(Taxon)) %>%
    mutate(PctSensTaxa = (SensTaxa / NumOfTaxa) * 100) %>%
    select(Stratum, StationID, SampleDate, Replicate, PctSensTaxa)

  writelog('Percentage of Sensitive taxa', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('*** DATA *** IBI Step 5 - IBI-step5.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi4_2, logfile = file.path(dirname(logfile), 'IBI-step5.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### Reference ranges for IBI metrics in Southern California Marine Bays
  writelog('*** REFERENCE RANGES *** IBI Reference Ranges - ibi_ref_ranges_table.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi_ref_ranges_table, logfile = file.path(dirname(logfile), 'ibi_ref_ranges_table.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### IBI category response ranges for Southern California Marine Bays
  writelog('*** CATEGORY RESPONSE RANGES *** IBI Category Response Ranges - ibi_category_response_table.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi_category_response_table, logfile = file.path(dirname(logfile), 'ibi_category_response_table.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  ### IBI Metrics:
  ibi_metrics <- ibi1 %>%
    full_join(ibi2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    full_join(ibi3_2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    full_join(ibi4_2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    replace(., is.na(.), 0) %>%

    mutate(Score = 0) %>%
    mutate(Score = if_else((NumOfTaxa < ibi_ref_ranges_table["NumOfTaxa",]$ref_low | NumOfTaxa > ibi_ref_ranges_table["NumOfTaxa",]$ref_high), Score + 1, Score)) %>%
    mutate(Score = if_else((NumOfMolluscTaxa < ibi_ref_ranges_table["NumOfMolluscTaxa",]$ref_low | NumOfMolluscTaxa > ibi_ref_ranges_table["NumOfMolluscTaxa",]$ref_high), Score + 1, Score)) %>%
    mutate(Score = if_else((NotomastusAbun < ibi_ref_ranges_table["NotomastusAbun",]$ref_low | NotomastusAbun > ibi_ref_ranges_table["NotomastusAbun",]$ref_high), Score + 1, Score)) %>%
    mutate(Score = if_else((PctSensTaxa < ibi_ref_ranges_table["PctSensTaxa",]$ref_low | PctSensTaxa > ibi_ref_ranges_table["PctSensTaxa",]$ref_high), Score + 1, Score)) %>%
    mutate(Category = case_when(Score == 0 ~ "Reference", Score == 1 ~ "Low Disturbance", Score == 2 ~ "Moderate Disturbance", (Score == 3 | Score == 4) ~ "High Disturbance")) %>%
    mutate(`Category Score` = case_when(Score == 0 ~ 1, Score == 1 ~ 2, Score == 2 ~ 3, (Score == 3 | Score == 4) ~ 4)) %>%
    mutate(Index = "IBI") %>%
    distinct()

  writelog('*** DATA *** Final IBI Metrics - IBI-final.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(ibi_metrics, logfile = file.path(dirname(logfile), 'IBI-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  writelog('\nEND: IBI function.\n', logfile = logfile, verbose = verbose)

  return(ibi_metrics)
}

