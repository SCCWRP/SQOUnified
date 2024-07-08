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
IBI <- function(BenthicData, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = F)
{

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\n### BEGIN: IBI function.\n', logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'benthic.ibi.input.RData'
  if (verbose) {
    save(BenthicData, file = file.path( dirname(logfile), rawinput.filename ))
  }
  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n#### Raw input to IBI:",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') ### This will load a dataframe called 'BenthicData' into your environment"),
    verbose = verbose
  )
  create_download_link(data = BenthicData, logfile = logfile, filename = 'IBI-RawInput.csv', linktext = 'Download Raw Input to IBI Function', verbose = verbose)

  # SQO List (New)
  writelog(
    "\n#### SQO List (New) which gets joined to raw input",
    logfile = logfile,
    data = sqo.list.new %>% head(25),
    verbose = verbose
  )
  create_download_link(data = sqo.list.new, logfile = logfile, filename = 'sqo.list.new.csv', linktext = 'Download sqo.list.new which gets joined to raw input', verbose = verbose)


  # Prepare the given data frame so that we can compute the IBI score and categories
  ibi_data <- BenthicData %>%
    filter(Exclude!="Yes") %>%
    left_join(sqo.list.new, by = c('Taxon'='TaxonName')) %>%
    mutate_if(is.numeric, list(~na_if(., -88))) %>%
    select('StationID','SampleDate', 'Replicate','Taxon','Abundance','Stratum', 'Phylum', 'IBISensitive', "Mollusc", "Crustacean") %>%
    mutate(n=if_else(Taxon=="NoOrganismsPresent", 0,1))
  # Write to the logs for preparing the given data frame for IBI score and categories
  writelog(
    '\n#### Prepare the given data frame so that we can compute the IBI score and categories',
    logfile = logfile,
    code = "
      ibi_data <- BenthicData %>%
        filter(Exclude != 'Yes') %>%
        left_join(sqo.list.new, by = c('Taxon' = 'TaxonName')) %>%
        mutate_if(is.numeric, list(~na_if(., -88))) %>%
        select('StationID', 'SampleDate', 'Replicate', 'Taxon', 'Abundance', 'Stratum', 'Phylum', 'IBISensitive', 'Mollusc', 'Crustacean') %>%
        mutate(n = if_else(Taxon == 'NoOrganismsPresent', 0, 1))
    ",
    data = ibi_data %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_data, logfile = logfile, filename = 'IBI-prepared-data.csv', linktext = 'Download prepared data for IBI score', verbose = verbose)




  # Group and get NumOfTaxa
  ibi1 <- ibi_data %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NumOfTaxa =sum(n))

  # Write to the logs for grouping and getting NumOfTaxa
  writelog(
    '\n#### Group and get NumOfTaxa',
    logfile = logfile,
    code = "
      ibi1 <- ibi_data %>%
        group_by(Stratum, StationID, SampleDate, Replicate) %>%
        summarise(NumOfTaxa = sum(n))
    ",
    data = ibi1 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi1, logfile = logfile, filename = 'IBI-grouped-data.csv', linktext = 'Download grouped data for IBI score', verbose = verbose)




  # Get Number of Mollusc Taxa
  ibi2 <- ibi_data %>%
    filter(Mollusc == "Mollusc") %>%
    group_by(Stratum, StationID, SampleDate,Replicate) %>%
    summarise(NumOfMolluscTaxa = length(Taxon))

  # Write to the logs for getting Number of Mollusc Taxa
  writelog(
    '\n#### Get Number of Mollusc Taxa',
    logfile = logfile,
    code = "
    ibi2 <- ibi_data %>%
      filter(Mollusc == 'Mollusc') %>%
      group_by(Stratum, StationID, SampleDate, Replicate) %>%
      summarise(NumOfMolluscTaxa = length(Taxon))
  ",
    data = ibi2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi2, logfile = logfile, filename = 'IBI-mollusc-taxa.csv', linktext = 'Download Number of Mollusc Taxa', verbose = verbose)



  # Get Notomastus abundance
  ibi3_2 <- ibi_data %>%
    filter(str_detect(Taxon,"Notomastus")) %>%
    group_by(Stratum, StationID, SampleDate, Replicate) %>%
    summarise(NotomastusAbun = sum(Abundance))

  # Write to the logs for getting Notomastus abundance
  writelog(
    '\n#### Get Notomastus abundance',
    logfile = logfile,
    code = "
    ibi3_2 <- ibi_data %>%
      filter(str_detect(Taxon, 'Notomastus')) %>%
      group_by(Stratum, StationID, SampleDate, Replicate) %>%
      summarise(NotomastusAbun = sum(Abundance))
  ",
    data = ibi3_2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi3_2, logfile = logfile, filename = 'IBI-notomastus-abundance.csv', linktext = 'Download Notomastus abundance', verbose = verbose)




  # Get Percentage of sensitive taxa
  ibi4_2 <- ibi_data %>%
    filter(IBISensitive=="S") %>%
    left_join(ibi1, by=c("Stratum", "StationID", "SampleDate", "Replicate")) %>%
    group_by(Stratum, StationID, SampleDate, Replicate, NumOfTaxa) %>%
    summarise(SensTaxa = length(Taxon)) %>%
    mutate(PctSensTaxa=(SensTaxa/NumOfTaxa)*100) %>%
    select(Stratum, StationID, SampleDate, Replicate, PctSensTaxa)

  # Write to the logs for getting Percentage of sensitive taxa
  writelog(
    '\n#### Get Percentage of sensitive taxa',
    logfile = logfile,
    code = "
      ibi4_2 <- ibi_data %>%
        filter(IBISensitive == 'S') %>%
        left_join(ibi1, by = c('Stratum', 'StationID', 'SampleDate', 'Replicate')) %>%
        group_by(Stratum, StationID, SampleDate, Replicate, NumOfTaxa) %>%
        summarise(SensTaxa = length(Taxon)) %>%
        mutate(PctSensTaxa = (SensTaxa / NumOfTaxa) * 100) %>%
        select(Stratum, StationID, SampleDate, Replicate, PctSensTaxa)
    ",
    data = ibi4_2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi4_2, logfile = logfile, filename = 'IBI-sensitive-taxa.csv', linktext = 'Download Percentage of sensitive taxa', verbose = verbose)



  # Reference ranges for IBI metrics in Southern California Marine Bays
  # [ Table 4.19 CASQO Technical Manual 3rd edition 2021 - page 68 ]
  ibi_ref_ranges_table <- data.frame(ref_low = c(13, 2, 0, 19),
                                     ref_high = c(99, 25, 59, 47.1))
  row.names(ibi_ref_ranges_table) <- c("NumOfTaxa", "NumOfMolluscTaxa", "NotomastusAbun", "PctSensTaxa")

  # Write to the logs for reference ranges for IBI metrics
  writelog(
    '\n#### Reference ranges for IBI metrics in Southern California Marine Bays\n[ Table 4.19 CASQO Technical Manual 3rd edition 2021 - page 68 ]',
    logfile = logfile,
    code = "
      ibi_ref_ranges_table <- data.frame(
        ref_low = c(13, 2, 0, 19),
        ref_high = c(99, 25, 59, 47.1)
      )
      row.names(ibi_ref_ranges_table) <- c('NumOfTaxa', 'NumOfMolluscTaxa', 'NotomastusAbun', 'PctSensTaxa')
    ",
    data = ibi_ref_ranges_table %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_ref_ranges_table, logfile = logfile, filename = 'IBI-ref-ranges.csv', linktext = 'Download reference ranges for IBI metrics', verbose = verbose)




  # IBI category response ranges for Southern California Marine Bays
  # [ Table 4.20 - CASQO Technical Manual page 68-69]
  ibi_category_response_table <- data.frame(ibi_score = as.factor(c(0, 1, 2, 3, 4)),
                                            category = as.factor(c("Reference",
                                                                   "Low Disturbance",
                                                                   "Moderate Disturbance",
                                                                   "High Disturbance",
                                                                   "High Disturbance")),
                                            category_score = as.factor(c(1, 2, 3, 4, 4)))
  # Write to the logs for IBI category response ranges
  writelog(
    '\n#### IBI category response ranges for Southern California Marine Bays\n[ Table 4.20 - CASQO Technical Manual page 68-69]',
    logfile = logfile,
    code = "
      ibi_category_response_table <- data.frame(
        ibi_score = as.factor(c(0, 1, 2, 3, 4)),
        category = as.factor(c('Reference', 'Low Disturbance', 'Moderate Disturbance', 'High Disturbance', 'High Disturbance')),
        category_score = as.factor(c(1, 2, 3, 4, 4))
      )
    ",
    data = ibi_category_response_table,
    verbose = verbose
  )
  create_download_link(data = ibi_category_response_table, logfile = logfile, filename = 'IBI-category-response-ranges.csv', linktext = 'Download IBI category response ranges', verbose = verbose)




  # IBI Metrics:
  # We stitch together all the necessary IBI metrics to determine the IBI index.
  # Each of the metrics is then compared to the tables listed above (Table 5.4 and Table 5.5) to determine the IBI score,
  # the IBI Category, and IBI Category Score
  ibi_metrics1 <- ibi1 %>%
    full_join(ibi2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    full_join(ibi3_2, by = c("Stratum", "SampleDate", "StationID", "Replicate")) %>%
    full_join(ibi4_2, by = c("Stratum", "SampleDate", "StationID", "Replicate"))

  # Write to the logs for stitching together IBI metrics
  writelog(
    '\n#### IBI Metrics:\nWe stitch together all the necessary IBI metrics to determine the IBI index.\nEach of the metrics is then compared to the tables listed above (Table 5.4 and Table 5.5) to determine the IBI score,\nthe IBI Category, and IBI Category Score',
    logfile = logfile,
    code = "
    ibi_metrics1 <- ibi1 %>%
      full_join(ibi2, by = c('Stratum', 'SampleDate', 'StationID', 'Replicate')) %>%
      full_join(ibi3_2, by = c('Stratum', 'SampleDate', 'StationID', 'Replicate')) %>%
      full_join(ibi4_2, by = c('Stratum', 'SampleDate', 'StationID', 'Replicate'))
  ",
    data = ibi_metrics1 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_metrics1, logfile = logfile, filename = 'IBI-metrics.csv', linktext = 'Download IBI metrics', verbose = verbose)



  # Replace NA with 0
  ibi_metrics2 <- ibi_metrics1 %>%
    replace(.,is.na(.),0)

  # Write to the logs for replacing NA with 0 in IBI metrics
  writelog(
    '\n#### Replace NA with 0 in IBI metrics',
    logfile = logfile,
    code = "
      ibi_metrics2 <- ibi_metrics1 %>%
        replace(., is.na(.), 0)
    ",
    data = ibi_metrics2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_metrics2, logfile = logfile, filename = 'IBI-metrics-replace-na.csv', linktext = 'Download IBI metrics with NA replaced', verbose = verbose)

  # Calculate Scores and Categorize
  ibi_final <- ibi_metrics2 %>%

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

  # Write to the logs for calculating scores and categorizing IBI
  writelog(
    '\n#### Calculate Scores and Categorize for IBI',
    logfile = logfile,
    code = "
      ibi_final <- ibi_metrics2 %>%

        # We replace any NAs with 0 so that we can compare the values to the tables listed above
        # The IBI score is set to zero before comparison the reference range.
        mutate(Score = 0) %>%
        # For each metric that is out of the reference range (above or below), the IBI score goes up by one.
        mutate(Score = if_else((NumOfTaxa < ibi_ref_ranges_table['NumOfTaxa',]$ref_low  | NumOfTaxa > ibi_ref_ranges_table['NumOfTaxa',]$ref_high),
                               Score + 1, Score)) %>%
        mutate(Score = if_else((NumOfMolluscTaxa < ibi_ref_ranges_table['NumOfMolluscTaxa',]$ref_low  | NumOfMolluscTaxa > ibi_ref_ranges_table['NumOfMolluscTaxa',]$ref_high),
                               Score + 1, Score)) %>%
        mutate(Score = if_else((NotomastusAbun < ibi_ref_ranges_table['NotomastusAbun',]$ref_low  | NotomastusAbun > ibi_ref_ranges_table['NotomastusAbun',]$ref_high),
                               Score + 1, Score)) %>%
        mutate(Score = if_else((PctSensTaxa < ibi_ref_ranges_table['PctSensTaxa',]$ref_low  | PctSensTaxa > ibi_ref_ranges_table['PctSensTaxa',]$ref_high),
                               Score + 1, Score)) %>%
        # The IBI score is then compared to condition category response ranges (Table 5.5) to determine the IBI category and category score.
        mutate(Category = case_when(Score == 0 ~ 'Reference', Score == 1 ~ 'Low Disturbance', Score == 2 ~ 'Moderate Disturbance', (Score == 3 | Score == 4) ~ 'High Disturbance')) %>%
        mutate(`Category Score` = case_when(Score == 0 ~ 1, Score == 1 ~ 2, Score == 2 ~ 3, (Score == 3 | Score == 4) ~ 4)) %>%
        mutate(Index = 'IBI') %>%
        distinct()
    ",
    data = ibi_final %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_final, logfile = logfile, filename = 'IBI-final-scores.csv', linktext = 'Download final IBI scores', verbose = verbose)



  writelog('\n### END: IBI function.\n', logfile = logfile, verbose = verbose)

  return(ibi_final)
}

