# GENERIC IBI (Alt) ---------------------------------------------------------------------------
#' Compute the Index of Biotic Integrity (IBI) and IBI condition category (generic version).
#'
#' @description
#'   The IBI is a multi-metric index that compares the values of four different metrics to the ranges expected under
#'   reference conditions. The score increases by one for each metric that is outside of the reference range. The four
#'   metrics are:
#'
#'   (1) the total number of taxa - measure of biodiversity
#'   (2) the total number of mollusc taxa - measure of sensitivity to eutrophication and potentially invasive taxa
#'   (3) the abundance of Notomastus sp. - measure of the presence of organic matter indicative taxa
#'   (4) the number of sensitive taxa - measure of the presence of pollution sensitive taxa designated by Thompson and Lowe, 2004
#'
#' @details
#'   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality Assessment
#'   Technical Support Manual. SCCWRP Technical Report 777.
#'
#'   Details on validation of the index can be found in Ranasinghe et al. 2009 Calibration and evaluation of five indicators
#'   of benthic community condition in two California bay and estuary habitats. Marine Pollution Bulletin 59:5-13.
#'
#'   Background concepts of the index can be found in Thompson and Lowe 2004. Assessment of macrobenthos response to sediment
#'   contamination in the San Francisco Estuary, California USA. Environmental Toxicology and Chemistry 23:2178-2187.
#'
#'   Reference ranges for IBI metrics in Southern California Marine Bays
#'   (Table 4.19 CASQO Technical Manual 3rd edition 2021 - page 68):
#'   \itemize{
#'     \item NumOfTaxa: 13 - 99
#'     \item NumOfMolluscTaxa: 2 - 25
#'     \item NotomastusAbun: 0 - 59
#'     \item PctSensTaxa: 19 - 47.1
#'   }
#'
#'   IBI condition categories:
#'   \itemize{
#'     \item Reference: score = 0 (category score 1)
#'     \item Low Disturbance: score = 1 (category score 2)
#'     \item Moderate Disturbance: score = 2 (category score 3)
#'     \item High Disturbance: score = 3 or 4 (category score 4)
#'   }
#'
#'   This version uses the SQO Excel tool look-up list as decided by the Bight 23 index code subcommittee.
#'
#' @param benthic_data a data frame containing benthic data and station information with at minimum:
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the sampling location;
#'
#'    \strong{\code{replicate}} - a numeric identifying the replicate number;
#'
#'    \strong{\code{sampledate}} - the date of sample collection;
#'
#'    \strong{\code{taxon}} - name of the organism. Use \code{NoOrganismsPresent} with 0 abundance for empty samples;
#'
#'    \strong{\code{abundance}} - number of individuals counted;
#'
#'    \strong{\code{exclude}} - "Yes" or "No" indicating if the taxon name is ambiguous.
#'
#' @param retrofit_taxonomy Logical. If TRUE (default), modern SCAMIT Edition 14 taxonomy in the
#'    submitted data is retrofitted back to SQO-compatible names via \code{\link{benthicdata_prep}}
#'    before the index is calculated. Set FALSE if the data has already been prepped/retrofitted.
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @usage
#' IBI(benthic_data)
#'
#' @examples
#' \dontrun{
#'   IBI(my_benthic_data)
#' }
#'
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect
#' @importFrom lubridate ymd
#'
#' @export
IBI <- function(benthic_data,
                            retrofit_taxonomy = TRUE,
                            logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'IBI_generic_log.Rmd'),
                            verbose = FALSE,
                            knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n## BEGIN: Generic IBI function.\n', logfile = logfile, verbose = verbose)

  # Standardize and (optionally) retrofit the submitted taxonomy to SQO-compatible names
  benthic_data <- benthicdata_prep(benthic_data, retrofit = retrofit_taxonomy, logfile = logfile, verbose = verbose)$benthic_data

  # Reference data (xl_tool.SoCalLUList) is available as a package dataset — see ?xl_tool.SoCalLUList

  #create an empty dataframe to populate with IBI scores
  ibi.out.null <- tibble(stationid = "dummy",
                         replicate = NaN,
                         sampledate = ymd("2000/01/1"),
                         index = "IBI",
                         score = NaN,
                         condition_category = NA,
                         condition_category_score = NA,
                         note = NA)

  #in case a sample had no animals, force it into the High Disturbance category
  # Collapse any duplicate NoOrganismsPresent rows within a sample so we don't
  # emit multiple defaunated rows per (stationid, sampledate, replicate).
  defaunated <- benthic_data %>%
    filter(taxon == "NoOrganismsPresent") %>%
    distinct(stationid, sampledate, replicate, .keep_all = TRUE) %>%
    mutate(index = "IBI", score = NaN, condition_category = "High Disturbance", condition_category_score = 4, note = "Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)

  # Prepare the given data frame so that we can compute the IBI score and categories
  ibi_data <- benthic_data %>%
    left_join(xl_tool.SoCalLUList, by = c("taxon" = "TaxonName")) %>%
    filter(taxon != "NoOrganismsPresent")

  #Export data so the user knows what is going to be used in subsequent calculations
  ibi_data.review <- ibi_data %>%
    mutate(Notomastus_flag = case_when(str_detect(taxon, "Notomastus") ~ 1,
                                       TRUE ~ 0)) %>%
    select(stationid, sampledate, replicate, taxon, abundance, exclude, SpeciesLevel, Mollusc, IBISensitive, Notomastus_flag)

  writelog(
    '### IBI Step 1 - Data to be analyzed with SQO designations\n',
    logfile = logfile,
    data = ibi_data.review %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_data.review, logfile = logfile, filename = 'IBI_generic-step1-data_with_designations.csv', linktext = 'Download IBI data with SQO designations', verbose = verbose)


  # Calculate taxa richness
  ibi1 <- ibi_data %>%
    filter(exclude == "No") %>%
    mutate(rich_flag = case_when(Phylum == "" ~ 0,
                                 is.na(Phylum) ~ 0,
                                 TRUE ~ 1)) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfTaxa = length(taxon), .groups = "drop_last")


  # Calculate mollusc taxa richness
  ibi2 <- ibi_data %>%
    filter(exclude == "No") %>%
    mutate(flag = (case_when(Mollusc == "Mollusc" ~ 1,
                             TRUE ~ 0))) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfMolluscTaxa = sum(flag), .groups = "drop_last") %>%
    ungroup()


  # calculate Notomastus spp. abundance
  ibi3 <- ibi_data %>%
    mutate(flag = case_when(str_detect(taxon, "Notomastus") ~ abundance,
                            TRUE ~ 0)) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NotomastusAbun = sum(flag), .groups = "drop_last") %>%
    ungroup()


  # Calculate % Sensitive Taxa
  ibi4 <- ibi_data %>%
    mutate(flag = case_when(IBISensitive == "S" ~ 1,
                            TRUE ~ 0)) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(sensitive_S = sum(flag), .groups = "drop_last") %>%
    ungroup() %>%
    left_join(ibi1, by = c("stationid", "sampledate", "replicate")) %>%
    mutate(PctSensTaxa = (sensitive_S / NumOfTaxa) * 100) %>%
    select(stationid, sampledate, replicate, PctSensTaxa)


  ### IBI Metrics:
  ibi_metrics <- ibi1 %>%
    full_join(ibi2, by = c("sampledate", "stationid", "replicate")) %>%
    full_join(ibi3, by = c("sampledate", "stationid", "replicate")) %>%
    full_join(ibi4, by = c("sampledate", "stationid", "replicate"))

  writelog(
    '### IBI Step 2 - IBI metric values\n',
    logfile = logfile,
    data = ibi_metrics %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi_metrics, logfile = logfile, filename = 'IBI_generic-step2-metric_values.csv', linktext = 'Download IBI metric values', verbose = verbose)


  ### Reference ranges for IBI metrics in Southern California Marine Bays
  ### [ Table 4.19 CASQO Technical Manual 3rd edition 2021 - page 68 ]
  ibi_ref_ranges_table <- data.frame(metric = c("NumOfTaxa", "NumOfMolluscTaxa", "NotomastusAbun", "PctSensTaxa"),
                                     ref_low = c(13, 2, 0, 19),
                                     ref_high = c(99, 25, 59, 47.1))


  # Calculate IBI scores
  ibi.scores <- ibi_metrics %>%
    pivot_longer(., cols = c(-stationid, -sampledate, -replicate), names_to = "metric", values_to = "value") %>%
    left_join(., ibi_ref_ranges_table, by = "metric") %>%
    mutate(out_of_range = case_when(value < ref_low | value > ref_high ~ 1,
                                    TRUE ~ 0)) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(score = sum(out_of_range), .groups = "drop_last") %>%
    ungroup() %>%
    mutate(index = "IBI", .before = score) %>%
    mutate(condition_category = case_when(score == 0 ~ "Reference",
                                          score == 1 ~ "Low Disturbance",
                                          score == 2 ~ "Moderate Disturbance",
                                          score %in% c(3, 4) ~ "High Disturbance"),
           condition_category_score = case_when(score == 0 ~ 1,
                                                score == 1 ~ 2,
                                                score == 2 ~ 3,
                                                score %in% c(3, 4) ~ 4))


  #gathering the station information for each site
  ibi.stations <- benthic_data %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct(stationid, sampledate, replicate, .keep_all = TRUE)

  # Drop any samples from defaunated that already have a calculated score, so
  # the bind_rows + join below can never emit two rows per sample.
  defaunated <- defaunated %>%
    anti_join(ibi.scores, by = c("stationid", "sampledate", "replicate"))

  ibi.out <- ibi.scores %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score) %>%
    bind_rows(ibi.out.null, ., defaunated) %>%
    full_join(ibi.stations, ., by = c("stationid", "sampledate", "replicate"))

  if (length(ibi.out$stationid) > 1) {
    ibi.out.2 <- ibi.out %>%
      filter(stationid != "dummy") %>%
      mutate(note = if_else(is.na(note), "none", note))
  } else {
    ibi.out.2 <- ibi.out %>%
      mutate(note = "IBI scores not calculated")
  }

  writelog(
    '### IBI Final - IBI Scores\n',
    logfile = logfile,
    data = ibi.out.2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = ibi.out.2, logfile = logfile, filename = 'IBI_generic-final_scores.csv', linktext = 'Download IBI scores', verbose = verbose)

  writelog('\n## END: Generic IBI function.\n', logfile = logfile, verbose = verbose)

  return(ibi.out.2)
}
