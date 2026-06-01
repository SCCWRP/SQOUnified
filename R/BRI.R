# GENERIC SQO BRI (Alt) -----------------------------------------------------------------------
#' Compute the benthic response index (BRI) score and BRI condition category (generic version).
#'
#' @description
#'   The BRI is the 4th root abundance weighted pollution tolerance score of the organisms present in a sample
#'   relative to the 4th root abundance of all taxa in the sample with a tolerance score (aka p-code, p-value).
#'   The higher the BRI score, the more degraded the sample.
#'
#' @details
#'   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality
#'   Assessment Technical Support Manual. SCCWRP Technical Report 777.
#'
#'   Details on validation of the index can be found in Ranasinghe et al. 2009 Calibration and evaluation of five
#'   indicators of benthic community condition in two California bay and estuary habitats. Marine Pollution
#'   Bulletin 59:5-13.
#'
#'   Background concepts of the index can be found in Smith et al 2001. Benthic Response Index for Assessing
#'   Infaunal Communities on the Southern California Mainland Shelf. Ecological Applications 11: 1073-1087.
#'
#'   The BRI score is calculated as:
#'   \deqn{ BRI = \frac{\sum \left(\sqrt[4]{\textrm{Abundance}} \times P\right)}{\sum \sqrt[4]{\textrm{Abundance}}} }
#'
#'   BRI condition categories for Southern California Marine Bays:
#'   \itemize{
#'     \item Reference: score < 39.96 (category score 1)
#'     \item Low Disturbance: 39.96 <= score < 49.15 (category score 2)
#'     \item Moderate Disturbance: 49.15 <= score < 73.27 (category score 3)
#'     \item High Disturbance: score >= 73.27 (category score 4)
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
#' BRI(benthic_data)
#'
#' @examples
#' \dontrun{
#'   BRI(my_benthic_data)
#' }
#'
#' @import dplyr
#' @importFrom tidyr drop_na
#' @importFrom lubridate ymd
#'
#' @export
BRI <- function(benthic_data,
                                retrofit_taxonomy = TRUE,
                                logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'BRI_generic_log.Rmd'),
                                verbose = FALSE,
                                knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n### BEGIN: Generic SQO BRI function.\n', logfile = logfile, verbose = verbose)

  # Standardize and (optionally) retrofit the submitted taxonomy to SQO-compatible names
  benthic_data <- benthicdata_prep(benthic_data, retrofit = retrofit_taxonomy, logfile = logfile, verbose = verbose)$benthic_data

  # Reference data (xl_tool.SoCalLUList) is available via R/sysdata.rda

  #create empty dataframe to populate w/ bri scores
  bri.out.null <- tibble(stationid = "dummy",
                         replicate = NaN,
                         sampledate = ymd("2000-01-1"),
                         index = "BRI",
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
    mutate(index = "BRI", score = NaN, condition_category = "High Disturbance", condition_category_score = 4, note = "Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)

  #matching p codes to taxa in the submitted data
  all.for.bri <- benthic_data %>%
    filter(taxon != "NoOrganismsPresent") %>%
    left_join(., select(xl_tool.SoCalLUList, TaxonName, ToleranceScore), by = c('taxon' = 'TaxonName'))

  #identify the taxa in the submitted data with a tolerance score
  taxa_w_pvalue <- all.for.bri %>%
    group_by(taxon, ToleranceScore) %>%
    summarise(Freq_of_Occ = length(stationid), .groups = "drop_last") %>%
    ungroup() %>%
    drop_na(ToleranceScore)

  writelog(
    '#### BRI Step 1 - Taxa with a tolerance score\n',
    logfile = logfile,
    data = taxa_w_pvalue %>% head(25),
    verbose = verbose
  )
  create_download_link(data = taxa_w_pvalue, logfile = logfile, filename = 'BRI_generic-step1-taxa_with_tolerance.csv', linktext = 'Download taxa with tolerance scores', verbose = verbose)

  #identify those taxa in the submitted data without a tolerance score
  taxa_wo_pvalue <- all.for.bri %>%
    group_by(taxon, ToleranceScore) %>%
    summarise(Freq_of_Occ = length(stationid), .groups = "drop_last") %>%
    ungroup() %>%
    filter(is.na(ToleranceScore)) %>%
    select(-ToleranceScore)

  writelog(
    '#### BRI Step 2 - Taxa without a tolerance score\n',
    logfile = logfile,
    data = taxa_wo_pvalue %>% head(25),
    verbose = verbose
  )
  create_download_link(data = taxa_wo_pvalue, logfile = logfile, filename = 'BRI_generic-step2-taxa_without_tolerance.csv', linktext = 'Download taxa without tolerance scores', verbose = verbose)


  bri.out <- all.for.bri %>%
    drop_na(ToleranceScore) %>%
    mutate(fourthroot_abun = abundance ** 0.25,
           tolerance_value = fourthroot_abun * ToleranceScore) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarize(numerator = sum(tolerance_value, na.rm = T), denomenator = sum(fourthroot_abun, na.rm = T), score = numerator / denomenator, .groups = "drop_last") %>%
    select(stationid, sampledate, replicate, score) %>%
    mutate(
      condition_category = case_when((score < 39.96) ~ "Reference",
                                     (score >= 39.96 & score < 49.15) ~ "Low Disturbance",
                                     (score >= 49.15 & score < 73.27) ~ "Moderate Disturbance",
                                     (score >= 73.27) ~ "High Disturbance"
      )) %>%
    mutate(
      condition_category_score = case_when((condition_category == "Reference") ~ 1,
                                           (condition_category == "Low Disturbance") ~ 2,
                                           (condition_category == "Moderate Disturbance") ~ 3,
                                           (condition_category == "High Disturbance") ~ 4),
      note = NA)

  bri.stations <- benthic_data %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct(stationid, sampledate, replicate, .keep_all = TRUE)

  # Drop any samples from defaunated that already have a calculated score, so
  # the bind_rows + join below can never emit two rows per sample.
  defaunated <- defaunated %>%
    anti_join(bri.out, by = c("stationid", "sampledate", "replicate"))

  bri.out.2 <- bri.out.null %>%
    bind_rows(bri.out, defaunated) %>%
    left_join(bri.stations, ., by = c("stationid", "sampledate", "replicate")) %>%
    mutate(note = case_when(is.na(condition_category) ~ "Unable to calculate BRI, likely no p-code taxa",
                            TRUE ~ note),
           index = "BRI")


  if (length(bri.out.2$stationid) > 1) {
    bri.out.3 <- bri.out.2 %>%
      filter(stationid != "dummy")
  } else {
    bri.out.3 <- bri.out.2 %>%
      mutate(note = "BRI scores not calculated")
  }

  writelog(
    '#### BRI Final - BRI Scores\n',
    logfile = logfile,
    data = bri.out.3 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = bri.out.3, logfile = logfile, filename = 'BRI_generic-final_scores.csv', linktext = 'Download BRI scores', verbose = verbose)

  writelog('\n### END: Generic SQO BRI function.\n', logfile = logfile, verbose = verbose)

  return(bri.out.3)
}
