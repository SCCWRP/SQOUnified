# GENERIC RIVPACS (Alt) -----------------------------------------------------------------------
#' Compute the RIVPACS (River Invertebrate Prediction and Classification System) index (generic version).
#'
#' @description
#'   The RIVPACS index is an Observed-to-Expected (O/E) index comparing the infaunal community observed
#'   in a sample to that which would be expected at the location under reference conditions, based on a
#'   predictive discriminant function model using latitude, longitude, and water depth.
#'
#' @details
#'   The RIVPACS model assigns each sample to reference site clusters based on abiotic predictors
#'   (latitude, longitude, depth) using Mahalanobis distance. Taxa expected at each site are predicted
#'   from reference community composition weighted by cluster membership probabilities. The O/E ratio
#'   compares observed richness to expected richness.
#'
#'   Samples flagged as outliers (outside the chi-squared critical value at the 0.01 level) receive
#'   a caution note indicating they may be outside the RIVPACS habitat model calibration range.
#'
#'   RIVPACS condition categories:
#'   \itemize{
#'     \item Reference: 0.90 < score < 1.10 (category score 1)
#'     \item Low Disturbance: 0.74 < score <= 0.90 or 1.10 <= score < 1.26 (category score 2)
#'     \item Moderate Disturbance: 0.32 < score <= 0.74 or score >= 1.26 (category score 3)
#'     \item High Disturbance: score <= 0.32 (category score 4)
#'   }
#'
#'   Details on the specifics of the calculation can be found in Bay et al. 2021. Sediment Quality
#'   Assessment Technical Support Manual. SCCWRP Technical Report 777.
#'
#'   Background concepts can be found in Van Sickle et al. 2006. Selecting discriminant function models
#'   for predicting the expected richness of aquatic macroinvertebrates. Freshwater Biology 51:359-372.
#'
#' @param BenthicData a data frame containing benthic data and station information with at minimum:
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
#'    \strong{\code{exclude}} - "Yes" or "No" indicating if the taxon name is ambiguous;
#'
#'    \strong{\code{latitude}} - latitude in decimal degrees;
#'
#'    \strong{\code{longitude}} - longitude in decimal degrees (negative for west);
#'
#'    \strong{\code{depth}} - station depth in meters.
#'
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @usage
#' alt.RIVPACS.generic(BenthicData)
#'
#' @examples
#' \dontrun{
#'   alt.RIVPACS.generic(my_benthic_data)
#' }
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace_all
#' @importFrom tibble column_to_rownames as_tibble
#'
#' @export
alt.RIVPACS.generic <- function(BenthicData,
                                logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'RIVPACS_generic_log.Rmd'),
                                verbose = FALSE,
                                knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n### BEGIN: Generic RIVPACS function.\n', logfile = logfile, verbose = verbose)

  # Reference data (socal.reference.*, socal.example.*) is available via R/sysdata.rda
  # SoCalRivpacs.2() is available in the package namespace from R/SoCalRivpacs2.R

  #rectifying field and dataframe naming conventions to match pre-existing RIVPACS code
  benthic_data <- BenthicData %>%
    rename(Taxa = taxon) %>%
    mutate(sample_id = paste(stationid, sampledate, replicate, sep = "_")) %>%
    filter(Taxa != "NoOrganismsPresent") #removing defaunated samples from analysis

  #selecting sample id information to associate w/ OE scores later
  sampleids <- benthic_data %>%
    select(sample_id, stationid, sampledate, replicate) %>%
    distinct()

  #in case a sample had no animals, force it into the High Disturbance category
  defaunated <- BenthicData %>%
    filter(taxon == "NoOrganismsPresent") %>%
    mutate(index = "Rivpacs", score = NaN, condition_category = "High Disturbance", condition_category_score = 4, note = "Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)

  # selecting location and depth information for discriminant function model
  scb.predictors <- benthic_data %>%
    select(sample_id, Latitude = latitude, Longitude = longitude, SampleDepth = depth) %>%
    distinct() %>%
    column_to_rownames("sample_id") %>%
    as.matrix()

  # selecting the observed taxa in each sample
  scb.taxa <- benthic_data %>%
    filter(exclude == "No") %>%
    select(sample_id, Taxa, Abundance = abundance) %>%
    mutate(Taxa = str_replace_all(Taxa, "[ ()]", "_")) # changing naming convention to match discriminant function

  #Putting the taxa data into a taxon-abundance matrix
  scb.taxa.2 <- scb.taxa %>%
    group_by(sample_id, Taxa) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    pivot_wider(id_cols = sample_id, names_from = Taxa,
                values_from = Abundance, values_fill = 0) %>%
    column_to_rownames("sample_id")

  writelog(
    '#### RIVPACS Step 1 - Predictors and taxa matrices prepared\n',
    logfile = logfile,
    verbose = verbose
  )

  #submitting abiotic predictors and biotic response data to the RIVPACS discriminant function
  rivpacs.mod <- SoCalRivpacs.2(observed.predictors = scb.predictors, observed.taxa = scb.taxa.2)

  # extract the observed taxa and expected taxa information from the RIVPACS function output
  oe.tab <- rivpacs.mod$oe.table %>%
    as_tibble() %>%
    left_join(., sampleids, by = "sample_id")

  # gathering the station information associated with each sample
  oe.stations <- BenthicData %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()

  # Assign condition categories based upon O:E scores and thresholds
  rivpacs.scores <- oe.tab %>%
    rename(score = O.over.E) %>%
    select(-sample_id) %>%
    mutate(
      note = case_when(outlier.01 == 'FAIL' ~ "Caution-Sample Outside RIVPACS habitat model",
                       TRUE ~ NA),
      index = "Rivpacs") %>%
    relocate(., stationid, sampledate, replicate, O, E, index, score, outlier.05, outlier.01, note)

  writelog(
    '#### RIVPACS Step 2 - O/E details\n',
    logfile = logfile,
    data = rivpacs.scores %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rivpacs.scores, logfile = logfile, filename = 'RIVPACS_generic-step1-OE_details.csv', linktext = 'Download RIVPACS O/E details', verbose = verbose)

  #clean up scores and add associated station information
  rivpacs.scores.2 <- rivpacs.scores %>%
    mutate(condition_category = case_when((score <= 0.32) ~ "High Disturbance",
                                          ((score > 0.32 & score <= 0.74) | (score >= 1.26)) ~ "Moderate Disturbance",
                                          ((score > 0.74 & score <= 0.90) | score >= 1.10 & score < 1.26) ~ "Low Disturbance",
                                          (score > 0.90 | score < 1.10) ~ "Reference"),
           condition_category_score = case_when(condition_category == "Reference" ~ 1,
                                                condition_category == "Low Disturbance" ~ 2,
                                                condition_category == "Moderate Disturbance" ~ 3,
                                                condition_category == "High Disturbance" ~ 4)) %>%
    bind_rows(defaunated, .) %>%
    select(-O, -E, -outlier.05, -outlier.01) %>%
    left_join(oe.stations, ., by = c("stationid", "sampledate", "replicate"))

  writelog(
    '#### RIVPACS Final - RIVPACS Scores\n',
    logfile = logfile,
    data = rivpacs.scores.2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rivpacs.scores.2, logfile = logfile, filename = 'RIVPACS_generic-final_scores.csv', linktext = 'Download RIVPACS scores', verbose = verbose)

  writelog('\n### END: Generic RIVPACS function.\n', logfile = logfile, verbose = verbose)

  return(rivpacs.scores.2)
}
