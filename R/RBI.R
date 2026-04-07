# GENERIC RBI (Alt) ---------------------------------------------------------------------------
#' Compute the Relative Benthic Index (RBI) score (generic version).
#'
#' @description
#'   The RBI is a multi-metric index calculated from the weighted sum of: (a) four community metrics related to
#'   biodiversity (total number of taxa, number of crustacean taxa, abundance of crustacean individuals, and number
#'   of mollusc taxa), (b) abundances of three positive indicator taxa, and (c) the presence of two negative
#'   indicator species.
#'
#' @details
#'   The metrics used in the RBI are:
#'   \enumerate{
#'     \item Total number of taxa - a measure of biodiversity
#'     \item Number of mollusc taxa - measure of taxa sensitive to eutrophication
#'     \item Number of crustacean taxa - measure of taxa sensitive to pesticides and low DO
#'     \item Crustacean abundance
#'     \item Number of individuals of \emph{Monocorophium insidiosum} - pollution sensitive amphipod
#'     \item Number of individuals of \emph{Asthenothaerus diegensis} - pollution sensitive bivalve
#'     \item Number of individuals of \emph{Goniada littorea} - pollution sensitive polychaete
#'     \item Presence of \emph{Capitella capitata} complex and Oligochaeta - pollution indicative taxa
#'   }
#'
#'   Observed values are scaled relative to values observed at reference sites in the southern California
#'   calibration data set. Scaled metrics are combined into three meta-metrics: Taxa Richness Weighted Value (TWV),
#'   Positive Indicator Taxa (PIT), and Negative Indicator Taxa (NIT). The final RBI score is:
#'
#'   \deqn{ RBI = \frac{(TWV + NIT + 2 \times PIT) - 0.03}{4.69} }
#'
#'   RBI condition categories:
#'   \itemize{
#'     \item Reference: score > 0.27 (category score 1)
#'     \item Low Disturbance: 0.16 < score <= 0.27 (category score 2)
#'     \item Moderate Disturbance: 0.09 <= score <= 0.16 (category score 3)
#'     \item High Disturbance: score < 0.09 (category score 4)
#'   }
#'
#'   Details on the specifics of the calculation of the index can be found in Bay et al. 2021. Sediment Quality
#'   Assessment Technical Support Manual. SCCWRP Technical Report 777.
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
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @usage
#' RBI(benthic_data)
#'
#' @examples
#' \dontrun{
#'   RBI(my_benthic_data)
#' }
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom lubridate ymd
#'
#' @export
RBI <- function(benthic_data,
                            logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'RBI_generic_log.Rmd'),
                            verbose = FALSE,
                            knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n### BEGIN: Generic RBI function.\n', logfile = logfile, verbose = verbose)

  # Reference data (xl_tool.SoCalLUList) is available via R/sysdata.rda

  #create an empty dataframe to populate with RBI scores
  rbi.out.null <- tibble(stationid = "dummy",
                         replicate = NaN,
                         sampledate = ymd("2000/01/1"),
                         index = "RBI",
                         score = NaN,
                         condition_category = NA,
                         condition_category_score = NA,
                         note = NA)

  #in case a sample had no animals, force it into the High Disturbance category
  defaunated <- benthic_data %>%
    filter(taxon == "NoOrganismsPresent") %>%
    mutate(index = "RBI", score = NaN, condition_category = "High Disturbance", condition_category_score = 4, note = "Defaunated Sample") %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score, note)


  # Prepare the given data frame so that we can compute the RBI score and categories
  rbi_data <- benthic_data %>%
    select(stationid, replicate, sampledate, taxon, abundance, exclude) %>%
    left_join(xl_tool.SoCalLUList, by = c("taxon" = "TaxonName")) %>%
    filter(taxon != "NoOrganismsPresent")

  #Export data so the user knows what is going to be used in subsequent calculations
  rbi_data.review <- rbi_data %>%
    mutate(rich_flag = case_when(Phylum == "" ~ 0,
                                 is.na(Phylum) ~ 0,
                                 TRUE ~ 1),
           pit_flag = if_else(taxon %in% c("Monocorophium insidiosum", "Asthenothaerus diegensis", "Goniada littorea"), 1, 0),
           nit_flag = if_else(taxon %in% c("Capitella capitata Cmplx", "Oligochaeta"), 1, 0)) %>%
    select(stationid, sampledate, replicate, taxon, abundance, exclude, SpeciesLevel, Mollusc, Crustacean, rich_flag, pit_flag, nit_flag)

  writelog(
    '#### RBI Step 1 - Data to be analyzed with SQO designations\n',
    logfile = logfile,
    data = rbi_data.review %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_data.review, logfile = logfile, filename = 'RBI_generic-step1-data_with_designations.csv', linktext = 'Download RBI data with SQO designations', verbose = verbose)


  ####Calculate the different metrics for the RBI
  # calculate taxa richness
  rbi1 <- rbi_data %>%
    filter(exclude == "No") %>%
    mutate(rich_flag = case_when(Phylum == "" ~ 0,
                                 is.na(Phylum) ~ 0,
                                 TRUE ~ 1)) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfTaxa = sum(rich_flag), .groups = "drop_last") %>%
    ungroup()

  # calculate mollusc taxa richness
  rbi2 <- rbi_data %>%
    filter(exclude == "No") %>%
    mutate(flag = (case_when(Mollusc == "Mollusc" ~ 1,
                             TRUE ~ 0))) %>%
    group_by(stationid, sampledate, replicate) %>%
    summarise(NumOfMolluscTaxa = sum(flag), .groups = "drop_last") %>%
    ungroup()

  # calculate crustacean richness
  rbi3 <- rbi_data %>%
    filter(exclude == "No") %>%
    mutate(flag = case_when(Crustacean == "Crustacean" ~ 1,
                            TRUE ~ 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(NumOfCrustaceanTaxa = sum(flag), .groups = "drop_last") %>%
    ungroup()

  # calculate crustacean abundance
  rbi4 <- rbi_data %>%
    mutate(flag = case_when(Crustacean == "Crustacean" ~ abundance,
                            TRUE ~ 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(CrustaceanAbun = sum(flag), .groups = "drop_last") %>%
    ungroup()

  # calculate abundance of M. insidiosum
  rbi5 <- rbi_data %>%
    mutate(flag = case_when(taxon == "Monocorophium insidiosum" ~ abundance,
                            TRUE ~ 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(M_insidiosumAbun = sum(flag), .groups = "drop_last") %>%
    ungroup()

  # calculate abundance of A. diegensis
  rbi6 <- rbi_data %>%
    mutate(flag = case_when(taxon == "Asthenothaerus diegensis" ~ abundance,
                            TRUE ~ 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(A_diegensisAbun = sum(flag), .groups = "drop_last") %>%
    ungroup()

  # calculate abundance of G. littorea
  rbi7 <- rbi_data %>%
    mutate(flag = case_when(taxon == "Goniada littorea" ~ abundance,
                            TRUE ~ 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(G_littoreaAbun = sum(flag), .groups = "drop_last") %>%
    ungroup()

  # calculate negative indicator taxa score (NIT)
  rbi8 <- rbi_data %>%
    mutate(badness = case_when(taxon %in% c("Capitella capitata Cmplx", "Oligochaeta") ~ -0.1,
                               TRUE ~ 0)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(NIT = sum(badness), .groups = "drop_last") %>%
    ungroup()

  # combine all the metrics into one data frame
  rbi_metrics <- rbi1 %>%
    dplyr::full_join(rbi2, by = c("stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi3, by = c("stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi4, by = c("stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi5, by = c("stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi6, by = c("stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi7, by = c("stationid", "replicate", "sampledate")) %>%
    dplyr::full_join(rbi8, by = c("stationid", "replicate", "sampledate"))

  writelog(
    '#### RBI Step 2 - Raw RBI metrics\n',
    logfile = logfile,
    data = rbi_metrics %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_metrics, logfile = logfile, filename = 'RBI_generic-step2-raw_metrics.csv', linktext = 'Download raw RBI metrics', verbose = verbose)

  # Scale observed values to the maxima observed in index calibration dataset
  rbi_scaled <- rbi_metrics %>%
    mutate(scaled_NumTaxa = (NumOfTaxa / 99),
           scaled_NumMolluscTaxa = (NumOfMolluscTaxa / 28),
           scaled_NumCrustaceanTaxa = (NumOfCrustaceanTaxa / 29),
           scaled_CrustaceanAbun = (CrustaceanAbun / 1693),
           # calculate Taxa Richness Weighted Value (TWV)
           TWV = (scaled_NumTaxa + scaled_NumMolluscTaxa + scaled_NumCrustaceanTaxa + (0.25 * scaled_CrustaceanAbun)),
           # scale the indicator taxa
           scaled_M_insid = (((M_insidiosumAbun)^(1/4)) / ((473)^(1/4))),
           scaled_A_dieg = (((A_diegensisAbun)^(1/4)) / ((27)^(1/4))),
           scaled_G_littor = (((G_littoreaAbun)^(1/4)) / ((15)^(1/4))),
           # calculate Positive Indicator Taxa (PIT)
           PIT = scaled_M_insid + scaled_A_dieg + scaled_G_littor,
           # integrate TWV, NIT, and PIT
           Raw_RBI = TWV + NIT + (2 * PIT),
           # scaling the RBI score
           score = ((Raw_RBI - 0.03) / 4.69),
           # RBI Categories based on RBI scores
           condition_category = case_when((score > 0.27) ~ "Reference",
                                          (score > 0.16 & score <= 0.27) ~ "Low Disturbance",
                                          (score >= 0.09 & score <= 0.16) ~ "Moderate Disturbance",
                                          (score < 0.09) ~ "High Disturbance"),
           # RBI Category Scores
           condition_category_score = case_when((condition_category == "Reference") ~ 1,
                                                (condition_category == "Low Disturbance") ~ 2,
                                                (condition_category == "Moderate Disturbance") ~ 3,
                                                (condition_category == "High Disturbance") ~ 4),
           index = "RBI")

  writelog(
    '#### RBI Step 3 - Scaled RBI metrics\n',
    logfile = logfile,
    data = rbi_scaled %>% select(-condition_category, -condition_category_score) %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_scaled %>% select(-condition_category, -condition_category_score), logfile = logfile, filename = 'RBI_generic-step3-scaled_metrics.csv', linktext = 'Download scaled RBI metrics', verbose = verbose)


  #gathering the station information for each site
  rbi.stations <- benthic_data %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()

  rbi.out <- rbi_scaled %>%
    select(stationid, sampledate, replicate, index, score, condition_category, condition_category_score) %>%
    bind_rows(rbi.out.null, ., defaunated) %>%
    full_join(rbi.stations, ., by = c("stationid", "sampledate", "replicate"))

  if (length(rbi.out$stationid) > 1) {
    rbi.out.2 <- rbi.out %>%
      filter(stationid != "dummy") %>%
      mutate(note = if_else(is.na(note), "none", note))
  } else {
    rbi.out.2 <- rbi.out %>%
      mutate(note = "RBI scores not calculated")
  }

  writelog(
    '#### RBI Final - RBI Scores\n',
    logfile = logfile,
    data = rbi.out.2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi.out.2, logfile = logfile, filename = 'RBI_generic-final_scores.csv', linktext = 'Download RBI scores', verbose = verbose)

  writelog('\n### END: Generic RBI function.\n', logfile = logfile, verbose = verbose)

  return(rbi.out.2)
}
