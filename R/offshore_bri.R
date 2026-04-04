# OFFSHORE BRI (Edition 14) -------------------------------------------------------------------
#' Compute the offshore Benthic Response Index (BRI) score using SCAMIT Edition 14 tolerance values.
#'
#' @description
#'   Calculates BRI scores for offshore benthic samples using depth-zone-specific pollution tolerance
#'   values (p-codes) from SCAMIT Edition 14. Unlike the bay/estuary BRI, this version accounts for
#'   depth-dependent community composition by applying separate tolerance scores in shallow (<25m),
#'   mid (35-110m), and deep (130-324m) depth zones. Samples in overlapping depth ranges (25-35m
#'   and 110-130m) receive the average of the two adjacent zone scores.
#'
#' @details
#'   The offshore BRI is calculated as the cube root (rather than 4th root) abundance-weighted pollution
#'   tolerance score:
#'
#'   \deqn{ BRI = \frac{\sum \left(\sqrt[3]{\textrm{Abundance}} \times P_{dz}\right)}{\sum \sqrt[3]{\textrm{Abundance}}} }
#'
#'   where \eqn{P_{dz}} is the depth-zone-specific tolerance value for each taxon.
#'
#'   Depth zones and their BRI score derivation:
#'   \itemize{
#'     \item Shallow: depth < 25m (shallow p-values only)
#'     \item Shallow-Mid overlap: 25-35m (average of shallow and mid scores)
#'     \item Mid: 35-110m (mid p-values only)
#'     \item Mid-Deep overlap: 110-130m (average of mid and deep scores)
#'     \item Deep: 130-324m (deep p-values only)
#'     \item Samples deeper than 324m are flagged as out of range
#'   }
#'
#'   Condition categories based on BRI score:
#'   \itemize{
#'     \item Reference: score < 25 (class 1)
#'     \item Marginal Deviation: 25 <= score < 34 (class 2)
#'     \item Biodiversity Loss: 34 <= score < 44 (class 3)
#'     \item Function Loss: 44 <= score < 72 (class 4)
#'     \item Defaunation: score >= 72 (class 5)
#'   }
#'
#'   This function was developed based on Smith et al. (2001) and refined by the Bight 2023
#'   index code sub-committee. Version 2.0 calculates BRI scores for overlapping depth zones
#'   by averaging scores from each adjacent depth zone.
#'
#'   Geographic warnings are issued for samples north of Point Conception (lat > 34.45) or
#'   south of the US-Mexico border (lat < 32.52), as the BRI is nominally applicable only
#'   to the Southern California Bight.
#'
#' @param BenthicData a data frame containing infauna abundance data with the following columns:
#'
#'    \strong{\code{StationID}} - an alpha-numeric identifier of the sampling location;
#'
#'    \strong{\code{SampleDate}} - the date of sample collection;
#'
#'    \strong{\code{Replicate}} - a numeric identifying the replicate number;
#'
#'    \strong{\code{Taxon}} - name of the organism (SCAMIT Edition 14 naming conventions).
#'        If no organisms were present, use \code{NoOrganismsPresent} with 0 abundance;
#'
#'    \strong{\code{Abundance}} - number of individuals counted for the specified taxon.
#'
#' @param StationData a data frame containing station information with the following columns:
#'
#'    \strong{\code{StationID}} - an alpha-numeric identifier matching \code{BenthicData};
#'
#'    \strong{\code{Depth}} - station depth in meters;
#'
#'    \strong{\code{Latitude}} - station latitude in decimal degrees;
#'
#'    \strong{\code{Longitude}} - station longitude in decimal degrees (negative for west).
#'
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @usage
#' BRI.Offshore(BenthicData, StationData)
#'
#' @examples
#' data(offshore_bri_sampledata) # loads offshore_infauna_sampledata and offshore_station_sampledata
#' BRI.Offshore(offshore_infauna_sampledata, offshore_station_sampledata)
#'
#' @import dplyr
#' @importFrom tidyr drop_na pivot_longer
#' @importFrom stats na.omit
#'
#' @export
BRI.Offshore <- function(BenthicData,
                         StationData,
                         logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'BRI_Offshore_log.Rmd'),
                         verbose = FALSE,
                         knitlog = FALSE)
{

  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n### BEGIN: Offshore BRI (Edition 14) function.\n', logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'offshore.bri.input.RData'
  if (verbose) {
    save(BenthicData, StationData, file = file.path(dirname(logfile), rawinput.filename))
  }

  writelog(
    '\nOffshore BRI calculation based on Smith et al. (2001), refined by the Bight 2023 index code sub-committee.',
    logfile = logfile,
    verbose = verbose
  )

  # Reference data (ed.14.ptaxa, pcodes) is available via R/sysdata.rda — loaded into
  # the package namespace automatically when the package is loaded.

  # ---- Standardize column names to lowercase ----
  names(BenthicData) <- tolower(names(BenthicData))
  names(StationData) <- tolower(names(StationData))

  # ---- Prep the data ----
  station_info <- StationData %>%
    mutate(
      stationid = as.character(stationid),
      depth = as.numeric(depth)
    )

  infauna <- BenthicData %>%
    mutate(
      stationid = as.character(stationid),
      abundance = as.numeric(abundance),
      sample_id = paste(stationid, sampledate, replicate, sep = "_")
    )

  # Join taxa names and counts to depth and geographic information
  taxa_to_calc <- infauna %>%
    select(sample_id, stationid, replicate, sampledate, taxon, abundance) %>%
    left_join(
      station_info %>% select(stationid, depth, latitude, longitude),
      by = "stationid"
    ) %>%
    arrange(sample_id, desc(abundance))

  writelog(
    '#### Offshore BRI Step 1 - Taxa joined to station info\n',
    logfile = logfile,
    data = taxa_to_calc %>% head(25),
    verbose = verbose
  )
  create_download_link(data = taxa_to_calc, logfile = logfile, filename = 'BRI_offshore-step1-taxa_to_analyze.csv', linktext = 'Download taxa to be analyzed', verbose = verbose)


  # ---- Match taxa to p-codes ----
  all.4.bri <- taxa_to_calc %>%
    left_join(ed.14.ptaxa, by = "taxon") %>%
    left_join(pcodes, by = "p_code") %>%
    select(-shallow_mid, -mid_deep)

  # Taxa with assigned p-codes
  with.pcode <- all.4.bri %>%
    distinct(taxon, p_code) %>%
    tidyr::drop_na(p_code) %>%
    arrange(taxon)

  writelog(
    '#### Offshore BRI Step 2 - Taxa with assigned p-codes\n',
    logfile = logfile,
    data = with.pcode %>% head(25),
    verbose = verbose
  )
  create_download_link(data = with.pcode, logfile = logfile, filename = 'BRI_offshore-step2-taxa_with_pcode.csv', linktext = 'Download taxa with p-codes', verbose = verbose)

  # Taxa without p-codes
  no.pcode <- all.4.bri %>%
    distinct(taxon, p_code) %>%
    filter(is.na(p_code)) %>%
    select(-p_code) %>%
    arrange(taxon)

  writelog(
    '#### Offshore BRI Step 3 - Taxa without p-codes\n',
    logfile = logfile,
    data = no.pcode %>% head(25),
    verbose = verbose
  )
  create_download_link(data = no.pcode, logfile = logfile, filename = 'BRI_offshore-step3-taxa_without_pcode.csv', linktext = 'Download taxa without p-codes', verbose = verbose)

  writelog(
    '#### Offshore BRI Step 4 - All taxa with p-codes by sample\n',
    logfile = logfile,
    data = all.4.bri %>% head(25),
    verbose = verbose
  )
  create_download_link(data = all.4.bri, logfile = logfile, filename = 'BRI_offshore-step4-all_taxa_by_sample.csv', linktext = 'Download all taxa by sample', verbose = verbose)


  # ---- Identify special cases ----

  # Defaunated samples
  defaunated <- all.4.bri %>%
    filter(taxon == "NoOrganismsPresent") %>%
    select(sample_id, stationid, sampledate, replicate) %>%
    mutate(
      bri_score = NaN,
      bri_cond = "Defaunation",
      bri_class = 5
    )

  defaunated.samps <- unique(defaunated$sample_id)

  # Samples too deep for BRI
  too.deep <- all.4.bri %>%
    filter(depth > 324) %>%
    distinct(sample_id, stationid, sampledate, replicate) %>%
    mutate(
      bri_score = NaN,
      bri_cond = "BRI not Applicable",
      bri_class = NaN
    )

  too.deep.samps <- unique(too.deep$sample_id)

  # ---- Calculate depth-zone-specific BRI scores ----

  # Helper: flag samples to exclude from calculation
  flag_exclusions <- function(df) {
    df %>%
      mutate(
        drop_flag = case_when(
          sample_id %in% too.deep.samps ~ 1,
          sample_id %in% defaunated.samps ~ 1,
          TRUE ~ 0
        )
      ) %>%
      filter(drop_flag == 0) %>%
      select(-drop_flag)
  }

  # Helper: calculate BRI score for a given depth zone column
  calc_zone_bri <- function(df, zone_col, score_name) {
    df %>%
      flag_exclusions() %>%
      tidyr::drop_na(!!sym(zone_col)) %>%
      mutate(cube_abun = (abundance)^(1/3)) %>%
      group_by(sample_id, stationid, sampledate, replicate) %>%
      mutate(
        tot_bri_abun = sum(abundance),
        tot_cube_abun = sum(cube_abun),
        tol_score = !!sym(zone_col) * cube_abun
      ) %>%
      ungroup() %>%
      group_by(sample_id, stationid, depth, sampledate, replicate, tot_bri_abun, tot_cube_abun) %>%
      summarise(numerator = sum(tol_score), .groups = "drop_last") %>%
      ungroup() %>%
      mutate(!!score_name := numerator / tot_cube_abun)
  }

  bri_scores.shallow <- calc_zone_bri(all.4.bri, "shallow", "shallow_bri_score")
  bri_scores.mid     <- calc_zone_bri(all.4.bri, "mid", "mid_bri_score")
  bri_scores.deep    <- calc_zone_bri(all.4.bri, "deep", "deep_bri_score")

  writelog(
    '\n#### Offshore BRI Step 5 - Depth-zone BRI scores calculated (shallow, mid, deep)\n',
    logfile = logfile,
    verbose = verbose
  )

  # ---- Combine depth zone scores ----
  join_cols <- c("sample_id", "stationid", "sampledate", "replicate", "depth")

  bri_scores.1 <- bri_scores.shallow %>%
    left_join(bri_scores.mid, by = join_cols, suffix = c("_s", "_m")) %>%
    left_join(bri_scores.deep, by = join_cols, suffix = c("", "_d")) %>%
    mutate(
      bri_score = case_when(
        depth < 25 ~ shallow_bri_score,
        depth >= 25 & depth <= 35 ~ (shallow_bri_score + mid_bri_score) / 2,
        depth > 35 & depth < 110 ~ mid_bri_score,
        depth >= 110 & depth <= 130 ~ (mid_bri_score + deep_bri_score) / 2,
        depth > 130 & depth <= 324 ~ deep_bri_score
      ),
      depth_zone = case_when(
        depth < 25 ~ "shallow",
        depth >= 25 & depth <= 35 ~ "shallow_mid",
        depth > 35 & depth < 110 ~ "mid",
        depth >= 110 & depth <= 130 ~ "mid_deep",
        depth > 130 & depth <= 324 ~ "deep",
        TRUE ~ "out_of_range"
      )
    )

  # ---- Assign condition categories ----
  bri_scores.2 <- bri_scores.1 %>%
    select(sample_id, stationid, sampledate, replicate, depth, depth_zone,
           shallow_bri_score, mid_bri_score, deep_bri_score, bri_score) %>%
    mutate(
      bri_cond = case_when(
        bri_score < 25  ~ "Reference",
        bri_score >= 25 & bri_score < 34 ~ "Marginal Deviation",
        bri_score >= 34 & bri_score < 44 ~ "Biodiversity Loss",
        bri_score >= 44 & bri_score < 72 ~ "Function Loss",
        bri_score >= 72 ~ "Defaunation"
      ),
      bri_class = case_when(
        bri_score < 25  ~ 1,
        bri_score >= 25 & bri_score < 34 ~ 2,
        bri_score >= 34 & bri_score < 44 ~ 3,
        bri_score >= 44 & bri_score < 72 ~ 4,
        bri_score >= 72 ~ 5
      )
    ) %>%
    ungroup() %>%
    bind_rows(defaunated, too.deep)

  # ---- Attach station info and add usage notes ----
  bri_w_station_info <- bri_scores.2 %>%
    left_join(
      station_info %>% select(stationid, latitude, longitude),
      by = "stationid"
    ) %>%
    mutate(
      note = case_when(
        sample_id %in% defaunated.samps ~ "No fauna in sample",
        latitude > 34.45  ~ "Caution - Sample outside of the geographic range of the BRI",
        latitude < 32.52  ~ "Caution - Sample outside of the geographic range of the BRI",
        depth > 200 & depth <= 324 ~ "Caution - Sample deeper than Bight recommendations but within BRI calibration",
        depth > 324 ~ "Do Not Use - Sample deeper than BRI calibration",
        TRUE ~ "None"
      )
    ) %>%
    relocate(depth, latitude, longitude, .after = replicate)

  writelog(
    '#### Offshore BRI Step 6 - Final BRI Scores\n',
    logfile = logfile,
    data = bri_w_station_info %>% head(25),
    verbose = verbose
  )
  create_download_link(data = bri_w_station_info, logfile = logfile, filename = 'BRI_offshore-final_scores.csv', linktext = 'Download final offshore BRI scores', verbose = verbose)

  writelog('\n### END: Offshore BRI (Edition 14) function.\n', logfile = logfile, verbose = verbose)

  return(bri_w_station_info)
}
