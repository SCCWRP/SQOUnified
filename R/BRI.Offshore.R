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
#'   In addition to the scores, the function returns the intermediate tables it builds along the way
#'   (taxa matched to p-codes, unmatched taxa with "did you mean" suggestions, and every taxon joined
#'   to its p-code by sample) so users can review exactly how their data were processed.
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
#'     \item Samples deeper than 324m are flagged "BRI not Applicable"
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
#'   index code sub-committee. Version 3.0 calculates BRI scores for overlapping depth zones by
#'   averaging the scores from each adjacent depth zone (rather than using overlap p-values), retains
#'   every in-range sample even when it has no taxa with a p-code in its depth zone, and flags samples
#'   whose organisms have no p-code with a caution note.
#'
#'   Geographic warnings are issued for samples north of Point Conception (lat > 34.45) or
#'   south of the US-Mexico border (lat < 32.52), as the BRI is nominally applicable only
#'   to the Southern California Bight. Samples whose taxon names cannot be matched to a p-code are
#'   returned in \code{taxa_without_pcode}; when the optional \code{fuzzyjoin} package is installed,
#'   a \code{did_you_mean_taxon} column suggests close matches (e.g. for likely misspellings).
#'
#' @param BenthicData a data frame containing infauna abundance and station information with the
#'    following columns (column names are case-insensitive):
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
#'    \strong{\code{Abundance}} - number of individuals counted for the specified taxon;
#'
#'    \strong{\code{Depth}} - station depth in meters;
#'
#'    \strong{\code{Latitude}} - station latitude in decimal degrees;
#'
#'    \strong{\code{Longitude}} - station longitude in decimal degrees (negative for west).
#'
#' @param output_format Quoted string controlling the shape of the \code{bri_scores} element of the
#'    returned list. \code{"wide"} (default) returns one row per sample with the score, condition
#'    category, class, and notes. \code{"long"} pivots to a tidy \code{index}/\code{score}/\code{category}/
#'    \code{category_score} layout (used by \code{\link{SQOUnified}}).
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @return
#'   A named \code{list} with the scores plus the intermediate tables users can review:
#'   \itemize{
#'     \item \code{bri_scores} - the BRI score, condition category, class, and usage notes for each
#'       sample. One row per sample when \code{output_format = "wide"}; a tidy index/score/category/
#'       category_score layout when \code{output_format = "long"}.
#'     \item \code{taxa_with_pcode} - distinct submitted taxa that matched a p-code.
#'     \item \code{taxa_without_pcode} - submitted taxa with no p-code, with fuzzy-matched
#'       \code{did_you_mean_taxon} suggestions when the \code{fuzzyjoin} package is installed.
#'     \item \code{all_taxa_by_sample} - every taxon joined to its p-code and depth-zone tolerance
#'       values, by sample.
#'   }
#'
#' @usage
#' BRI.Offshore(BenthicData, output_format = 'wide')
#'
#' @examples
#' data(offshore_bri_sampledata)
#' result <- BRI.Offshore(offshore_bri_sampledata)
#' result$bri_scores          # the BRI scores by station and replicate
#' result$taxa_without_pcode  # submitted taxa that could not be matched to a p-code
#'
#' @import dplyr
#' @importFrom tidyr drop_na pivot_longer
#' @importFrom lubridate as_date
#'
#' @export
BRI.Offshore <- function(BenthicData,
                         output_format = 'wide',
                         logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'BRI_Offshore_log.Rmd'),
                         verbose = FALSE,
                         knitlog = FALSE)
{

  # Backwards compatibility: old two-argument call BRI.Offshore(benthic_data, station_data)
  if (is.data.frame(output_format)) {
    shared_cols <- intersect(tolower(names(BenthicData)), tolower(names(output_format)))
    names(BenthicData)    <- tolower(names(BenthicData))
    names(output_format)  <- tolower(names(output_format))
    BenthicData  <- dplyr::left_join(BenthicData, output_format, by = shared_cols)
    output_format <- 'wide'
  }

  if (!output_format %in% c('wide','long')) stop('Invalid output format - ', output_format, ' - Acceptable values are "wide" or "long"')


  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n## BEGIN: Offshore BRI (Edition 14) function.\n', logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'offshore.bri.input.RData'
  if (verbose) {
    save(BenthicData, file = file.path(dirname(logfile), rawinput.filename))
  }

  writelog(
    '\nOffshore BRI calculation based on Smith et al. (2001), refined by the Bight 2023 index code sub-committee.',
    logfile = logfile,
    verbose = verbose
  )

  # Reference data (ed.14.ptaxa, pcodes) are exported package datasets (see ?ed.14.ptaxa, ?pcodes).
  # ed.14.ptaxa maps each taxon to a p_code; pcodes maps each p_code to its depth-zone tolerance values.

  # ---- Standardize column names to lowercase ----
  names(BenthicData) <- tolower(names(BenthicData))

  # ---- Prep the data ----
  infauna <- BenthicData %>%
    mutate(
      stationid = as.character(stationid),
      abundance = as.numeric(abundance),
      depth     = as.numeric(depth),
      sampledate = as_date(sampledate),
      sample_id = paste(stationid, sampledate, replicate, sep = "_")
    )

  taxa_to_calc <- infauna %>%
    select(sample_id, stationid, replicate, sampledate, taxon, abundance,
           depth, latitude, longitude) %>%
    arrange(sample_id, desc(abundance))

  writelog(
    '### Offshore BRI Step 1 - Taxa joined to station info\n',
    logfile = logfile,
    data = taxa_to_calc %>% head(25),
    verbose = verbose
  )
  create_download_link(data = taxa_to_calc, logfile = logfile, filename = 'BRI_offshore-step1-taxa_to_analyze.csv', linktext = 'Download taxa to be analyzed', verbose = verbose)


  # ---- Match taxa to p-codes ----
  # The pre-averaged shallow_mid / mid_deep transition-zone p-values are dropped; overlap-zone
  # scores are instead derived by averaging the two adjacent zone scores (see Step 5).
  all.4.bri <- taxa_to_calc %>%
    left_join(ed.14.ptaxa, by = "taxon") %>%
    left_join(pcodes, by = "p_code") %>%
    select(-shallow_mid, -mid_deep)

  # Taxa with assigned p-codes
  with.pcode <- all.4.bri %>%
    distinct(taxon, p_code) %>%
    drop_na(p_code) %>%
    arrange(taxon)

  writelog(
    '### Offshore BRI Step 2 - Taxa with assigned p-codes\n',
    logfile = logfile,
    data = with.pcode %>% head(25),
    verbose = verbose
  )
  create_download_link(data = with.pcode, logfile = logfile, filename = 'BRI_offshore-step2-taxa_with_pcode.csv', linktext = 'Download taxa with p-codes', verbose = verbose)

  # Taxa without p-codes. When fuzzyjoin (Suggests) is installed, attach close-match suggestions so
  # users can spot likely misspellings; otherwise return the submitted taxa with empty suggestions.
  no.pcode <- all.4.bri %>%
    distinct(taxon, p_code) %>%
    filter(is.na(p_code)) %>%
    filter(taxon != "NoOrganismsPresent") %>%
    select(-p_code) %>%
    arrange(taxon)

  if (requireNamespace("fuzzyjoin", quietly = TRUE)) {
    no.pcode <- no.pcode %>%
      fuzzyjoin::stringdist_left_join(ed.14.ptaxa, by = "taxon", max_dist = 2) %>%
      select(submitted_taxon = taxon.x, did_you_mean_taxon = taxon.y, p_code)
  } else {
    no.pcode <- no.pcode %>%
      mutate(did_you_mean_taxon = NA_character_, p_code = NA_character_) %>%
      select(submitted_taxon = taxon, did_you_mean_taxon, p_code)
  }

  writelog(
    '### Offshore BRI Step 3 - Taxa without p-codes (with did-you-mean suggestions)\n',
    logfile = logfile,
    data = no.pcode %>% head(25),
    verbose = verbose
  )
  create_download_link(data = no.pcode, logfile = logfile, filename = 'BRI_offshore-step3-taxa_without_pcode.csv', linktext = 'Download taxa without p-codes', verbose = verbose)

  writelog(
    '### Offshore BRI Step 4 - All taxa with p-codes by sample\n',
    logfile = logfile,
    data = all.4.bri %>% head(25),
    verbose = verbose
  )
  create_download_link(data = all.4.bri, logfile = logfile, filename = 'BRI_offshore-step4-all_taxa_by_sample.csv', linktext = 'Download all taxa by sample', verbose = verbose)


  # ---- Identify special-case samples ----

  # Defaunated samples (no organisms present) -> worst condition (class 5). depth is carried so the
  # downstream notes and final output report it.
  defaunated <- all.4.bri %>%
    filter(taxon == "NoOrganismsPresent") %>%
    distinct(sample_id, stationid, sampledate, replicate, depth) %>%
    mutate(
      tot_bri_abun_dz = 0,
      bri_score = NaN,
      bri_cond = "Defaunation",
      bri_class = 5
    )

  defaunated.samps <- unique(defaunated$sample_id)

  # Samples too deep for the BRI (>324m). anti_join keeps a sample that is both defaunated and too
  # deep classified as defaunated (the more informative label).
  too.deep <- all.4.bri %>%
    filter(depth > 324) %>%
    distinct(sample_id, stationid, sampledate, replicate, depth) %>%
    anti_join(defaunated, by = "sample_id") %>%
    mutate(
      tot_bri_abun_dz = NaN,
      bri_score = NaN,
      bri_cond = "BRI not Applicable",
      bri_class = NaN
    )

  too.deep.samps <- unique(too.deep$sample_id)

  # ---- Calculate depth-zone-specific BRI scores ----

  # Helper: drop defaunated / too-deep samples from the scoring calculation
  flag_exclusions <- function(df) {
    df %>%
      filter(!(sample_id %in% too.deep.samps), !(sample_id %in% defaunated.samps))
  }

  # Helper: calculate the BRI score and total p-coded abundance for a given depth-zone p-value column
  calc_zone_bri <- function(df, zone_col) {
    df %>%
      flag_exclusions() %>%
      filter(!is.na(.data[[zone_col]])) %>%
      mutate(
        cube_abun = abundance^(1/3),
        tol_score = .data[[zone_col]] * cube_abun
      ) %>%
      group_by(sample_id, stationid, sampledate, replicate, depth) %>%
      summarise(
        tot_bri_abun = sum(abundance),
        zone_bri_score = sum(tol_score) / sum(cube_abun),
        .groups = "drop"
      )
  }

  # Skeleton of every in-range, non-defaunated sample. Left-joining the zone scores onto this
  # guarantees a sample is never dropped just because it has no taxa with a p-code in a given zone.
  sample_skeleton <- all.4.bri %>%
    flag_exclusions() %>%
    distinct(sample_id, stationid, sampledate, replicate, depth)

  bri_scores.shallow <- calc_zone_bri(all.4.bri, "shallow") %>%
    rename(shallow_bri_score = zone_bri_score, tot_abun_shallow = tot_bri_abun)
  bri_scores.mid <- calc_zone_bri(all.4.bri, "mid") %>%
    rename(mid_bri_score = zone_bri_score, tot_abun_mid = tot_bri_abun)
  bri_scores.deep <- calc_zone_bri(all.4.bri, "deep") %>%
    rename(deep_bri_score = zone_bri_score, tot_abun_deep = tot_bri_abun)

  writelog(
    '\n### Offshore BRI Step 5 - Depth-zone BRI scores calculated (shallow, mid, deep)\n',
    logfile = logfile,
    verbose = verbose
  )

  # ---- Combine depth-zone scores and resolve the score for each sample's depth zone ----
  join_cols <- c("sample_id", "stationid", "sampledate", "replicate", "depth")

  zone_scores <- sample_skeleton %>%
    left_join(bri_scores.shallow, by = join_cols) %>%
    left_join(bri_scores.mid, by = join_cols) %>%
    left_join(bri_scores.deep, by = join_cols) %>%
    mutate(
      depth_zone = case_when(
        depth < 25 ~ "shallow",
        depth >= 25 & depth <= 35 ~ "shallow_mid",
        depth > 35 & depth < 110 ~ "mid",
        depth >= 110 & depth <= 130 ~ "mid_deep",
        depth > 130 & depth <= 324 ~ "deep",
        TRUE ~ "out_of_range"
      ),
      bri_score = case_when(
        depth < 25 ~ shallow_bri_score,
        depth >= 25 & depth <= 35 ~ (shallow_bri_score + mid_bri_score) / 2,
        depth > 35 & depth < 110 ~ mid_bri_score,
        depth >= 110 & depth <= 130 ~ (mid_bri_score + deep_bri_score) / 2,
        depth > 130 & depth <= 324 ~ deep_bri_score
      ),
      # Total p-coded abundance for the sample's depth zone (overlap zones average the two, rounded
      # up to keep the reported count an integer). Drives the "no organisms with P-code" caution.
      tot_bri_abun_dz = case_when(
        depth < 25 ~ tot_abun_shallow,
        depth >= 25 & depth <= 35 ~ ceiling((tot_abun_shallow + tot_abun_mid) / 2),
        depth > 35 & depth < 110 ~ tot_abun_mid,
        depth >= 110 & depth <= 130 ~ ceiling((tot_abun_mid + tot_abun_deep) / 2),
        depth > 130 & depth <= 324 ~ tot_abun_deep
      )
    )

  # ---- Assign condition categories ----
  bri_scores.2 <- zone_scores %>%
    select(sample_id, stationid, sampledate, replicate, depth, depth_zone,
           shallow_bri_score, mid_bri_score, deep_bri_score, tot_bri_abun_dz, bri_score) %>%
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
    bind_rows(defaunated, too.deep)

  # ---- Attach station info and add usage notes ----
  bri_w_station_info <- bri_scores.2 %>%
    left_join(
      infauna %>% distinct(stationid, latitude, longitude),
      by = "stationid"
    ) %>%
    mutate(
      note = case_when(
        sample_id %in% defaunated.samps ~ "No fauna in sample",
        tot_bri_abun_dz == 0 ~ "Caution - No organisms with P-code in Sample",
        latitude > 34.45  ~ "Caution - Sample outside of the geographic range of the BRI",
        latitude < 32.52  ~ "Caution - Sample outside of the geographic range of the BRI",
        depth > 200 & depth <= 324 ~ "Caution - Sample deeper than Bight recommendations but within BRI calibration",
        depth > 324 ~ "Do Not Use - Sample deeper than BRI calibration",
        TRUE ~ "None"
      )
    ) %>%
    # tot_bri_abun_dz was only needed to derive the caution note above
    select(-tot_bri_abun_dz) %>%
    relocate(depth, latitude, longitude, .after = replicate)

  writelog(
    '### Offshore BRI Step 6 - Final BRI Scores\n',
    logfile = logfile,
    data = bri_w_station_info %>% head(25),
    verbose = verbose
  )
  create_download_link(data = bri_w_station_info, logfile = logfile, filename = 'BRI_offshore-final_scores.csv', linktext = 'Download final offshore BRI scores', verbose = verbose)

  writelog('\n## END: Offshore BRI (Edition 14) function.\n', logfile = logfile, verbose = verbose)


  # output_format controls the shape of the main scores table only
  if (output_format == 'long') {
    bri_scores_out <- bri_w_station_info %>%
      select(
        stationid,
        sampledate,
        replicate,
        depth,
        notes = note,
        OFFSHORE_BRI_score = bri_score,
        OFFSHORE_BRI_category = bri_cond,
        OFFSHORE_BRI_category_score = bri_class
      ) %>%
      pivot_longer(
        cols = -c(stationid, sampledate, replicate, depth, notes),
        names_to = c("index", ".value"),
        names_pattern = "^(.+?)_(score|category_score|category)$"
      ) %>%
      mutate(
        index = if_else(
          index == 'OFFSHORE_BRI',
          'Offshore BRI',
          index
        )
      )
  } else {
    bri_scores_out <- bri_w_station_info
  }

  # Knit this offshore BRI log to HTML on a direct call (no-op when knitlog = FALSE, e.g. when
  # SQOUnified() drives this and merges the log into its consolidated report instead).
  knit.log(logfile, verbose = verbose, knitlog = knitlog)

  # Return the scores along with the intermediate tables so users can review how the index
  # processed their data (e.g. result$taxa_without_pcode).
  return(
    list(
      bri_scores         = bri_scores_out,
      taxa_with_pcode    = with.pcode,
      taxa_without_pcode = no.pcode,
      all_taxa_by_sample = all.4.bri
    )
  )
}
