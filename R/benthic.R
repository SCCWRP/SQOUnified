# GENERIC SQO BENTHIC LINE OF EVIDENCE (Alt) --------------------------------------------------
#' Compute the SQO Benthic Line of Evidence (BLOE) integrated score (generic version).
#'
#' @description
#'   Wrapper function that calculates all four SQO benthic indices (BRI, IBI, RBI, RIVPACS) plus
#'   M-AMBI, integrates them into a BLOE score, and performs SCAMIT Edition 14 to SQO taxonomy conversion.
#'
#' @details
#'   This function performs the following steps:
#'   \enumerate{
#'     \item Converts submitted SCAMIT Edition 14 taxonomy back to SQO-compatible taxonomy via one-to-one
#'           name swaps, daughter-to-parent rollups, and complex many-to-one resolutions
#'     \item Identifies taxa not on the SQO look-up list (with optional fuzzy matching suggestions
#'           if the \code{fuzzyjoin} package is installed)
#'     \item Calculates each benthic index: BRI, IBI, RBI, RIVPACS, and M-AMBI
#'     \item Integrates the four traditional SQO indices (BRI, IBI, RBI, RIVPACS) into a BLOE score
#'           by taking the ceiling of the median condition category score
#'     \item Returns BLOE scores with M-AMBI results attached
#'   }
#'
#'   IMPORTANT: The indices contained in this calculator are only intended for Southern California
#'   embayments (i.e., not offshore).
#'
#'   The SQO taxa lookup list from the online Excel SQO calculator is used (per Bight 23 index code
#'   subcommittee decision). The exclude code notation is used to eliminate ambiguous taxa.
#'
#'   BLOE condition categories (based on median of four index category scores, rounded up):
#'   \itemize{
#'     \item Reference: BLOE score = 1
#'     \item Low Disturbance: BLOE score = 2
#'     \item Moderate Disturbance: BLOE score = 3
#'     \item High Disturbance: BLOE score = 4
#'   }
#'
#' @param BenthicData a data frame containing benthic data and station information with at minimum:
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the sampling location;
#'
#'    \strong{\code{replicate}} - a numeric identifying the replicate number;
#'
#'    \strong{\code{sampledate}} - the date of sample collection;
#'
#'    \strong{\code{taxon}} - name of the organism (SCAMIT Edition 14 naming preferred).
#'        Use \code{NoOrganismsPresent} with 0 abundance for empty samples;
#'
#'    \strong{\code{abundance}} - number of individuals counted;
#'
#'    \strong{\code{exclude}} - "Yes" or "No" indicating if the taxon name is ambiguous;
#'
#'    \strong{\code{latitude}} - latitude in decimal degrees;
#'
#'    \strong{\code{longitude}} - longitude in decimal degrees (negative for west);
#'
#'    \strong{\code{depth}} - station depth in meters;
#'
#'    \strong{\code{salinity}} - bottom water salinity in PSU (used for BLOE applicability notes and M-AMBI).
#'
#' @param EG_Ref_values Optional. A data frame with Ecological Group assignments for M-AMBI.
#'    If NULL (default), uses the standard EG list.
#' @param EG_Scheme Quoted string specifying the EG value column to use for M-AMBI. Default "Hybrid".
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @usage
#' SQO_BLOE.generic(BenthicData)
#'
#' @examples
#' \dontrun{
#'   SQO_BLOE.generic(my_benthic_data)
#' }
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stringr str_flatten_comma str_detect str_flatten
#' @importFrom fuzzyjoin stringdist_left_join
#'
#' @export
benthic.sqo <- function(BenthicData,
                              EG_Ref_values = NULL,
                              EG_Scheme = "Hybrid",
                              logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'BLOE_generic_log.Rmd'),
                              verbose = FALSE,
                              knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n# BEGIN: Generic SQO BLOE function.\n', logfile = logfile, verbose = verbose)

  # Reference data (SoCal.SQO.ed14.link, ed14.rollups, ed.14.complex, xl_tool.SoCalLUList)
  # is available via R/sysdata.rda
  # Individual index functions (alt.SQO.BRI.generic, alt.IBI.generic, alt.RBI.generic,
  # alt.RIVPACS.generic, alt.MAMBI.generic) are available in the package namespace

  # Standardize the submitted data and retrofit the modern (Ed14) taxonomy back to
  # SQO-compatible names. benthicdata_prep() handles column standardization, the
  # Ed14 -> SQO taxonomy retrofitting, and identification of taxa not on the SQO
  # look-up list. Its taxonomy-change and unmatched-taxa logs are written within.
  prepped <- benthicdata_prep(BenthicData, retrofit = TRUE, logfile = logfile, verbose = verbose)

  # Retrofitted data is what the four traditional SQO indices are calculated on
  BenthicData.3 <- prepped$benthic_data

  # Taxa not on the SQO look-up list (returned to the caller for review)
  unmatched_taxa <- prepped$unmatched_taxa

  # M-AMBI uses modern (ed14) taxonomy, so it gets the standardized-but-not-retrofitted data
  BenthicData <- prepped$standardized


  #### Run each of the individual benthic indices ####

  writelog('\n## BLOE Step 3 - Calculating individual indices\n', logfile = logfile, verbose = verbose)

  # Data has already been retrofitted by benthicdata_prep above, so the indices
  # are told not to retrofit again (retrofit_taxonomy = FALSE).
  bri.scores.x <- BRI(BenthicData.3, retrofit_taxonomy = FALSE, logfile = logfile, verbose = verbose)
  writelog('\nSQO BRI Complete\n', logfile = logfile, verbose = verbose)

  ibi.scores.x <- IBI(BenthicData.3, retrofit_taxonomy = FALSE, logfile = logfile, verbose = verbose)
  writelog('\nSQO IBI Complete\n', logfile = logfile, verbose = verbose)

  rbi.scores.x <- RBI(BenthicData.3, retrofit_taxonomy = FALSE, logfile = logfile, verbose = verbose)
  writelog('\nSQO RBI Complete\n', logfile = logfile, verbose = verbose)

  rivpacs.scores.x <- RIVPACS(BenthicData.3, retrofit_taxonomy = FALSE, logfile = logfile, verbose = verbose)
  writelog('\nSQO RIVPACS Complete\n', logfile = logfile, verbose = verbose)

  # Note: MAMBI uses modern (ed14) taxonomy, not the SQO-retrofitted taxonomy
  mambi.scores.x <- MAMBI(BenthicData, EG_Ref_values = EG_Ref_values, EG_Scheme = EG_Scheme,
                                       logfile = logfile, verbose = verbose)
  writelog('\nSQO M-AMBI Complete\n', logfile = logfile, verbose = verbose)


  #### Integrate individual index scores into BLOE ####

  all.sqo.scores.x <- bind_rows(bri.scores.x, ibi.scores.x, rbi.scores.x, rivpacs.scores.x) %>%
    arrange(., stationid, sampledate, replicate, index)

  writelog(
    '## BLOE Step 4 - All benthic index scores\n',
    logfile = logfile,
    data = all.sqo.scores.x %>% head(25),
    verbose = verbose
  )
  create_download_link(data = all.sqo.scores.x, logfile = logfile, filename = 'BLOE-step3-all_index_scores.csv', linktext = 'Download all index scores', verbose = verbose)


  BLOE.scores.x <- all.sqo.scores.x %>%
    group_by(stationid, sampledate, replicate) %>%
    # Integrate by calculating the ceiling of the median condition category score
    mutate(BLOE_score = ceiling(median(condition_category_score, na.rm = TRUE)),
           Note = str_flatten(note, collapse = ",")) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(-score, -condition_category, -note),
                names_from = index, values_from = condition_category_score) %>%
    mutate(notes = case_when(is.na(salinity) & is.na(Note) ~ "Salinity value unknown-confirm salinity >=27 PSU for BLOE to be accurate.",
                             is.na(salinity) & !is.na(Note) ~ paste("Salinity value unknown-confirm salinity >=27 PSU for BLOE to be accurate.", Note, sep = "; "),
                             salinity < 27 & is.na(Note) ~ "Caution, salinity value <27 PSU - BLOE may not be accurate. Consider M-AMBI",
                             salinity < 27 & !is.na(Note) ~ paste("Caution, salinity value <27 PSU - BLOE may not be accurate. Consider M-AMBI", Note, sep = "; "),
                             salinity >= 27 & !is.na(Note) ~ Note,
                             salinity >= 27 & is.na(Note) ~ "None",
                             TRUE ~ "None"),
           BLOE_category = case_when(BLOE_score == 1 ~ "Reference",
                                     BLOE_score == 2 ~ "Low Disturbance",
                                     BLOE_score == 3 ~ "Moderate Disturbance",
                                     BLOE_score == 4 ~ "High Disturbance")
    ) %>%
    relocate(stationid, sampledate, replicate, BLOE_score, BLOE_category, BRI_cond = BRI, IBI_cond = IBI, RBI_cond = RBI, Rivpacs_cond = Rivpacs) %>%
    select(-Note)

  # Attach M-AMBI scores to BLOE output
  mambi.scores.sqoformat <- mambi.scores.x %>%
    mutate(MAMBI_cond = case_when(SQO_mambi_condition == "Reference" ~ 1,
                                  SQO_mambi_condition == "Low Disturbance" ~ 2,
                                  SQO_mambi_condition == "Moderate Disturbance" ~ 3,
                                  SQO_mambi_condition == "High Disturbance" ~ 4,
                                  TRUE ~ NA), .after = SQO_mambi_condition)

  BLOE.scores.w.MAMBI <- BLOE.scores.x %>%
    left_join(., select(mambi.scores.sqoformat, stationid, sampledate, replicate, SQO_mambi_condition, MAMBI_cond, note),
              by = c("stationid", "sampledate", "replicate")) %>%
    mutate(notes = case_when(note == "None" ~ notes,
                             is.na(note) ~ notes,
                             TRUE ~ paste(notes, note, sep = "; "))) %>%
    select(-note) %>%
    rename(MAMBI_SQO_Cat = SQO_mambi_condition) %>%
    relocate(MAMBI_SQO_Cat, MAMBI_cond, .after = Rivpacs_cond)


  writelog(
    '## BLOE Final - Integrated BLOE scores with M-AMBI\n',
    logfile = logfile,
    data = BLOE.scores.w.MAMBI %>% head(25),
    verbose = verbose
  )
  create_download_link(data = BLOE.scores.w.MAMBI, logfile = logfile, filename = 'BLOE-final_scores.csv', linktext = 'Download BLOE scores', verbose = verbose)

  writelog('\n# END: Generic SQO BLOE function.\n', logfile = logfile, verbose = verbose)

  # Return list with all results for inspection
  results <- list("sqo_bloe" = BLOE.scores.x,
                  "sqo_bloe_w_mambi" = BLOE.scores.w.MAMBI,
                  "sqo_bri" = bri.scores.x,
                  "ibi" = ibi.scores.x,
                  "rbi" = rbi.scores.x,
                  "rivpacs" = rivpacs.scores.x,
                  "mambi" = mambi.scores.sqoformat,
                  "sqo_potential_mismatches" = unmatched_taxa)
  
  sqo_bloe <- results$sqo_bloe
  bri <- results$sqo_bri
  ibi <- results$ibi
  rbi <- results$rbi
  rivpacs <- results$rivpacs



  # Build a comments column with depth, salinity, and stratum for human review
  sqo_bloe <- sqo_bloe %>%
    select(
      stationid, 
      sampledate, 
      replicate, 
      depth,
      salinity,
      stratum, 
      # same thing for the integrated benthic assessment score
      BLOE_score,
      BLOE_category_score = BLOE_score, 
      
      BLOE_category,
      notes
    ) 

  bri <- bri %>% 
    select(
      stationid, 
      sampledate, 
      replicate, 
      BRI_score = score,
      BRI_category_score = condition_category_score,
      BRI_category = condition_category
    )
    
  ibi <- ibi %>% 
    select(
      stationid, 
      sampledate, 
      replicate, 
      IBI_score = score,
      IBI_category_score = condition_category_score,
      IBI_category = condition_category
    )
    
  rbi <- rbi %>% 
    select(
      stationid, 
      sampledate, 
      replicate, 
      RBI_score = score,
      RBI_category_score = condition_category_score,
      RBI_category = condition_category
    )
    
  rivpacs <- rivpacs %>% 
    select(
      stationid, 
      sampledate, 
      replicate, 
      RIVPACS_score = score,
      RIVPACS_category_score = condition_category_score,
      RIVPACS_category = condition_category
    )
    



  joined <- sqo_bloe %>% 
    left_join(bri, by = c('stationid', 'sampledate', 'replicate')) %>%
    left_join(ibi, by = c('stationid', 'sampledate', 'replicate')) %>%
    left_join(rbi, by = c('stationid', 'sampledate', 'replicate')) %>%
    left_join(rivpacs, by = c('stationid', 'sampledate', 'replicate')) 

  # Pivot longer: BLOE, BRI, IBI, RBI, RIVPACS, MAMBI each become a row
  # with their _score, _category_score, and _category values in dedicated columns
  allsqo_long <- joined %>%
    pivot_longer(
      cols = -c(stationid, sampledate, replicate, depth, salinity, stratum, notes),
      names_to = c("index", ".value"),
      names_pattern = "^(.+?)_(score|category_score|category)$"
    ) %>%
    mutate(
      index = case_when(
        index == "BLOE" ~ "Integrated Benthic LOE SQO Assessment Score",
        TRUE ~ index
      )
    )

  

  results$all_benthic_sqo_scores_long <- allsqo_long

  return(results)
}


# BENTHIC DATA PREP ---------------------------------------------------------------------------
#' Standardize and retrofit submitted benthic data for SQO index calculation.
#'
#' @description
#'   Utility function that prepares submitted benthic data for the traditional SQO benthic indices
#'   (BRI, IBI, RBI, RIVPACS). It standardizes column names and types, normalizes missing sentinel
#'   values for salinity and depth, and---optionally---retrofits modern SCAMIT Edition 14 taxonomy
#'   back to the SQO-compatible names that the indices recognize. It also identifies taxa that are
#'   not on the SQO look-up list.
#'
#' @details
#'   When \code{retrofit = TRUE} the following taxonomy conversion steps are applied:
#'   \enumerate{
#'     \item One-to-one name changes (swap new names back to old versions)
#'     \item Daughter-to-parent rollups (combine daughter taxa to the higher taxonomic level the
#'           indices recognize)
#'     \item Complex many-to-one resolutions (a single ed14 taxon that mapped to multiple SQO taxa)
#'   }
#'   When \code{retrofit = FALSE} only the standardization steps are performed and the submitted
#'   taxonomy is passed through unchanged.
#'
#'   This function is called at the beginning of \code{\link{BRI}}, \code{\link{IBI}},
#'   \code{\link{RBI}}, and \code{\link{RIVPACS}}. It is intentionally NOT used by \code{MAMBI}
#'   (which uses modern ed14 taxonomy) or by \code{BRI.Offshore}.
#'
#'   Reference data (\code{SoCal.SQO.ed14.link}, \code{ed14.rollups}, \code{ed.14.complex},
#'   \code{xl_tool.SoCalLUList}) is available via \code{R/sysdata.rda}.
#'
#' @param BenthicData a data frame containing benthic data and station information with at minimum:
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the sampling location;
#'
#'    \strong{\code{replicate}} - a numeric identifying the replicate number;
#'
#'    \strong{\code{sampledate}} - the date of sample collection;
#'
#'    \strong{\code{taxon}} - name of the organism (SCAMIT Edition 14 naming preferred).
#'        Use \code{NoOrganismsPresent} with 0 abundance for empty samples;
#'
#'    \strong{\code{abundance}} - number of individuals counted;
#'
#'    \strong{\code{exclude}} - "Yes" or "No" indicating if the taxon name is ambiguous.
#'
#'    \code{salinity} and \code{depth} are normalized (sentinel values -88 and -99 set to NA) when present.
#'
#' @param retrofit Logical. If TRUE (default), modern Ed14 taxonomy is retrofitted back to
#'    SQO-compatible names. If FALSE, the submitted taxonomy is passed through unchanged.
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{benthic_data} - the prepared (and, when \code{retrofit = TRUE}, retrofitted) data frame for index calculation;
#'     \item \code{standardized} - the standardized data frame prior to taxonomy retrofitting;
#'     \item \code{taxa_changes} - a data frame detailing the taxonomy changes made (NULL when \code{retrofit = FALSE});
#'     \item \code{unmatched_taxa} - taxa not found on the SQO look-up list (with fuzzy-match suggestions if \code{fuzzyjoin} is installed).
#'   }
#'
#' @usage
#' benthicdata_prep(BenthicData, retrofit = TRUE)
#'
#' @examples
#' \dontrun{
#'   benthicdata_prep(my_benthic_data)
#' }
#'
#' @import dplyr
#' @importFrom stringr str_flatten_comma str_detect
#' @importFrom fuzzyjoin stringdist_left_join
#'
#' @export
benthicdata_prep <- function(BenthicData,
                             retrofit = TRUE,
                             logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'benthicdata_prep_log.Rmd'),
                             verbose = FALSE,
                             knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n## BEGIN: benthicdata_prep function.\n', logfile = logfile, verbose = verbose)

  # Reference data (SoCal.SQO.ed14.link, ed14.rollups, ed.14.complex, xl_tool.SoCalLUList)
  # is available via R/sysdata.rda

  # Standardize column names and ensure correct types
  names(BenthicData) <- tolower(names(BenthicData))

  BenthicData <- BenthicData %>%
    mutate(stationid = as.character(stationid),
           abundance = as.numeric(abundance))

  # Standardize missing sentinel values for salinity and depth (only if present)
  if ("salinity" %in% names(BenthicData)) {
    BenthicData <- BenthicData %>%
      mutate(salinity = ifelse(salinity %in% c(-88, -99), NA, salinity))
  }
  if ("depth" %in% names(BenthicData)) {
    BenthicData <- BenthicData %>%
      mutate(depth = ifelse(depth %in% c(-88, -99), NA, depth))
  }

  # Keep a copy of the standardized (non-retrofitted) data for consumers that
  # expect modern (ed14) taxonomy (e.g., M-AMBI).
  standardized <- BenthicData

  if (retrofit) {

    #### Retrofitting modern taxonomy back to SQO-compatible names ####

    # Step 1: One-to-one name changes (swap new names back to old versions)
    one.to.one <- SoCal.SQO.ed14.link %>%
      select(original_sqo_taxon, ed_14_taxon, change, type) %>%
      filter(type %in% c("one-to-one", "convention", "worms", "removed"))

    BenthicData.0 <- BenthicData %>%
      left_join(., one.to.one, by = c("taxon" = "ed_14_taxon")) %>%
      relocate(original_sqo_taxon, .before = taxon) %>%
      mutate(taxon.2 = if_else(is.na(original_sqo_taxon), taxon, original_sqo_taxon),
             change_type = if_else(is.na(type), "", type),
             taxa_changed = if_else(is.na(original_sqo_taxon), "", taxon)) %>%
      select(-taxon, -type, -change, -original_sqo_taxon)


    # Step 2: Combine daughter taxa to the higher taxonomic level the indices recognize
    BenthicData.1 <- BenthicData.0 %>%
      relocate(abundance, taxon.2, .before = stationid) %>%
      left_join(., ed14.rollups, by = c("taxon.2" = "ed.14_daughters")) %>%
      mutate(taxon.3 = if_else(is.na(sqo.name), taxon.2, sqo.name),
             # rolled_up = if_else(is.na(sqo.name), "no", "yes"), .before = taxon.2) %>%
             rolled_up=case_when(is.na(sqo.name)~"no", taxon.3==taxon.2~"no", TRUE~"yes"), .before = taxon.2) %>%
      select(-sqo.name) %>%
      group_by(stationid, replicate, sampledate, taxon.3) %>%
      mutate(abundance.2 = sum(abundance),
             # When rolling up taxa we may combine taxa w/ different exclude codes
             # If any of the combined taxa has "No" (i.e., it is unique), give that precedence
             all.excludes = str_flatten_comma(exclude),
             exclude.prob.flag = case_when(str_detect(all.excludes, "No, Yes") ~ 1,
                                           str_detect(all.excludes, "Yes, No") ~ 1,
                                           TRUE ~ 0),
             exclude.2 = case_when(rolled_up == "yes" & exclude.prob.flag == 1 ~ "No",
                                   TRUE ~ exclude),
             .before = abundance) %>%
      ungroup() %>%
      mutate(change_type.2 = if_else(rolled_up == "yes", paste("rolled up to ", join_level, sep = ""), change_type),
             taxa_changed.2 = if_else(rolled_up == "yes", taxon.2, taxa_changed),
             .before = stationid) %>%
      select(-all.excludes, -exclude.prob.flag, -exclude) %>%
      rename(exclude = exclude.2)


    # Step 3: Complex changes where a single ed14 taxon was multiple taxa on the SQO LU list
    ed.14.complex.2 <- ed.14.complex %>%
      filter(Priority == "yes")

    BenthicData.2 <- BenthicData.1 %>%
      left_join(., ed.14.complex.2, by = c("taxon.3" = "Ed.14.Taxon")) %>%
      mutate(taxon.4 = case_when(Priority == "yes" ~ Original.SQO.Taxon,
                                 TRUE ~ taxon.3), .before = taxon.3) %>%
      mutate(change_type.3 = case_when(Priority == "yes" ~ "complex",
                                       # TRUE ~ change_type.2), .before = change_type.2)
                                       TRUE~change_type.2),
                                       taxa_changed.3=if_else(change_type.3=="complex", taxon.3, taxa_changed.2),
                                       .before = change_type.2)

    # Step 4: Interim output detailing all changes made to submitted data
    taxa.changes <- BenthicData.2 %>%
      group_by(stationid, replicate, sampledate, taxon.4, change_type.3) %>%
      mutate(
        taxa_changed.4 = str_flatten_comma(taxa_changed.3),
        .before = abundance) %>%
        ungroup() %>%
        mutate(
          taxon_submitted = case_when(
            change_type %in% c("one-to-one", "convention") ~ taxa_changed,
            TRUE~taxon.2
          ),
          .after = taxon.4
        ) %>%
      select(stationid, sampledate, replicate, taxon_used = taxon.4, taxon_submitted = taxon.2, changed_taxa = taxa_changed.4,
             abundance_used = abundance.2,
             abundance_submitted = abundance, type_of_change = change_type.3)

    writelog(
      '### BLOE Step 1 - Ed14 to SQO taxonomy changes\n',
      logfile = logfile,
      data = taxa.changes %>% head(25),
      verbose = verbose
    )
    create_download_link(data = taxa.changes, logfile = logfile, filename = 'BLOE-step1-taxonomy_changes.csv', linktext = 'Download taxonomy changes', verbose = verbose)


    # Step 5: Create final dataframe for index calculation
    BenthicData.3 <- BenthicData.2 %>%
      select(-abundance, -taxon.2, -taxon.3, -rolled_up, -change, -change_type, -change_type.2, -change_type.3,
             -taxa_changed, -taxa_changed.2, -join_level, -Original.SQO.Taxon, -type, -Priority) %>%
      distinct() %>%
      rename(taxon = taxon.4, abundance = abundance.2) %>%
      relocate(taxon, abundance, .after = sampledate)

  } else {

    writelog('\nTaxonomy retrofitting skipped (retrofit = FALSE); submitted taxonomy used as-is.\n', logfile = logfile, verbose = verbose)

    BenthicData.3 <- BenthicData
    taxa.changes <- NULL

  }

  # Identify unmatched taxa (not on SQO look-up list)
  unmatched_taxa <- BenthicData.3 %>%
    left_join(., xl_tool.SoCalLUList, by = c("taxon" = "TaxonName")) %>%
    select(exclude, stationid, replicate, taxon, Phylum, Class, Order, Family, ToleranceScore, IBISensitive, Mollusc, Crustacean) %>%
    mutate(unmatched_flag = case_when(exclude != "Yes" & is.na(Phylum) & is.na(ToleranceScore) & is.na(IBISensitive) & is.na(Mollusc) & is.na(Crustacean) ~ "unmatched",
                                     TRUE ~ "matching")) %>%
    filter(unmatched_flag == "unmatched") %>%
    select(stationid, replicate, taxon) %>%
    filter(taxon != "NoOrganismsPresent")

  # Fuzzy match suggestions if fuzzyjoin is available
  if (requireNamespace("fuzzyjoin", quietly = TRUE)) {
    unmatched_taxa <- fuzzyjoin::stringdist_left_join(
      unmatched_taxa,
      select(xl_tool.SoCalLUList, TaxonName),
      by = c("taxon" = "TaxonName"),
      max_dist = 2
    ) %>%
      rename(not_on_LU_list = taxon, did_you_mean = TaxonName)
  } else {
    unmatched_taxa <- unmatched_taxa %>%
      rename(not_on_LU_list = taxon)
  }

  writelog(
    '### BLOE Step 2 - Unmatched taxa\n',
    logfile = logfile,
    data = unmatched_taxa %>% head(25),
    verbose = verbose
  )
  create_download_link(data = unmatched_taxa, logfile = logfile, filename = 'BLOE-step2-unmatched_taxa.csv', linktext = 'Download unmatched taxa', verbose = verbose)

  writelog('\n## END: benthicdata_prep function.\n', logfile = logfile, verbose = verbose)

  return(
    list(
      benthic_data   = BenthicData.3,
      standardized   = standardized,
      taxa_changes   = taxa.changes,
      unmatched_taxa = unmatched_taxa
    )
  )
}
