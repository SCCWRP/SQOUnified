#' Compute the SQO index scores.
#' Based on the Technical Manual - June 2021 edition
#' https://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf
#'
#' When this function is called, it will calculate scores based on whatever data is passed into it.
#' If Chem, Tox, and Benthic are provided, integrated scores will be calculated based on the 3 LOE's.
#'
#'
#' @param benthic a data frame. This data frame must contain the following
#'  information with these headings:
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
#'    \code{Taxon} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#'    \code{Abundance} - the number of each Species observed in a sample;
#'
#'    \code{Salinity} - the salinity observed at the location in PSU, ideally at time of sampling;
#'
#'    \code{Stratum} - The stratum under which the station falls (Bays, Estuaries, etc);
#'
#'    \code{Exclude} - Yes or No;
#'
#'
#' @param chem a dataframe with the following columns:
#'
#'    \code{StationID}
#'
#'    \code{AnalyteName}
#'
#'    \code{Result}
#'
#'    \code{RL}
#'
#'    \code{MDL}
#'
#'    \code{units} (optional) -  Metals should be in mg/dry kg (mg/kg dw) and all organic constituents should be in ug/dry kg (ug/kg dw).
#'
#'    \code{fieldrep} (optional) - data will be filtered to where fieldrep = 1
#'
#'    \code{labrep} (optional) - data will be filtered to where labrep = 1
#'
#'    \code{sampletypecode} (optional) - data will be filtered to where sampletypecode = Result (to avoid including data from QA/QC samples)
#'
#'
#'
#' @param toxresults a dataframe with the following columns: stationid, toxbatch, species, sampletypecode
#'    matrix, labrep, result. This data must also include the control samples
#'    (stationcode 0000, sampletypecode CNEG etc)
#'
#'    The input dataframe is structured as follows
#'
#'    \strong{\code{lab}} -  (optional) The laboratory which performed the test
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the location;
#'
#'    \strong{\code{toxbatch}} - the toxbatch id - used to join with the control sample
#'
#'    \strong{\code{species}} - The Genus and species of the animale that was tested
#'
#'    \strong{\code{sampletypecode}} - The sampletype used Grab, CNEG etc. Control samples must be included
#'
#'    \strong{\code{matrix}} - (optional) Whole Sediment, Sediment Water Interface, etc. Be sure to not include Reference Toxicants
#'
#'    \strong{\code{labrep}} - There should be 5 per station, species pair
#'
#'    \strong{\code{result}} - the percentage that survived the test, or had normal development
#'
#'
#' @examples
#' SQOUnified(benthic = benthic_sampledata, chem = chem_sampledata, tox = tox_sampledata)
#'
#' @importFrom dplyr case_when full_join select rename mutate arrange
#' @importFrom purrr map
#' @importFrom tidyr spread


#' @export
SQOUnified <- function(benthic = NULL, chem = NULL, tox = NULL, offshore_benthic = NULL, offshore_stations = NULL, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'log.Rmd' ), verbose = F, logtitle = 'Unified SQO Logs', knitlog = F) {

  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, title = logtitle)


  #load("data/site_assessment_criteria.RData")

  if (all(is.null(c(benthic,chem,tox,offshore_benthic, offshore_stations)))){
    stop("
      No data was provided.
      Please provide benthic, chemistry and toxicity data to get the integrated site assessments
    ")
  }
  # check the data coming in before anything
  checkdata(benthic, chem, tox, offshore_benthic, offshore_stations, logfile = logfile, verbose = verbose)



  # Compute ALL SQO scores
  writelog('Compute ALL SQO scores', logfile = logfile, verbose = verbose)

  # Each line of evidence writes its own self-contained log (with its intermediate tables) into a
  # subfolder. We record them here and, at the end, merge them into this one consolidated report
  # which is knit to a single self-contained HTML (see the "Assemble the consolidated report" block).
  report.sections <- list()

  # ---- Toxicity ----
  if (!is.null(tox)) {

    toxlogfile <- file.path( dirname(logfile), 'Toxicity', 'toxlog.Rmd' )
    toxlibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(toxlogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = toxlibs, title = 'Toxicity SQO Logs')

    tox <- tox.sqo(tox, logfile = toxlogfile, verbose = verbose, knitlog = FALSE) %>%
      mutate(LOE = 'Toxicity') %>%
      select(StationID, LOE, Index, Score, Category, `Category Score`)

    report.sections[[length(report.sections) + 1]] <- list(title = 'Toxicity', logfile = toxlogfile, subdir = 'Toxicity')

  } else {
    tox <- data.frame(
      StationID = c(),
      LOE = c(),
      Index = c(),
      Score = c(),
      Category = c(),
      `Category Score` = c()
    )
  }

  # ---- Chemistry ----
  if (!is.null(chem)) {

    chemlogfile <- file.path( dirname(logfile), 'Chemistry', 'chemlog.Rmd' )
    chemlibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(chemlogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = chemlibs, title = 'Chemistry SQO Logs')

    chem <- chem.sqo(chem, logfile = chemlogfile, verbose = verbose, knitlog = FALSE) %>%
      mutate(LOE = 'Chemistry') %>%
      select(StationID, LOE, Index, Score, Category, `Category Score`)

    report.sections[[length(report.sections) + 1]] <- list(title = 'Chemistry', logfile = chemlogfile, subdir = 'Chemistry')

  } else {
    chem = data.frame(
      StationID = c(),
      LOE = c(),
      Index = c(),
      Score = c(),
      Category = c(),
      `Category Score` = c()
    )
  }

  

  # ---- Benthic ----
  if (!is.null(benthic)) {

    benthiclogfile <- file.path( dirname(logfile), 'Benthic', 'benthiclog.Rmd' )
    benthiclibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(benthiclogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = benthiclibs, title = 'Benthic SQO Logs')

    # Work in progress here - probably have to extract all the dataframes from the list to extract their scores
    benthic <- benthic.sqo(benthic, logfile = benthiclogfile, verbose = verbose, knitlog = FALSE)$all_benthic_sqo_scores_long %>%
      mutate(LOE = 'Benthic') %>%
      rename(
        StationID = stationid,
        Replicate = replicate,
        Index = index,
        Score = score,
        Category = category,
        `Category Score` = category_score
      ) %>%
      select(StationID, Replicate,  LOE, Index, Score, Category, `Category Score`) %>%
      # David says only keep replicate 1.
      # 999 times out of 1000 there should only be one replicate
      filter(
        Replicate == 1
      ) %>%
      # The way i see it now is, Replicate becomes a useless field.
      # The other datatypes dont use it
      # SampleDate is not very useful either, not in the Bight program at least
      # The fact that the station is called B18 shows it was sampled in 2018, which is all we care about
      # besides, I queried the tbl_infaunalabundance_initial table,
      # there are not two sampledates for the same station, unlike SMC, which has permanent station names,
      # and sampledates are used to distinguish the station at different times
      # The Bight program however, which only samples every 5 years, puts the year of the sample in the stationid
      select(-c(Replicate))

    report.sections[[length(report.sections) + 1]] <- list(title = 'Benthic', logfile = benthiclogfile, subdir = 'Benthic')

  } else {
    benthic = data.frame(
      StationID = c(),
      LOE = c(),
      Index = c(),
      Score = c(),
      Category = c(),
      `Category Score` = c()
    )
  }


  # ---- Offshore BRI ----
  if ( (!is.null(offshore_benthic)) && (!is.null(offshore_stations)) ) {

    obrilogfile <- file.path( dirname(logfile), 'Offshore BRI', 'offshorebri_log.Rmd' )
    obrilibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(obrilogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = obrilibs, title = 'Offshore BRI SQO Logs')

    # BRI.Offshore takes a single combined data frame; join benthic and station data here.
    offshore_combined <- dplyr::left_join(
      offshore_benthic %>% dplyr::rename_with(tolower),
      offshore_stations %>% dplyr::rename_with(tolower),
      by = "stationid"
    )

    # BRI.Offshore returns a list (scores + interim tables); pull the long-format scores table
    offshore_bri <- BRI.Offshore(offshore_combined, output_format = 'long', logfile = obrilogfile, verbose = verbose, knitlog = FALSE)$bri_scores %>%
      # Keep only replicate 1, mirroring the bay/estuary benthic rule (David: keep replicate 1;
      # there is almost always only one). This gives each station a single Benthic row so the
      # spread() in the integrated assessment cannot hit duplicate StationID + LOE keys.
      filter(replicate == 1) %>%
      mutate(LOE = 'Benthic') %>%
      select(
          StationID = stationid,
          LOE,
          Index = index,
          Score = score,
          Category = category,
          `Category Score` = category_score
      )

    report.sections[[length(report.sections) + 1]] <- list(title = 'Offshore BRI', logfile = obrilogfile, subdir = 'Offshore BRI')

  } else {
    offshore_bri <- data.frame(
      StationID = c(),
      LOE = c(),
      Index = c(),
      Score = c(),
      Category = c(),
      `Category Score` = c()
    )
  }

  # For integrated site assessment, each station needs one Benthic, Chemistry, and Toxicity category.
  # Stations with bay/estuary benthic data use the BLOE category.
  # Stations with only offshore BRI use the offshore BRI category mapped to standard benthic categories.
  # Build the wide, one-row-per-station table of LOE categories joined to the site assessment
  # criteria. The Benthic category comes from either the bay/estuary BLOE (benthic.sqo) or the
  # offshore BRI (both carry LOE = 'Benthic'); a station needs all three LOEs to receive an
  # integrated assessment. This wide table is also displayed in the consolidated report below.
  integrated_wide <- bind_rows(benthic, chem, tox, offshore_bri) %>%
    filter(
      grepl("Assessment|Offshore", Index)
    ) %>%
    # Map offshore BRI categories to standard benthic LOE categories
    mutate(
      Category = case_when(
        Category == "Marginal Deviation" ~ "Low Disturbance",
        Category == "Biodiversity Loss"  ~ "Moderate Disturbance",
        Category == "Function Loss"      ~ "High Disturbance",
        Category == "Defaunation"        ~ "High Disturbance",
        TRUE ~ Category
      )
    ) %>%
    select(
      StationID, LOE, Category
    ) %>%
    group_by(
      StationID
    ) %>%
    spread(
      LOE, Category
    ) %>%
    # An entire line of evidence can be absent from a run (e.g. only benthic and/or offshore
    # data was supplied), in which case spread() never creates that column and the criteria
    # join below would error on the missing column. Guarantee all three category columns exist;
    # any station missing one or more LOEs then simply gets NA and so no integrated assessment.
    {
      missing_loes <- setdiff(c("Benthic", "Chemistry", "Toxicity"), names(.))
      if (length(missing_loes) > 0) .[missing_loes] <- NA_character_
      .
    } %>%
    left_join(
      site_assessment_criteria,
      by = c("Benthic","Chemistry","Toxicity")
    ) %>%
    ungroup()

  # Reshape the wide table into the long score format used for the returned data frame
  integrated <- integrated_wide %>%
    select(
      StationID, `Site Assessment`
    ) %>%
    mutate(
      LOE = "Integrated",
      Index = "Site Assessment",
      Score = NA_real_,
      `Category Score` = NA_real_
    ) %>%
    select(
      StationID,
      LOE,
      Index,
      Score,
      Category = `Site Assessment`,
      `Category Score`
    )

  out <- bind_rows(
    benthic, chem, tox, offshore_bri, integrated
  ) %>%
  arrange(
    StationID, LOE, Index
  )

  writelog('Done computing ALL SQO scores', logfile = logfile, verbose = verbose)

  # ---- Assemble the consolidated report ----
  # Merge each line-of-evidence log into this single report (in a sensible reading order) and knit
  # it to one self-contained HTML. Because each section's tables are read from the snapshots
  # writelog() saved, nothing is recomputed at knit time. The per-LOE intermediate CSVs remain in
  # their respective subfolders.
  if (length(report.sections) > 0) {
    # Sections appear in the consolidated report in the order their line-of-evidence blocks run
    # above (currently Toxicity, then Chemistry, then Benthic, then Offshore BRI). To change the
    # report order, reorder those blocks - there is no separate ordering list to keep in sync.
    for (sec in report.sections) {
      append_log_section(logfile, sec$logfile, subdir = sec$subdir, title = sec$title, verbose = verbose)
    }

    # ---- Integrated site assessment ----
    # The integration runs in this master function (not a per-LOE sub-log), so its code and result
    # tables are written directly into the consolidated report here, after the per-LOE sections.
    writelog('\n# Integrated Site Assessment\n', logfile = logfile, verbose = verbose)
    writelog(
      paste(
        'The overall site assessment integrates the Benthic, Chemistry, and Toxicity category for',
        'each station. The Benthic category comes from either the bay/estuary benthic line of',
        'evidence (benthic.sqo) or the offshore BRI (BRI.Offshore). A station must have all three',
        'lines of evidence to receive an integrated site assessment; stations missing any line of',
        'evidence appear with an NA site assessment.'
      ),
      logfile = logfile,
      verbose = verbose
    )

    writelog(
      '\n## Site assessment criteria\n\nLook-up table mapping each combination of Benthic, Chemistry, and Toxicity categories to an overall site assessment (CASQO Technical Manual, Bay et al. 2021).',
      logfile = logfile,
      data = site_assessment_criteria,
      verbose = verbose,
      pageLength = 15
    )

    writelog(
      '\n## Integrated categories by station\n',
      logfile = logfile,
      code = '
        # Reduce each line of evidence to its site-level category, one row per station. The Benthic
        # category comes from either benthic.sqo or BRI.Offshore (both carry LOE = "Benthic", with
        # the offshore BRI categories mapped onto the standard benthic categories first). Absent
        # lines of evidence are filled with NA so every station still appears.
        integrated_wide <- bind_rows(benthic, chem, tox, offshore_bri) %>%
          filter(grepl("Assessment|Offshore", Index)) %>%
          mutate(
            Category = case_when(
              Category == "Marginal Deviation" ~ "Low Disturbance",
              Category == "Biodiversity Loss"  ~ "Moderate Disturbance",
              Category == "Function Loss"      ~ "High Disturbance",
              Category == "Defaunation"        ~ "High Disturbance",
              TRUE ~ Category
            )
          ) %>%
          select(StationID, LOE, Category) %>%
          group_by(StationID) %>%
          spread(LOE, Category) %>%
          # Fill any line of evidence absent from this run with NA so the criteria join never
          # fails on a missing column (single-LOE / offshore-only runs).
          {
            missing_loes <- setdiff(c("Benthic", "Chemistry", "Toxicity"), names(.))
            if (length(missing_loes) > 0) .[missing_loes] <- NA_character_
            .
          } %>%
          left_join(site_assessment_criteria, by = c("Benthic", "Chemistry", "Toxicity")) %>%
          ungroup()

        # The joined "Site Assessment" column is the overall integrated category for each station
      ',
      data = integrated_wide %>% select(StationID, Benthic, Chemistry, Toxicity, `Site Assessment`),
      verbose = verbose
    )
    # Write the integrated assessment scores to a CSV in the master log directory and add a
    # download link. Written directly to the master log (not via append_log_section), so the
    # './...' link resolves relative to the consolidated report's own location.
    create_download_link(
      data = integrated_wide %>% select(StationID, Benthic, Chemistry, Toxicity, `Site Assessment`),
      logfile = logfile,
      filename = 'integrated-site-assessment-scores.csv',
      linktext = 'Download integrated site assessment scores',
      verbose = verbose
    )

    knit.log(logfile, verbose = verbose, knitlog = knitlog)
  }

  return(out)


}
