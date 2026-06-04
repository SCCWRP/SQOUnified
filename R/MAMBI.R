# GENERIC M-AMBI (Alt) ------------------------------------------------------------------------
#' Compute the multivariate AMBI (M-AMBI) index score (generic version).
#'
#' @description
#'   This function calculates the United States version of the multivariate AMBI (M-AMBI) index scores
#'   following Pelletier et al. 2018, which is in turn built upon the work of Sigovini et al. 2013
#'   and Muxica et al. 2007.
#'
#' @details
#'   The M-AMBI integrates three metrics via factor analysis:
#'   \itemize{
#'     \item AMBI score - abundance-weighted ecological group tolerance score
#'     \item Species Richness (S) - for saline samples; or percent Oligochaetes for tidal freshwater
#'     \item Shannon-Wiener Diversity (H') - base-2 Shannon diversity
#'   }
#'
#'   Samples are classified by salinity zone and compared against good/bad standards from
#'   Pelletier et al. 2018 for each zone. Factor analysis (PCA with varimax rotation) is used
#'   to produce M-AMBI scores via the Ecological Quality Ratio (EQR) transformation.
#'
#'   Salinity zones:
#'   \itemize{
#'     \item TF: Tidal Freshwater (0 - 0.2 PSU) - uses AMBI, H', and percent Oligochaetes
#'     \item OH: Oligohaline (0.2 - 5 PSU)
#'     \item MH: Mesohaline (5 - 18 PSU)
#'     \item PH/WPH: Polyhaline (18 - 30 PSU)
#'     \item EH/WEH: Euhaline (30 - 40 PSU)
#'     \item HH: Hyperhaline (> 40 PSU)
#'   }
#'
#'   M-AMBI condition categories (SQO):
#'   \itemize{
#'     \item Reference: score >= 0.578
#'     \item Low Disturbance: 0.483 <= score < 0.578
#'     \item Moderate Disturbance: 0.387 < score < 0.483
#'     \item High Disturbance: score <= 0.387
#'   }
#'
#'   AMBI applicability guidelines follow Borja and Muxica 2005:
#'   \itemize{
#'     \item <= 20 percent unassigned abundance: AMBI applicable
#'     \item 20-50 percent: Use with care
#'     \item > 50 percent: Not recommended
#'   }
#'
#' @param benthic_data a data frame containing benthic data and station information with at minimum:
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the sampling location;
#'
#'    \strong{\code{replicate}} - a numeric identifying the replicate number;
#'
#'    \strong{\code{sampledate}} - the date of sample collection;
#'
#'    \strong{\code{latitude}} - latitude in decimal degrees;
#'
#'    \strong{\code{longitude}} - longitude in decimal degrees (negative for west);
#'
#'    \strong{\code{salinity}} - bottom water salinity in PSU;
#'
#'    \strong{\code{taxon}} - name of the organism (SCAMIT ed14 format preferred).
#'        Use \code{NoOrganismsPresent} with 0 abundance for empty samples;
#'
#'    \strong{\code{abundance}} - number of individuals counted;
#'
#'    \strong{\code{exclude}} - "Yes" or "No" indicating if the taxon name is ambiguous.
#'
#' @param EG_Ref_values a data frame with Ecological Group assignments. If NULL (default),
#'    uses the standard EG list based on Gillett et al. 2015, updated with each Bight cycle.
#' @param EG_Scheme a quoted string specifying which EG value column to use. Default is "Hybrid".
#'    Other options: "US_East", "US_West", "US_Gulf", "US", "Standard".
#' @param logfile Path to a logfile. Default is an RMarkdown file in a timestamped logs directory.
#' @param verbose Logical. If TRUE, detailed logging output is produced. Default FALSE.
#' @param knitlog Logical. If TRUE, the log file is knitted to HTML upon completion. Default FALSE.
#'
#' @usage
#' MAMBI(benthic_data)
#'
#' @examples
#' \dontrun{
#'   MAMBI(my_benthic_data)
#' }
#'
#' @import dplyr
#' @import vegan
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stringr str_detect str_replace
#' @importFrom purrr map list_rbind
#'
#' @export
MAMBI <- function(benthic_data,
                              EG_Ref_values = NULL,
                              EG_Scheme = "Hybrid",
                              logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'MAMBI_generic_log.Rmd'),
                              verbose = FALSE,
                              knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n## BEGIN: Generic M-AMBI function.\n', logfile = logfile, verbose = verbose)

  # Reference data (Saline_Standards, TidalFresh_Standards, us.mambi.eg.values.04_23_24) are available as package datasets
  # EQR() is available in the package namespace from R/EQR.R

  # If no custom EG reference values provided, use the standard list
  if (is.null(EG_Ref_values)) {
    EG_Ref_values <- us.mambi.eg.values.04_23_24
  }

  Input_File.0 <- benthic_data %>%
    mutate(
      Species_ended_in_sp = (str_detect(taxon, " sp$")),
      taxon = (str_replace(taxon, " sp$", ""))) %>%
    #splitting out west coast sites to facilitate use of different criteria
    mutate(
      Coast = (ifelse(longitude <= -115, "West", "Gulf-East"))) %>%
    # associating a sample with one of the different salinity zones
    mutate(
      SalZone = case_when(
        salinity > 30 & salinity <= 40 & Coast == "Gulf-East" ~ "EH",
        salinity > 18 & salinity <= 30 & Coast == "Gulf-East" ~ "PH",
        salinity > 5 & salinity <= 18 ~ "MH",
        salinity > 0.2 & salinity <= 5 ~ "OH",
        salinity >= 0 & salinity <= 0.2 ~ "TF",
        salinity > 40 ~ "HH",
        salinity > 30 & salinity <= 40 & Coast == "West" ~ "WEH",
        salinity > 18 & salinity <= 30 & Coast == "West" ~ "WPH"
      )
    )

  EG_to_use <- EG_Ref_values %>%
    select(all_of(EG_Scheme), Taxon, Exclude) %>%
    rename(EG = all_of(EG_Scheme)) %>%
    relocate(Taxon, Exclude, EG) %>%
    mutate(EG = ifelse(Taxon == "Oligochaeta", "V", EG))

  #in case a sample had no animals, force it into the High Disturbance category
  defaunated <- Input_File.0 %>%
    filter(taxon == "NoOrganismsPresent") %>%
    mutate(mambi_score = 0, Orig_mambi_condition = "Bad", SQO_mambi_condition = "High Disturbance", ambi_score = 7, H = NaN, S = NaN, oligo_pct = NaN,
           use_ambi = "Yes", use_mambi = "Yes", note = "Defaunated Sample") %>%
    select(stationid, sampledate, replicate, latitude, longitude, mambi_score, Orig_mambi_condition, SQO_mambi_condition,
           ambi_score, H, S, oligo_pct, use_mambi, use_ambi, note)

  # dropping azoic samples from analysis
  Input_File <- Input_File.0 %>% filter(taxon != "NoOrganismsPresent")

  Sample.info <- Input_File %>%
    select(stationid, replicate, sampledate, latitude, longitude, salinity, Coast, SalZone) %>%
    distinct()

  Input_File2 <- Input_File %>%
    filter(!is.na(SalZone))

  EG.Assignment <- Input_File %>%
    left_join(., EG_to_use, by = c("taxon" = "Taxon")) %>%
    group_by(stationid, replicate, sampledate) %>%
    mutate(tot_abun = sum(abundance)) %>%
    ungroup() %>%
    mutate(rel_abun = ((abundance / tot_abun) * 100))

  # Export interim file with all taxa with an assigned EG value
  taxa_w_EG <- EG.Assignment %>%
    filter(EG %in% c("I", "II", "III", "IV", "V")) %>%
    group_by(taxon, EG) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop_last")

  writelog(
    '### M-AMBI Step 1 - Taxa with EG values\n',
    logfile = logfile,
    data = taxa_w_EG %>% head(25),
    verbose = verbose
  )
  create_download_link(data = taxa_w_EG, logfile = logfile, filename = 'MAMBI_generic-step1-taxa_with_EG.csv', linktext = 'Download taxa with EG values', verbose = verbose)

  # Export interim file with all taxa without an assigned EG value
  taxa_wo_EG <- EG.Assignment %>%
    filter(is.na(EG) | EG == "") %>%
    group_by(taxon) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop_last")

  writelog(
    '### M-AMBI Step 2 - Taxa without EG values\n',
    logfile = logfile,
    data = taxa_wo_EG %>% head(25),
    verbose = verbose
  )
  create_download_link(data = taxa_wo_EG, logfile = logfile, filename = 'MAMBI_generic-step2-taxa_without_EG.csv', linktext = 'Download taxa without EG values', verbose = verbose)


  # Calculating applicability of the index to each sample
  AMBI.applicability <- EG.Assignment %>%
    mutate(EG_Test = ifelse(is.na(EG), "NoEG", "YesEG")) %>%
    pivot_wider(id_cols = c(stationid, replicate, sampledate), names_from = EG_Test, values_from = rel_abun, values_fill = 0, values_fn = sum) %>%
    mutate(use_ambi = case_when(
      NoEG <= 20 ~ "Yes",
      NoEG > 20 & NoEG <= 50 ~ "With Care",
      NoEG > 50 ~ "Not Recommended",
      is.na(NoEG) ~ "Yes"))

  MAMBI.applicability <- AMBI.applicability %>%
    left_join(., Sample.info, by = c("stationid", "replicate", "sampledate")) %>%
    mutate(use_mambi = case_when(is.na(SalZone) ~ "No - No Salinity Value",
                                 use_ambi == "With Care" ~ "Caution - Sparse AMBI Coverage",
                                 use_ambi == "Not Reccommended" ~ "Not Recommended - Poor AMBI Coverage",
                                 TRUE ~ "Yes")) %>%
    select(stationid, replicate, sampledate, use_ambi, use_mambi)


  # establish the salinity zones in the data set
  Sal_range.dataset <- unique(Input_File2$SalZone)


  ######Saline calcs ################
  # Calculating AMBI for each sample
  AMBI.Scores <- EG.Assignment %>%
    group_by(stationid, replicate, sampledate, tot_abun, EG) %>%
    summarise(sum_rel = sum(rel_abun), .groups = "drop_last") %>%
    ungroup() %>%
    replace_na(list(EG = "NoEG")) %>%
    mutate(
      EG_Score = case_when(
        EG == "I" ~ sum_rel * 0,
        EG == "II" ~ sum_rel * 1.5,
        EG == "III" ~ sum_rel * 3,
        EG == "IV" ~ sum_rel * 4.5,
        EG == "V" ~ sum_rel * 6,
        EG == "NoEG" ~ 0)) %>%
    mutate(EG_Score = ifelse(tot_abun == 0, 7, EG_Score)) %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(ambi_score = (sum(EG_Score, na.rm = TRUE) / 100), .groups = "drop_last") %>%
    ungroup()

  #Calculating taxa richness for each sample
  Rich <- Input_File %>%
    filter(exclude != "Yes") %>%
    group_by(stationid, replicate, sampledate) %>%
    summarise(S = length(taxon), .groups = "drop_last") %>%
    ungroup()

  # Calculating base 2 shannon wiener diversity for each sample
  Divy <- Input_File %>%
    filter(exclude != "Yes") %>%
    group_by(stationid, replicate, sampledate, taxon) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    pivot_wider(id_cols = c(stationid, replicate, sampledate), names_from = taxon, values_from = abundance, values_fill = 0) %>%
    mutate(H = vegan::diversity(select(., 4:ncol(.)), index = "shannon", base = 2)) %>%
    select(stationid, replicate, sampledate, H)

  # Combining the three metrics
  metrics <- AMBI.Scores %>%
    left_join(., Rich, by = c("stationid", "replicate", "sampledate")) %>%
    left_join(., Divy, by = c("stationid", "replicate", "sampledate"))

  # Combining the metrics with the station information
  metrics.1 <- Sample.info %>%
    left_join(., metrics, by = c("stationid", "replicate", "sampledate")) %>%
    select(stationid, replicate, sampledate, ambi_score, S, H, SalZone)

  # Isolating samples for which salinity has not been submitted
  no.SalZone.data <- metrics.1 %>%
    filter(is.na(SalZone)) %>%
    left_join(., Sample.info, by = c("stationid", "replicate", "sampledate", "SalZone")) %>%
    select(stationid, replicate, sampledate, ambi_score, S, H, latitude, longitude, SalZone) %>%
    mutate(mambi_score = NA_real_, Orig_mambi_condition = NA_character_, SQO_mambi_condition = NA_character_)

  # Adding the standard good/bad endpoints for each salinity zone
  metrics.2 <- metrics.1 %>%
    filter(!is.na(SalZone)) %>%
    bind_rows(., Saline_Standards)

  # Run factor analysis across the different salinity zones
  saline.mambi <- purrr::map(Sal_range.dataset, function(sal)
  {
    sal.df <- filter(metrics.2, SalZone == sal)
    METRICS.tot <- select(sal.df, ambi_score, S, H)

    options(warn = -1)
    METRICS.fa2 <- princomp(METRICS.tot, cor = T, covmat = cov(METRICS.tot))
    options(warn = 0)
    METRICS.fa2.load <- loadings(METRICS.fa2) %*% diag(METRICS.fa2$sdev)
    METRICS.fa2.load.varimax <- loadings(varimax(METRICS.fa2.load))
    METRICS.scores2 <- scale(METRICS.tot) %*% METRICS.fa2.load.varimax
    colnames(METRICS.scores2) <- c("x", "y", "z")
    METRICS.tr <- METRICS.scores2

    # Transforming factor analysis loadings into M-AMBI scores using the EQR function
    eqr <- EQR(METRICS.tr)
    colnames(eqr) <- c("mambi_score")
    eqr <- data.frame(eqr)

    # Adding back in station information and classifying scores
    results <- sal.df %>%
      bind_cols(., eqr) %>%
      left_join(., Sample.info, by = c("stationid", "replicate", "sampledate", "SalZone")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, SalZone, ambi_score, S, H, mambi_score) %>%
      filter(!stationid %in% Saline_Standards$stationid, SalZone != "TF") %>%
      mutate(
        Orig_mambi_condition = case_when(
          mambi_score < 0.2 ~ "Bad",
          mambi_score >= 0.2 & mambi_score < 0.39 ~ "Poor",
          mambi_score >= 0.39 & mambi_score < 0.53 ~ "Moderate",
          mambi_score >= 0.53 & mambi_score < 0.77 ~ "Good",
          mambi_score >= 0.77 ~ "High"
        ),
        SQO_mambi_condition = case_when(
          mambi_score <= 0.387 ~ "High Disturbance",
          mambi_score > 0.387 & mambi_score < 0.483 ~ "Moderate Disturbance",
          mambi_score >= 0.483 & mambi_score < 0.578 ~ "Low Disturbance",
          mambi_score >= 0.578 ~ "Reference")
      )
  }) %>% list_rbind()


  ###################
  # if there are tidal freshwater samples, run TF subroutine
  if (any(Sal_range.dataset == "TF"))
  {
    TF.EG.Assignment <- EG.Assignment %>% filter(SalZone == "TF")
    TF.EG_Ref_values <- us.mambi.eg.values.04_23_24 %>% select(., Taxon, Exclude, EG = all_of(EG_Scheme), Oligochaeta)

    # calculate AMBI scores for each sample
    TF.AMBI.Scores <- TF.EG.Assignment %>%
      group_by(stationid, replicate, sampledate, tot_abun, EG, rel_abun) %>%
      summarise(sum_rel = sum(rel_abun), .groups = "drop_last") %>%
      ungroup() %>%
      replace_na(list(EG = "NoEG")) %>%
      mutate(
        EG_Score = case_when(
          EG == "I" ~ sum_rel * 0,
          EG == "II" ~ sum_rel * 1.5,
          EG == "III" ~ sum_rel * 3,
          EG == "IV" ~ sum_rel * 4.5,
          EG == "V" ~ sum_rel * 6,
          EG == "NoEG" ~ 0)) %>%
      group_by(stationid, replicate, sampledate) %>%
      summarise(ambi_score = (sum(EG_Score) / 100), .groups = "drop_last") %>%
      ungroup()

    # calculate % oligochaetes in each sample
    TF.Oligos <- Input_File %>%
      group_by(stationid, replicate, sampledate) %>%
      mutate(tot_abun = sum(abundance)) %>%
      ungroup() %>%
      left_join(., TF.EG_Ref_values, by = c("taxon" = "Taxon")) %>%
      filter(Oligochaeta == "Yes", SalZone == "TF") %>%
      group_by(stationid, replicate, sampledate) %>%
      summarise(oligo_pct = sum(abundance / tot_abun) * 100, .groups = "drop_last") %>%
      ungroup()

    # calculate shannon wiener diversity (base 2) for each sample
    TF.Divy <- Input_File %>%
      filter(SalZone == "TF") %>%
      group_by(stationid, replicate, sampledate, taxon) %>%
      summarise(abundance = sum(abundance), .groups = "drop") %>%
      pivot_wider(id_cols = c(stationid, replicate, sampledate), names_from = taxon, values_from = abundance, values_fill = 0) %>%
      mutate(H = diversity((select(., 4:(ncol(.)))), index = "shannon", base = 2)) %>%
      select(., stationid, replicate, sampledate, H)

    # join the three metrics together
    TF.metrics <- TF.AMBI.Scores %>%
      left_join(., TF.Divy, by = c("stationid", "replicate", "sampledate")) %>%
      left_join(., TF.Oligos, by = c("stationid", "replicate", "sampledate"))

    # add on sample information
    TF.metrics.1 <- Sample.info %>%
      filter(SalZone == "TF") %>%
      left_join(., TF.metrics, by = c("stationid", "replicate", "sampledate")) %>%
      select(stationid, replicate, sampledate, ambi_score, H, oligo_pct, SalZone)

    # add the best and worst standard values for tidal freshwater
    TF.metrics.2 <- bind_rows(TF.metrics.1, TidalFresh_Standards)

    TF.METRICS.tot <- TF.metrics.2 %>%
      select(ambi_score, H, oligo_pct)

    # factor analysis for tidal freshwater
    options(warn = -1)
    TF.METRICS.fa2 <- princomp(TF.METRICS.tot, cor = T, covmat = cov(TF.METRICS.tot))
    options(warn = 0)
    TF.METRICS.fa2.load <- loadings(TF.METRICS.fa2) %*% diag(TF.METRICS.fa2$sdev)
    TF.METRICS.fa2.load.varimax <- loadings(varimax(TF.METRICS.fa2.load))
    TF.METRICS.scores2 <- scale(TF.METRICS.tot) %*% TF.METRICS.fa2.load.varimax
    colnames(TF.METRICS.scores2) <- c("x", "y", "z")
    TF.METRICS.tr <- TF.METRICS.scores2

    # EQR transformation
    TF.eqr <- EQR(TF.METRICS.tr)
    colnames(TF.eqr) <- c("mambi_score")
    TF.eqr <- data.frame(TF.eqr)

    # classify TF scores
    TF.mambi <- TF.metrics.2 %>%
      bind_cols(., TF.eqr) %>%
      left_join(., Sample.info, by = c("stationid", "replicate", "SalZone", "sampledate")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, ambi_score, H, oligo_pct, mambi_score) %>%
      filter(!stationid %in% TidalFresh_Standards$stationid) %>%
      mutate(
        Orig_mambi_condition = case_when(
          mambi_score < 0.2 ~ "Bad",
          mambi_score >= 0.2 & mambi_score < 0.39 ~ "Poor",
          mambi_score >= 0.39 & mambi_score < 0.53 ~ "Moderate",
          mambi_score >= 0.53 & mambi_score < 0.77 ~ "Good",
          mambi_score >= 0.77 ~ "High"),
        SQO_mambi_condition = case_when(
          mambi_score <= 0.387 ~ "High Disturbance",
          mambi_score > 0.387 & mambi_score < 0.483 ~ "Moderate Disturbance",
          mambi_score >= 0.483 & mambi_score < 0.578 ~ "Low Disturbance",
          mambi_score >= 0.578 ~ "Reference"))

    # combine saline and TF results
    saline.mambi.2 <- saline.mambi %>%
      mutate(oligo_pct = NA) %>%
      select(stationid, replicate, sampledate, latitude, longitude, ambi_score, S, H, oligo_pct, mambi_score,
             Orig_mambi_condition, SQO_mambi_condition, SalZone)

    TF.mambi.2 <- TF.mambi %>%
      mutate(S = NA, SalZone = "TF") %>%
      select(stationid, replicate, sampledate, latitude, longitude, ambi_score, S, H, oligo_pct, mambi_score,
             Orig_mambi_condition, SQO_mambi_condition, SalZone)

    Overall.Results <- bind_rows(saline.mambi.2, TF.mambi.2, no.SalZone.data) %>%
      left_join(., MAMBI.applicability, by = c("stationid", "replicate", "sampledate")) %>%
      mutate(note = case_when(is.na(SalZone) ~ "No salinity data - cannot calculate M-AMBI",
                              use_mambi == "Yes" ~ "None",
                              use_mambi == "Caution - Sparse AMBI Coverage" ~ "Check sample - M-AMBI ok, but limited taxa coverage",
                              use_mambi == "Not Recommended - Poor AMBI Coverage" ~ "M-AMBI not recommended - poor taxa coverage")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, mambi_score, Orig_mambi_condition, SQO_mambi_condition,
             ambi_score, H, S, oligo_pct, use_mambi, use_ambi, note) %>%
      bind_rows(., defaunated) %>%
      mutate(index = "M-AMBI", .before = latitude)
  }
  else
  {
    #if there are no tidal freshwater samples, only the saline calculations are needed
    Overall.Results <- saline.mambi %>%
      bind_rows(., no.SalZone.data) %>%
      left_join(., MAMBI.applicability, by = c("stationid", "replicate", "sampledate")) %>%
      mutate(oligo_pct = NaN,
             note = case_when(is.na(SalZone) ~ "No salinity data - cannot calculate M-AMBI",
                              use_mambi == "Yes" ~ "None",
                              use_mambi == "Caution - Sparse AMBI Coverage" ~ "Check sample - M-AMBI ok, but limited taxa coverage",
                              use_mambi == "Not Recommended - Poor AMBI Coverage" ~ "M-AMBI not recommended - poor taxa coverage")) %>%
      select(stationid, replicate, sampledate, latitude, longitude, mambi_score, Orig_mambi_condition, SQO_mambi_condition,
             ambi_score, H, S, oligo_pct, use_mambi, use_ambi, note) %>%
      bind_rows(., defaunated) %>%
      mutate(index = "M-AMBI", .before = latitude)
  }

  #gathering all of the site/sample information
  station.info <- benthic_data %>%
    select(-taxon, -abundance, -exclude) %>%
    distinct()

  Overall.Results.2 <- Overall.Results %>%
    left_join(., station.info, by = c("stationid", "sampledate", "replicate", "latitude", "longitude"))

  writelog(
    '### M-AMBI Final - M-AMBI Scores\n',
    logfile = logfile,
    data = Overall.Results.2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = Overall.Results.2, logfile = logfile, filename = 'MAMBI_generic-final_scores.csv', linktext = 'Download M-AMBI scores', verbose = verbose)

  writelog('\n## END: Generic M-AMBI function.\n', logfile = logfile, verbose = verbose)

  return(Overall.Results.2)
}


# This function is private, only used by MAMBI
EQR <- function(data, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), '-log.txt') ), verbose = F) {
  segm <- data[nrow(data),] - data[(nrow(data)-1),]
  vett <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
  for (k in 1: ncol(data)) {vett[, k] <- data[(nrow(data)-1), k]}
  vett <- data - vett
  ris <- round((vett %*% segm / sqrt(sum(segm*segm))) / sqrt(sum(segm*segm)), 3)
  return(ris)
}
