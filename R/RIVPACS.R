#' CRiver Invertebrate Prediction and Classification System (RIVPACS) Index and RIVPACS Condition Category
#' (SoCal only)
#'
#' @description
#'   For more information concerning RIVPACS, consult the CASQO Technical Manual page 80
#'
#' @param benthic_data data frame stored in the R environment. Note that this data frame MUST contain the following
#'                    information with these headings:
#'
#'                         \code{StationID} - an alpha-numeric identifier of the location;
#'
#'                         \code{Replicate} - a numeric identifying the replicate number of samples taken at the location;
#'
#'                         \code{SampleDate} - the date of sample collection;
#'
#'                         \code{Latitude} - latitude in decimal degrees;
#'
#'                         \code{Longitude} - longitude in decimal degrees. Make sure there is a negative sign for the Western coordinates;
#'
#'                         \code{Taxon} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#'                         \code{Abundance} - the number of each Species observed in a sample;
#'
#'                         \code{Salinity} - the salinity observed at the location in PSU, ideally at time of sampling;
#'
#'                         \code{Stratum} - ;
#'
#'                         \code{Exclude} - ;
#'
#' @usage
#' RIVPACS(benthic_data)
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' data(benthic_sampledata)
#' RIVPACS(benthic_sampledata)
#'
#' @export
#'

#---- RIVPACS WRAPPER FUNCTION ----
# This is what we will use for RIVPACS
RIVPACS <- function(benthic_data, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = F){

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\n## BEGIN: RIVPACS function.\n', logfile = logfile, verbose = verbose)


  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'benthic.rivpacs.input.RData'
  if (verbose) {
    save(benthic_data, file = file.path( dirname(logfile), rawinput.filename ))
  }

  # Create code block and download link to RIVPACS input
  writelog(
    'Input to RIVPACS - RIVPACS-step0.csv',
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') ### This will load a dataframe called 'benthic_data' into your environment"),
    data = benthic_data,
    verbose = verbose
  )
  create_download_link(data = benthic_data, logfile = logfile, filename = 'RIVPACS-step0.csv', linktext = 'Download RIVPACS initial input', verbose = verbose)





  # Split to SoCal and SFBay.
  ## We are only working with SoCal data so we don't need to do this!


  # SCB Predictors - needs to be logged
  scb.predictors <- data.frame(Latitude = benthic_data$Latitude,
                               Longitude = benthic_data$Longitude,
                               SampleDepth = benthic_data$SampleDepth) %>%
    dplyr::distinct()
  # Create code block and download link to BRI input
  writelog(
    'SCB Predictors',
    logfile = logfile,
    code = '
      scb.predictors <- data.frame(Latitude = benthic_data$Latitude,
                               Longitude = benthic_data$Longitude,
                               SampleDepth = benthic_data$SampleDepth) %>%
        dplyr::distinct()
    ',
    data = scb.predictors,
    verbose = verbose
  )
  create_download_link(data = scb.predictors, logfile = logfile, filename = 'RIVPACS-SCB.Predictors.csv', linktext = 'Download RIVPACS initial input', verbose = verbose)




  # RIVPACS Data Prep step 1 - rename taxa to Taxon
  benthic_data <- benthic_data %>% dplyr::rename(Taxa = Taxon)
  # Write to the logs for RIVPACS Data Prep step 1
  writelog(
    '\nRIVPACS Data Prep step 1 - rename taxa to Taxon',
    logfile = logfile,
    code = "
      benthic_data <- benthic_data %>% dplyr::rename(Taxa = Taxon)",
    data = benthic_data,
    verbose = verbose
  )
  create_download_link(data = benthic_data, logfile = logfile, filename = 'RIVPACS-DataPrep-step1.csv', linktext = 'Download RIVPACS Data Prep step 1', verbose = verbose)



  # RIVPACS data prep step 2 - get distinct records on StationID, Latitude, Longitude, SampleDepth - also set row names to StationID
  scb.taxa <- benthic_data %>% dplyr::select(StationID, Latitude, Longitude, SampleDepth) %>%
    dplyr::distinct()

  # Write to the logs for RIVPACS data prep step 2
  writelog(
    '\nRIVPACS data prep step 2 - get distinct records on StationID, Latitude, Longitude, SampleDepth',
    logfile = logfile,
    code = "
      scb.taxa <- benthic_data %>%
        dplyr::select(StationID, Latitude, Longitude, SampleDepth) %>%
        dplyr::distinct()
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step2.csv', linktext = 'Download RIVPACS Data Prep step 2', verbose = verbose)



  # RIVPACS data prep step 3 - set up the scb predictors
  row.names(scb.predictors) <- scb.taxa$StationID
  scb.predictors <- as.matrix(scb.predictors)

  # Write to the logs for RIVPACS data prep step 3
  writelog(
    '\nRIVPACS data prep step 3 - set up the scb predictors',
    logfile = logfile,
    code = "
      row.names(scb.predictors) <- scb.taxa$StationID
      scb.predictors <- as.matrix(scb.predictors)
    ",
    data = scb.predictors,
    verbose = verbose
  )
  create_download_link(data = scb.predictors, logfile = logfile, filename = 'RIVPACS-DataPrep-step3.csv', linktext = 'Download RIVPACS Data Prep step 3', verbose = verbose)




  # RIVPACS data prep step 4 - Filter to replicate one and get distinct values on StationID Taxa and Abundance
  scb.taxa <- benthic_data %>%
    dplyr::filter(Replicate == 1) %>%
    dplyr::select(StationID, Taxa, Abundance) %>%
    dplyr::distinct()
  # Write to the logs for RIVPACS data prep step 4
  writelog(
    '\nRIVPACS data prep step 4 - Filter to replicate one and get distinct values on StationID, Taxa, and Abundance',
    logfile = logfile,
    code = "
      scb.taxa <- benthic_data %>%
        dplyr::filter(Replicate == 1) %>%
        dplyr::select(StationID, Taxa, Abundance) %>%
        dplyr::distinct()
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step4.csv', linktext = 'Download RIVPACS Data Prep step 4', verbose = verbose)




  # RIVPACS Data prep step 5 - remove certain special characters from taxa name
  scb.taxa$Taxa <- gsub(" ", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub("(", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub(")", "_", scb.taxa$Taxa, fixed = TRUE)

  # Write to the logs for RIVPACS data prep step 5
  writelog(
    '\nRIVPACS Data prep step 5 - remove certain special characters from taxa name',
    logfile = logfile,
    code = "
      scb.taxa$Taxa <- gsub(' ', '_', scb.taxa$Taxa, fixed = TRUE)
      scb.taxa$Taxa <- gsub('(', '_', scb.taxa$Taxa, fixed = TRUE)
      scb.taxa$Taxa <- gsub(')', '_', scb.taxa$Taxa, fixed = TRUE)
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step5.csv', linktext = 'Download RIVPACS Data Prep step 5', verbose = verbose)




  # RIVPACS Data prep step 6 - pivot the data out wide and make it a data.frame
  scb.taxa <- scb.taxa %>%
    tidyr::pivot_wider(id_cols = "StationID", names_from = "Taxa",
                       values_from = "Abundance", values_fn = list(Abundance = list))
  scb.taxa <- as.data.frame(scb.taxa)

  # Write to the logs for RIVPACS data prep step 6
  writelog(
    '\nRIVPACS Data prep step 6 - pivot the data out wide and make it a data.frame',
    logfile = logfile,
    code = "
      scb.taxa <- scb.taxa %>%
        tidyr::pivot_wider(id_cols = 'StationID', names_from = 'Taxa',
                           values_from = 'Abundance', values_fn = list(Abundance = list))
      scb.taxa <- as.data.frame(scb.taxa)
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step6.csv', linktext = 'Download RIVPACS Data Prep step 6', verbose = verbose)



  # RIVPACS Data prep step 7
  scb.taxa <- scb.taxa[, -1]
  # Write to the logs for RIVPACS data prep step 7
  writelog(
    '\nRIVPACS Data prep step 7',
    logfile = logfile,
    code = "
      scb.taxa <- scb.taxa[, -1]
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step7.csv', linktext = 'Download RIVPACS Data Prep step 7', verbose = verbose)


  # RIVPACS data prep step 8 - remove Abundance. from column names
  colnames(scb.taxa) <- gsub("Abundance.", "", colnames(scb.taxa))
  # Write to the logs for RIVPACS data prep step 8
  writelog(
    '\nRIVPACS data prep step 8 - remove Abundance. from column names',
    logfile = logfile,
    code = "
      colnames(scb.taxa) <- gsub('Abundance.', '', colnames(scb.taxa))
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step8.csv', linktext = 'Download RIVPACS Data Prep step 8', verbose = verbose)




  # RIVPACS data prep step 9 - Replace NAs with zero.
  scb.taxa[scb.taxa == "NULL"] <- 0
  scb.taxa = as.data.frame(lapply(scb.taxa, as.numeric))
  row.names(scb.taxa) <- row.names(scb.predictors)
  # Write to the logs for RIVPACS data prep step 9
  writelog(
    '\nRIVPACS data prep step 9 - Replace NAs with zero.',
    logfile = logfile,
    code = "
      scb.taxa[scb.taxa == 'NULL'] <- 0
      scb.taxa <- as.data.frame(lapply(scb.taxa, as.numeric))
      row.names(scb.taxa) <- row.names(scb.predictors)
    ",
    data = scb.taxa,
    verbose = verbose
  )
  create_download_link(data = scb.taxa, logfile = logfile, filename = 'RIVPACS-DataPrep-step9.csv', linktext = 'Download RIVPACS Data Prep step 9', verbose = verbose)



  # RIVPACS calculations. By default the functions use the example user data.
  socal <- SoCalRivpacs(observed.predictors = scb.predictors, observed.taxa = scb.taxa, logfile = logfile, verbose = verbose)
  # Write to the logs for RIVPACS calculations
  writelog(
    '\nRIVPACS calculations. By default the functions use the example user data.',
    logfile = logfile,
    code = "
      socal <- SoCalRivpacs(observed.predictors = scb.predictors, observed.taxa = scb.taxa, logfile = logfile, verbose = verbose)
    ",
    data = socal,
    verbose = verbose
  )
  create_download_link(data = socal, logfile = logfile, filename = 'RIVPACS-calculations.csv', linktext = 'Download RIVPACS calculations', verbose = verbose)




  # the stations column of the oe table dataframe was being returned as a factor. Need to make that a character
  socal$oe.table <- socal$oe.table %>%
    mutate_if(is.factor,as.character)

  # Write to the logs for converting stations column to character in oe table dataframe
  writelog(
    '\nThe stations column of the oe table dataframe was being returned as a factor. Need to make that a character',
    logfile = logfile,
    code = "
      socal$oe.table <- socal$oe.table %>%
        mutate_if(is.factor, as.character)
    ",
    data = socal$oe.table,
    verbose = verbose
  )
  create_download_link(data = socal$oe.table, logfile = logfile, filename = 'RIVPACS-socal.oe-table.csv', linktext = 'Download RIVPACS socal Observed/Expected table', verbose = verbose)



  # Get the distinct values in benthic data based on StationID Replicate SampleDate Stratum
  benthic_data <- benthic_data %>%
    dplyr::select(StationID, Replicate, SampleDate, Stratum) %>%
    dplyr::distinct()

  # Write to the logs for getting distinct values in benthic data based on StationID, Replicate, SampleDate, Stratum
  writelog(
    '\nGet the distinct values in benthic data based on StationID, Replicate, SampleDate, Stratum',
    logfile = logfile,
    code = "
      benthic_data <- benthic_data %>%
        dplyr::select(StationID, Replicate, SampleDate, Stratum) %>%
        dplyr::distinct()
    ",
    data = benthic_data,
    verbose = verbose
  )
  create_download_link(data = benthic_data, logfile = logfile, filename = 'benthic-data-distinct.csv', linktext = 'Download distinct benthic data', verbose = verbose)






  # Calculate RIVPACS Scores
  # Riv step 0 - select stations and Observed/Expected
  riv0 <- socal$oe.table %>%
    dplyr::select(stations, O.over.E)
  # Write to the logs for RIVPACS Scores calculation step 0
  writelog(
    '\nCalculate RIVPACS Scores\nRiv step 0 - select stations and Observed/Expected',
    logfile = logfile,
    code = "
    riv0 <- socal$oe.table %>%
      dplyr::select(stations, O.over.E)
  ",
    data = riv0,
    verbose = verbose
  )
  create_download_link(data = riv0, logfile = logfile, filename = 'RIVPACS-Scores-step0.csv', linktext = 'Download RIVPACS Scores step 0', verbose = verbose)




  # Riv step 1 - join with benthic data
  riv1 <- riv0 %>%
    dplyr::rename(StationID = stations, Score = O.over.E) %>%
    dplyr::full_join(benthic_data) %>%
    dplyr::mutate(Index = "RIVPACS")

  # Write to the logs for RIVPACS Scores calculation step 1
  writelog(
    '\nRiv step 1 - join with benthic data',
    logfile = logfile,
    code = "
      riv1 <- riv0 %>%
        dplyr::rename(StationID = stations, Score = O.over.E) %>%
        dplyr::full_join(benthic_data) %>%
        dplyr::mutate(Index = 'RIVPACS')
    ",
    data = riv1,
    verbose = verbose
  )
  create_download_link(data = riv1, logfile = logfile, filename = 'RIVPACS-Scores-step1.csv', linktext = 'Download RIVPACS Scores step 1', verbose = verbose)






  # Get the scores based on the thresholds page 73 table 4.25 - https://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf
  rivpacs.score <- riv1 %>%
    dplyr::mutate(Category = case_when((Score <= 0.32) ~ "High Disturbance",
                                       ((Score > 0.32 & Score <= 0.74) | (Score >= 1.26)) ~ "Moderate Disturbance",
                                       ((Score > 0.74 & Score <= 0.90) | Score >= 1.10 & Score < 1.26) ~ "Low Disturbance",
                                       (Score > 0.90 | Score < 1.10) ~ "Reference")) %>%
    dplyr::mutate(`Category Score` = case_when(Category == "Reference" ~ 1,
                                               Category == "Low Disturbance" ~ 2,
                                               Category == "Moderate Disturbance" ~ 3,
                                               Category == "High Disturbance" ~ 4)) %>%
    dplyr::select(StationID, SampleDate, Replicate, Stratum, Index, Score, Category, `Category Score`)
  # Write to the logs for getting scores based on thresholds
  writelog(
    '\nGet the scores based on the thresholds page 73 table 4.25 - https://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf',
    logfile = logfile,
    code = "
      rivpacs.score <- riv1 %>%
        dplyr::mutate(Category = case_when(
          Score <= 0.32 ~ 'High Disturbance',
          (Score > 0.32 & Score <= 0.74) | (Score >= 1.26) ~ 'Moderate Disturbance',
          (Score > 0.74 & Score <= 0.90) | (Score >= 1.10 & Score < 1.26) ~ 'Low Disturbance',
          Score > 0.90 & Score < 1.10 ~ 'Reference'
        )) %>%
        dplyr::mutate(`Category Score` = case_when(
          Category == 'Reference' ~ 1,
          Category == 'Low Disturbance' ~ 2,
          Category == 'Moderate Disturbance' ~ 3,
          Category == 'High Disturbance' ~ 4
        )) %>%
        dplyr::select(StationID, SampleDate, Replicate, Stratum, Index, Score, Category, `Category Score`)
    ",
    data = rivpacs.score,
    verbose = verbose
  )
  create_download_link(data = rivpacs.score, logfile = logfile, filename = 'RIVPACS-Scores.csv', linktext = 'Download RIVPACS Scores', verbose = verbose)






  writelog('\n## END: RIVPACS function.\n', logfile = logfile, verbose = verbose)


  return(rivpacs.score)
}


# ---- ORIGINAL RIVPACS FUNCTION ----
# This is the original RIVPACS function. Called by the above function.
SoCalRivpacs <- function(Pcutoff = 0.5,
                         reference.groups = socal.reference.groups,
                         observed.predictors = socal.example.habitat,
                         reference.taxa = socal.reference.taxa,
                         group.means = socal.reference.group.means,
                         reference.cov = socal.reference.covariance,
                         observed.taxa = socal.example.taxa,
                         logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
                         verbose = F
                         ) {

  writelog('\n### BEGIN: So Cal RIVPACS function.\n', logfile = logfile, verbose = verbose)

  # set up the hyphen log prefix - which hasnt yet worked as i want it to
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)



  # Save the reference groups to an RData file
  tmp.filename <- 'benthic.rivpacs.socal.reference.groups.RData'
  if (verbose) {
    save(reference.groups, file = file.path( dirname(logfile), tmp.filename ))
  }

  # Create code block and download link to reference.groups
  writelog(
    'SoCal RIVPACS Reference Groups',
    logfile = logfile,
    code = paste0("load('", tmp.filename, "') ### This will load a dataframe called 'reference.groups' into your environment"),
    data = reference.groups,
    verbose = verbose
  )
  create_download_link(data = reference.groups, logfile = logfile, filename = 'SoCalRIVPACS-reference.groups.csv', linktext = 'Download RIVPACS reference.groups', verbose = verbose)



  # Save the observed.predictors to an RData file
  tmp.filename <- 'benthic.rivpacs.socal.observed.predictors.RData'
  if (verbose) {
    save(observed.predictors, file = file.path( dirname(logfile), tmp.filename ))
  }
  # Create code block and download link to observed.predictors
  writelog(
    'SoCal RIVPACS Observed Predictors',
    logfile = logfile,
    code = paste0("load('", tmp.filename, "') ### This will load a dataframe called 'observed.predictors' into your environment"),
    data = observed.predictors,
    verbose = verbose
  )
  create_download_link(data = observed.predictors, logfile = logfile, filename = 'SoCalRIVPACS-observed.predictors.csv', linktext = 'Download RIVPACS observed.predictors', verbose = verbose)



  # Save the reference taxa to an RData file
  tmp.filename <- 'benthic.rivpacs.socal.reference.taxa.RData'
  if (verbose) {
    save(reference.taxa, file = file.path(dirname(logfile), tmp.filename))
  }
  # Create code block and download link to reference.taxa
  writelog(
    'SoCal RIVPACS Reference Taxa',
    logfile = logfile,
    code = paste0("load('", tmp.filename, "') ### This will load a dataframe called 'reference.taxa' into your environment"),
    data = reference.taxa,
    verbose = verbose
  )
  create_download_link(data = reference.taxa, logfile = logfile, filename = 'SoCalRIVPACS-reference.taxa.csv', linktext = 'Download RIVPACS reference.taxa', verbose = verbose)


  # Save the group means to an RData file
  tmp.filename <- 'benthic.rivpacs.socal.group.means.RData'
  if (verbose) {
    save(group.means, file = file.path(dirname(logfile), tmp.filename))
  }
  # Create code block and download link to group.means
  writelog(
    'SoCal RIVPACS Group Means',
    logfile = logfile,
    code = paste0("load('", tmp.filename, "') ### This will load a dataframe called 'group.means' into your environment"),
    data = group.means,
    verbose = verbose
  )
  create_download_link(data = group.means, logfile = logfile, filename = 'SoCalRIVPACS-group.means.csv', linktext = 'Download RIVPACS group.means', verbose = verbose)


  # Save the reference covariance to an RData file
  tmp.filename <- 'benthic.rivpacs.socal.reference.cov.RData'
  if (verbose) {
    save(reference.cov, file = file.path(dirname(logfile), tmp.filename))
  }
  # Create code block and download link to reference.cov
  writelog(
    'SoCal RIVPACS Reference Covariance',
    logfile = logfile,
    code = paste0("load('", tmp.filename, "') ### This will load a dataframe called 'reference.cov' into your environment"),
    data = reference.cov,
    verbose = verbose
  )
  create_download_link(data = reference.cov, logfile = logfile, filename = 'SoCalRIVPACS-reference.cov.csv', linktext = 'Download RIVPACS reference.cov', verbose = verbose)



  # Save the observed taxa to an RData file
  tmp.filename <- 'benthic.rivpacs.socal.observed.taxa.RData'
  if (verbose) {
    save(observed.taxa, file = file.path(dirname(logfile), tmp.filename))
  }
  # Create code block and download link to observed.taxa
  writelog(
    'SoCal RIVPACS Observed Taxa',
    logfile = logfile,
    code = paste0("load('", tmp.filename, "') ### This will load a dataframe called 'observed.taxa' into your environment"),
    data = observed.taxa,
    verbose = verbose
  )
  create_download_link(data = observed.taxa, logfile = logfile, filename = 'SoCalRIVPACS-observed.taxa.csv', linktext = 'Download RIVPACS observed.taxa', verbose = verbose)


  # Pcutoff is the probability cutoff
  # Log that so the user can see which value is being used
  writelog( paste0("\n SoCal RIVPACS Pcutoff: " , Pcutoff), logfile = logfile, verbose = verbose )


  # Names of predictor variables.
  predictor.variables <- c("Latitude", "Longitude", "SampleDepth")
  # Write to the logs for names of predictor variables
  writelog(
    '\nNames of predictor variables.',
    logfile = logfile,
    code = "
      predictor.variables <- c('Latitude', 'Longitude', 'SampleDepth')
    ",
    data = predictor.variables,
    verbose = verbose
  )
  create_download_link(data = predictor.variables, logfile = logfile, filename = 'predictor-variables.csv', linktext = 'Download predictor variables', verbose = verbose)





  # ----- Define function - Format Observed Data -----
  FormatObservedData <- function(
    logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
    verbose = F
  ) {

    # set up the hyphen log prefix - which hasnt yet worked as i want it to
    hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

    # Align observed (user) data columns with reference data columns. Columns in same
    # order. Observed data may have a different number of taxa (columns) than
    # reference data.

    # Convert observed.taxa to presence/absence (0/1)
    tmp.pa <- observed.taxa
    tmp.pa[tmp.pa > 0] <- 1
    # Write to the logs for converting observed.taxa to presence/absence (0/1)
    writelog(
      '\nConvert observed.taxa to presence/absence (0/1)',
      logfile = logfile,
      code = "
        tmp.pa <- observed.taxa
        tmp.pa[tmp.pa > 0] <- 1
      ",
      data = tmp.pa,
      verbose = verbose
    )
    create_download_link(data = tmp.pa, logfile = logfile, filename = 'observed-taxa-presence-absence.csv', linktext = 'Download observed taxa presence/absence', verbose = verbose)


    # Align rows using predictor variables.
    tmp.pa <- tmp.pa[row.names(observed.predictors), ]    # !!! is this required???
    # Write to the logs for aligning rows using predictor variables
    writelog(
      '\nAlign rows using predictor variables.',
      logfile = logfile,
      code = "
        tmp.pa <- tmp.pa[row.names(observed.predictors), ]
      ",
      data = tmp.pa,
      verbose = verbose
    )
    create_download_link(data = tmp.pa, logfile = logfile, filename = 'aligned-observed-taxa.csv', linktext = 'Download aligned observed taxa', verbose = verbose)


    # Container matrix.
    n.observed.sites <- dim(tmp.pa)[1]
    # Write to the logs for container matrix
    writelog(
      '\nContainer matrix.',
      logfile = logfile,
      code = "
        n.observed.sites <- dim(tmp.pa)[1]
      ",
      data = n.observed.sites,
      verbose = verbose
    )
    create_download_link(data = n.observed.sites, logfile = logfile, filename = 'container-matrix.csv', linktext = 'Download container matrix', verbose = verbose)


    # get number of reference taxa
    n.reference.taxa <- dim(reference.taxa)[2]
    # Write to the logs for getting number of reference taxa
    writelog(
      '\nGet number of reference taxa',
      logfile = logfile,
      code = "
        n.reference.taxa <- dim(reference.taxa)[2]
      ",
      data = n.reference.taxa,
      verbose = verbose
    )
    create_download_link(data = n.reference.taxa, logfile = logfile, filename = 'reference-taxa.csv', linktext = 'Download reference taxa count', verbose = verbose)


    # Observed Taxa presence absence matrix
    observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa),
                               nrow = n.observed.sites, ncol = n.reference.taxa,
                               dimnames = list(rownames(tmp.pa), names(reference.taxa)))
    # Write to the logs for Observed Taxa presence absence matrix
    writelog(
      '\nObserved Taxa presence absence matrix',
      logfile = logfile,
      code = "
        observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa),
                                   nrow = n.observed.sites, ncol = n.reference.taxa,
                                   dimnames = list(rownames(tmp.pa), names(reference.taxa)))
      ",
      data = observed.taxa.pa,
      verbose = verbose
    )
    create_download_link(data = observed.taxa.pa, logfile = logfile, filename = 'observed-taxa-pa-matrix.csv', linktext = 'Download Observed Taxa PA matrix', verbose = verbose)



    # Fill container with observed data.
    col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])
    for(i in 1:n.reference.taxa) {
      if(!is.na(col.match[i])) observed.taxa.pa[, i] <- tmp.pa[, col.match[i]]
    }
    # Write to the logs for filling container with observed data
    writelog(
      '\nFill container with observed data.',
      logfile = logfile,
      code = "
        col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])
      ",
      data = col.match,
      verbose = verbose
    )
    create_download_link(data = col.match, logfile = logfile, filename = 'container-observed-data.csv', linktext = 'Download container observed data', verbose = verbose)



    # The matrix observed.taxa.pa contains the observed.scores used for O/E.
    return(observed.taxa.pa)

  }
  # Write to the logs for defining the FormatObservedData function
  writelog(
    '\nDefine function - Format Observed Data',
    logfile = logfile,
    code = "
      FormatObservedData <- function() {

        # set up the hyphen log prefix - which hasnt yet worked as i want it to
        hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

        # Align observed (user) data columns with reference data columns. Columns in same
        # order. Observed data may have a different number of taxa (columns) than
        # reference data.

        # Convert observed.taxa to presence/absence (0/1)
        tmp.pa <- observed.taxa
        tmp.pa[tmp.pa > 0] <- 1

        # Align rows using predictor variables.
        tmp.pa <- tmp.pa[row.names(observed.predictors), ]    # !!! is this required???

        # Container matrix.
        n.observed.sites <- dim(tmp.pa)[1]

        # get number of reference taxa
        n.reference.taxa <- dim(reference.taxa)[2]

        # Observed Taxa presence absence matrix
        observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa),
                                   nrow = n.observed.sites, ncol = n.reference.taxa,
                                   dimnames = list(rownames(tmp.pa), names(reference.taxa)))


        # Fill container with observed data.
        col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])
        for(i in 1:n.reference.taxa) {
          if(!is.na(col.match[i])) observed.taxa.pa[, i] <- tmp.pa[, col.match[i]]
        }

        # The matrix observed.taxa.pa contains the observed.scores used for O/E.
        return(observed.taxa.pa)
      }
    ",
    verbose = verbose
  )


  # Call Format observed data function
  observed.data <- FormatObservedData(logfile = logfile, verbose = verbose)
  # Write to the logs for calling FormatObservedData function
  writelog(
    '\nCall FormatObservedData function',
    logfile = logfile,
    code = "
      observed.data <- FormatObservedData()
    ",
    data = observed.data,
    verbose = verbose
  )
  create_download_link(data = observed.data, logfile = logfile, filename = 'formatted-observed-data.csv', linktext = 'Download formatted observed data', verbose = verbose)




  # ----- Define Calculate Expected Data -----
  CalculateExpectedData <- function(
    logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
    verbose = F
  ) {

    # set up the hyphen log prefix - which hasnt yet worked as i want it to
    hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)



    # Calculate probability of sites belonging to groups. Follow RIVPACS assumption
    # of weighting the group probabilities by reference group size. Flags outlier
    # sites, using the chi-squared statistic.


    # In Calcluate Expected Data Function - Definitions.
    n.predictor.variables <- length(predictor.variables)
    group.size <- table(reference.groups)
    n.groups <- length(group.size)

    # Write to the logs for Calculate Expected Data Function - Definitions
    writelog(
      '\nIn Calculate Expected Data Function - Definitions',
      logfile = logfile,
      code = "
        n.predictor.variables <- length(predictor.variables)
        group.size <- table(reference.groups)
        n.groups <- length(group.size)
      ",
      verbose = verbose
    )
    writelog( paste0("\nn.predictor.variables: " , n.predictor.variables), logfile = logfile, verbose = verbose )
    writelog( paste0("\nn.groups: " , n.groups), logfile = logfile, verbose = verbose )
    writelog(
      '\nGroup Size:',
      logfile = logfile,
      data = group.size,
      verbose = verbose
    )



    # Chi-squared values for flagging outlier samples.
    degrees.freedom <- min(c(n.predictor.variables, (n.groups - 1)))
    crit.01 <- qchisq(0.99, df = degrees.freedom)
    crit.05 <- qchisq(0.95, df = degrees.freedom)

    # Write to the logs for Chi-squared values for flagging outlier samples
    writelog(
      '\nDegrees of Freedom and Chi-squared values for flagging outlier samples.',
      logfile = logfile,
      code = "
        degrees.freedom <- min(c(n.predictor.variables, (n.groups - 1)))
        crit.01 <- qchisq(0.99, df = degrees.freedom)
        crit.05 <- qchisq(0.95, df = degrees.freedom)

        print('crit.01')
        print(crit.01)
        print('crit.05')
        print(crit.05)
      ",
      data = degrees.freedom,
      verbose = verbose
    )
    create_download_link(data = degrees.freedom, logfile = logfile, filename = 'degress-freedom.csv', linktext = 'Download Degrees of Freedom values', verbose = verbose)



    # Container for probabilities.
    n.observed.sites.filtered <- dim(observed.predictors)[[1]]
    # Write to the logs for container for probabilities
    writelog(
      '\nContainer for probabilities.',
      logfile = logfile,
      code = "
        n.observed.sites.filtered <- dim(observed.predictors)[[1]]
      ",
      data = n.observed.sites.filtered,
      verbose = verbose
    )
    create_download_link(data = data.frame(n.observed.sites.filtered), logfile = logfile, filename = 'container-for-probabilities.csv', linktext = 'Download container for probabilities', verbose = verbose)



    # Group Probabilities
    group.probabilities <- matrix(rep(0, n.observed.sites.filtered * n.groups),
                                  nrow = n.observed.sites.filtered,
                                  dimnames = list(dimnames(observed.predictors)[[1]],
                                                  dimnames(group.means)[[1]]))
    # Write to the logs for Group Probabilities
    writelog(
      '\nGroup Probabilities',
      logfile = logfile,
      code = "
        group.probabilities <- matrix(rep(0, n.observed.sites.filtered * n.groups),
                                      nrow = n.observed.sites.filtered,
                                      dimnames = list(dimnames(observed.predictors)[[1]],
                                                      dimnames(group.means)[[1]]))
      ",
      data = group.probabilities,
      verbose = verbose
    )
    create_download_link(data = group.probabilities, logfile = logfile, filename = 'group-probabilities-initial.csv', linktext = 'Download group probabilities (initial)', verbose = verbose)





    # Container for outlier flags and minimum distance.
    outlier.flag <- data.frame(outlier.05 = rep(0, n.observed.sites.filtered),
                               outlier.01 = rep(0, n.observed.sites.filtered),
                               min.distance = rep(0, n.observed.sites.filtered),
                               row.names = dimnames(observed.predictors)[[1]])
    # Write to the logs for container for outlier flags and minimum distance
    writelog(
      '\nContainer for outlier flags and minimum distance.',
      logfile = logfile,
      code = "
        outlier.flag <- data.frame(outlier.05 = rep(0, n.observed.sites.filtered),
                                   outlier.01 = rep(0, n.observed.sites.filtered),
                                   min.distance = rep(0, n.observed.sites.filtered),
                                   row.names = dimnames(observed.predictors)[[1]])
      ",
      data = outlier.flag,
      verbose = verbose
    )
    create_download_link(data = outlier.flag, logfile = logfile, filename = 'outlier-flags-and-distance-initial.csv', linktext = 'Download outlier flags and distance (iniital)', verbose = verbose)




    # calculate group membership probabilities for each sample and find outliers.
    for(i in 1:n.observed.sites.filtered) {

      # Squared Mahalanobis distance from current sample to each group mean.
      distances <- mahalanobis(group.means,
                               observed.predictors[i,],
                               reference.cov,
                               inverted = TRUE)

      group.probabilities[i,] <- group.size * exp(-0.5 * distances) # see Clarke et al. (2000)
      group.probabilities[i,] <- group.probabilities[i, ] / sum(group.probabilities[i, ])

      # Outlier criteria is minimum distance.
      outlier.flag$min.distance[i] <- min(distances)

      # Check for outliers. Each sample is either a pass (0) or fail (1).
      if(outlier.flag$min.distance[i] > crit.05) outlier.flag[i, "outlier.05"] <- 1
      if(outlier.flag$min.distance[i] > crit.01) outlier.flag[i, "outlier.01"] <- 1

    }
    # Write to the logs for calculating group membership probabilities and finding outliers
    writelog(
      '\nCalculate group membership probabilities for each sample and find outliers.',
      logfile = logfile,
      code = "
        for(i in 1:n.observed.sites.filtered) {

          # Squared Mahalanobis distance from current sample to each group mean.
          distances <- mahalanobis(group.means,
                                   observed.predictors[i,],
                                   reference.cov,
                                   inverted = TRUE)

          group.probabilities[i,] <- group.size * exp(-0.5 * distances) # see Clarke et al. (2000)
          group.probabilities[i,] <- group.probabilities[i, ] / sum(group.probabilities[i, ])

          # Outlier criteria is minimum distance.
          outlier.flag$min.distance[i] <- min(distances)

          # Check for outliers. Each sample is either a pass (0) or fail (1).
          if(outlier.flag$min.distance[i] > crit.05) outlier.flag[i, 'outlier.05'] <- 1
          if(outlier.flag$min.distance[i] > crit.01) outlier.flag[i, 'outlier.01'] <- 1

        }
      ",
      verbose = verbose
    )
    create_download_link(data = group.probabilities, logfile = logfile, filename = 'group-probabilities.csv', linktext = 'Download group probabilities', verbose = verbose)
    create_download_link(data = outlier.flag, logfile = logfile, filename = 'outlier-flag.csv', linktext = 'Download outlier flag', verbose = verbose)




    # Occurrence frequencies of all taxa in the reference groups.
    freq.in.group <- apply(reference.taxa, 2,
                           function(x){tapply(x, reference.groups, function(y){sum(y) / length(y)})})
    # Write to the logs for occurrence frequencies of all taxa in the reference groups
    writelog(
      '\nOccurrence frequencies of all taxa in the reference groups.',
      logfile = logfile,
      code = "
        freq.in.group <- apply(reference.taxa, 2,
                               function(x) { tapply(x, reference.groups, function(y) { sum(y) / length(y) }) })
      ",
      data = freq.in.group,
      verbose = verbose
    )
    create_download_link(data = freq.in.group, logfile = logfile, filename = 'occurrence-frequencies.csv', linktext = 'Download occurrence frequencies', verbose = verbose)


    # Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).
    predicted.prob.all <- group.probabilities %*% freq.in.group
    # Write to the logs for Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4)
    writelog(
      '\nMatrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).',
      logfile = logfile,
      code = "
        predicted.prob.all <- group.probabilities %*% freq.in.group
      ",
      data = predicted.prob.all,
      verbose = verbose
    )
    create_download_link(data = predicted.prob.all, logfile = logfile, filename = 'predicted-prob-all.csv', linktext = 'Download predicted probabilities', verbose = verbose)


    # predicted.prob.all are the predicted (expected) probabilites.
    expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)
    # Write to the logs for predicted probabilities
    writelog(
      '\npredicted.prob.all are the predicted (expected) probabilities.',
      logfile = logfile,
      code = "
        expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)
      ",
      verbose = verbose
    )
    create_download_link(data = predicted.prob.all, logfile = logfile, filename = 'expected-predicted-prob.csv', linktext = 'Download predicted probabilities', verbose = verbose)
    create_download_link(data = outlier.flag, logfile = logfile, filename = 'expected-outliers.csv', linktext = 'Download outliers', verbose = verbose)
    create_download_link(data = data.frame(n.observed.sites.filtered), logfile = logfile, filename = 'expected-n-sites.csv', linktext = 'Download number of observed sites', verbose = verbose)


    writelog('\nEND: Calculate Expected Data function.\n', logfile = logfile, verbose = verbose)

    return(expected.data)

  }
  # Write to the logs for defining CalculateExpectedData function
  writelog(
    '\nDefine CalculateExpectedData function',
    logfile = logfile,
    code = "
      CalculateExpectedData <- function() {

        # set up the hyphen log prefix - which hasn't yet worked as I want it to
        hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

        # Calculate probability of sites belonging to groups. Follow RIVPACS assumption
        # of weighting the group probabilities by reference group size. Flags outlier
        # sites, using the chi-squared statistic.

        # In Calculate Expected Data Function - Definitions.
        n.predictor.variables <- length(predictor.variables)
        group.size <- table(reference.groups)
        n.groups <- length(group.size)

        # Chi-squared values for flagging outlier samples.
        degrees.freedom <- min(c(n.predictor.variables, (n.groups - 1)))
        crit.01 <- qchisq(0.99, df = degrees.freedom)
        crit.05 <- qchisq(0.95, df = degrees.freedom)

        # Container for probabilities.
        n.observed.sites.filtered <- dim(observed.predictors)[[1]]

        # Group Probabilities
        group.probabilities <- matrix(rep(0, n.observed.sites.filtered * n.groups),
                                      nrow = n.observed.sites.filtered,
                                      dimnames = list(dimnames(observed.predictors)[[1]],
                                                      dimnames(group.means)[[1]]))

        # Container for outlier flags and minimum distance.
        outlier.flag <- data.frame(outlier.05 = rep(0, n.observed.sites.filtered),
                                   outlier.01 = rep(0, n.observed.sites.filtered),
                                   min.distance = rep(0, n.observed.sites.filtered),
                                   row.names = dimnames(observed.predictors)[[1]])

        # Calculate group membership probabilities for each sample and find outliers.
        for(i in 1:n.observed.sites.filtered) {

          # Squared Mahalanobis distance from current sample to each group mean.
          distances <- mahalanobis(group.means,
                                   observed.predictors[i,],
                                   reference.cov,
                                   inverted = TRUE)

          group.probabilities[i,] <- group.size * exp(-0.5 * distances) # see Clarke et al. (2000)
          group.probabilities[i,] <- group.probabilities[i, ] / sum(group.probabilities[i, ])

          # Outlier criteria is minimum distance.
          outlier.flag$min.distance[i] <- min(distances)

          # Check for outliers. Each sample is either a pass (0) or fail (1).
          if(outlier.flag$min.distance[i] > crit.05) outlier.flag[i, 'outlier.05'] <- 1
          if(outlier.flag$min.distance[i] > crit.01) outlier.flag[i, 'outlier.01'] <- 1

        }

        # Occurrence frequencies of all taxa in the reference groups.
        freq.in.group <- apply(reference.taxa, 2,
                               function(x){tapply(x, reference.groups, function(y){sum(y) / length(y)})})

        # Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).
        predicted.prob.all <- group.probabilities %*% freq.in.group

        # predicted.prob.all are the predicted (expected) probabilities.
        expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)

        return(expected.data)

      }
    ",
    verbose = verbose
  )


  # Call CalculateExpectedData function
  expected.data <- CalculateExpectedData(logfile = logfile, verbose = verbose)
  # Write to the logs for calling CalculateExpectedData function
  writelog(
    '\nCall CalculateExpectedData function',
    logfile = logfile,
    code = "
      expected.data <- CalculateExpectedData()
    ",
    verbose = verbose
  )
  create_download_link(data = expected.data$predicted, logfile = logfile, filename = 'expected-predicted-prob.csv', linktext = 'Download predicted probabilities (returned from CalculateExpectedData function)', verbose = verbose)
  create_download_link(data = expected.data$outliers, logfile = logfile, filename = 'expected-outliers.csv', linktext = 'Download outliers (returned from CalculateExpectedData function)', verbose = verbose)
  create_download_link(data = data.frame(expected.data$n), logfile = logfile, filename = 'expected-n-sites.csv', linktext = 'Download number of observed sites (returned from CalculateExpectedData function)', verbose = verbose)




  # Define Calculate Scores Function -----
  CalculateScores <- function(
    logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
    verbose = F
  ) {

    # set up the hyphen log prefix - which hasnt yet worked as i want it to
    hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

    writelog('\nBEGIN: Calculate Scores function.\n', logfile = logfile, verbose = verbose)

    # Bray-Curtis dissimilarity
    observed.score <- vector(mode = "numeric", length = expected.data$n)
    expected.score <- vector(mode = "numeric", length = expected.data$n)
    BC <- vector(mode = "numeric", length = expected.data$n) # Bray-Curtis dissimilarity
    # Write to the logs for Bray-Curtis dissimilarity
    writelog(
      '\nBray-Curtis dissimilarity',
      logfile = logfile,
      code = "
        observed.score <- vector(mode = 'numeric', length = expected.data$n)
        expected.score <- vector(mode = 'numeric', length = expected.data$n)
        BC <- vector(mode = 'numeric', length = expected.data$n) # Bray-Curtis dissimilarity
      ",
      data = data.frame(observed.score, expected.score, BC),
      verbose = verbose
    )

    for(i in 1:expected.data$n) {
      tryCatch(
        {
          # predicted probabilities for current sample
          predicted.prob <- expected.data$predicted[i, ]

          # Write to the logs for: predicted probabilities for current sample
          writelog(
            '\npredicted probabilities for current sample',
            logfile = logfile,
            code = "
              predicted.prob <- expected.data$predicted[i, ]
              print('predicted.prob')
              print(predicted.prob)
            ",
            verbose = verbose
          )

          # subset of taxa with probabilities >= Pcutoff
          taxa.subset <- names(predicted.prob)[predicted.prob >= Pcutoff]
          # Write to the logs for: subset of taxa with probabilities >= Pcutoff
          writelog(
            '\nsubset of taxa with probabilities >= Pcutoff',
            logfile = logfile,
            code = "
              taxa.subset <- names(predicted.prob)[predicted.prob >= Pcutoff]
              print('taxa.subset')
              print(taxa.subset)
            ",
            verbose = verbose
          )


          # probabilites for subset of included taxa
          expected.prob <- predicted.prob[taxa.subset]
          # Write to the logs for: probabilities for subset of included taxa
          writelog(
            '\nprobabilities for subset of included taxa',
            logfile = logfile,
            code = "
              expected.prob <- predicted.prob[taxa.subset]
              print('expected.prob')
              print(expected.prob)
            ",
            verbose = verbose
          )


          # observed presence/absence for those taxa
          observed.pa <- observed.data[i, taxa.subset]
          # Write to the logs for: observed presence/absence for those taxa
          writelog(
            '\nobserved presence/absence for those taxa',
            logfile = logfile,
            code = "
              observed.pa <- observed.data[i, taxa.subset]
              print('observed.pa')
              print(observed.pa)
            ",
            verbose = verbose
          )


          # observed richness (O)
          observed.score[i] <- sum(observed.pa)
          # Write to the logs for: observed richness (O)
          writelog(
            '\nobserved richness (O)',
            logfile = logfile,
            code = "
              observed.score[i] <- sum(observed.pa)
              print('observed.score[i]')
              print(observed.score[i])
            ",
            verbose = verbose
          )


          # expected richness (E)
          expected.score[i] <- sum(expected.prob)
          # Write to the logs for: expected richness (E)
          writelog(
            '\nexpected richness (E)',
            logfile = logfile,
            code = "
              expected.score[i] <- sum(expected.prob)
              print('expected.score[i]')
              print(expected.score[i])
            ",
            verbose = verbose
          )


          # BC value
          BC[i] <- sum(abs(observed.pa - expected.prob)) /
            (observed.score[i] + expected.score[i])
          # Write to the logs for: BC value
          writelog(
            '\nBC value',
            logfile = logfile,
            code = "
              BC[i] <- sum(abs(observed.pa - expected.prob)) / (observed.score[i] + expected.score[i])
              print('BC[i]')
              print(BC[i])
            ",
            verbose = verbose
          )

        },
        error = function(e) {

          # observed richness (O)
          observed.score[i] <- NA_real_
          # Write to the logs for: observed richness (O)
          writelog(
            '\nobserved richness (O)',
            logfile = logfile,
            code = "
              observed.score[i] <- NA_real_
              print('observed.score[i]')
              print(observed.score[i])
            ",
            verbose = verbose
          )


          # expected richness (E)
          expected.score[i] <- NA_real_
          # Write to the logs for: expected richness (E)
          writelog(
            '\nexpected richness (E)',
            logfile = logfile,
            code = "
              expected.score[i] <- NA_real_
              print('expected.score[i]')
              print(expected.score[i])
            ",
            verbose = verbose
          )


          # BC value
          BC[i] <- NA_real_
          # Write to the logs for: BC value
          writelog(
            '\n# BC value',
            logfile = logfile,
            code = "
              BC[i] <- NA_real_
              print('BC[i]')
              print(BC[i])
            ",
            verbose = verbose
          )

        }
      )
    }


    # Get the stats dataframe
    O.over.E <- observed.score/expected.score
    stats <- data.frame(stations = row.names(observed.predictors),
                        O = observed.score,
                        E = round(expected.score, digits = 4),
                        O.over.E = round(O.over.E, digits = 4))
    # Write to the logs for getting the stats dataframe
    writelog(
      '\nGet the stats dataframe',
      logfile = logfile,
      code = "
        O.over.E <- observed.score / expected.score
        stats <- data.frame(
          stations = row.names(observed.predictors),
          O = observed.score,
          E = round(expected.score, digits = 4),
          O.over.E = round(O.over.E, digits = 4)
        )
      ",
      data = stats,
      verbose = verbose
    )
    create_download_link(data = stats, logfile = logfile, filename = 'stats-dataframe.csv', linktext = 'Download stats dataframe', verbose = verbose)



    # create outlier columns on that stats dataframe
    stats$outlier.05 <- expected.data$outliers$outlier.05
    stats$outlier.01 <- expected.data$outliers$outlier.01
    # Write to the logs for creating outlier columns on the stats dataframe
    writelog(
      '\nCreate outlier columns on the stats dataframe',
      logfile = logfile,
      code = "
        stats$outlier.05 <- expected.data$outliers$outlier.05
        stats$outlier.01 <- expected.data$outliers$outlier.01
      ",
      data = stats,
      verbose = verbose
    )
    create_download_link(data = stats, logfile = logfile, filename = 'stats-dataframe-with-outliers.csv', linktext = 'Download stats dataframe with outliers', verbose = verbose)


    # Convert to "PASS" or "FAIL"
    stats$outlier.05[stats$outlier.05 == 0] <- "PASS"
    stats$outlier.05[stats$outlier.05 == 1] <- "FAIL"
    stats$outlier.01[stats$outlier.01 == 0] <- "PASS"
    stats$outlier.01[stats$outlier.01 == 1] <- "FAIL"
    # Write to the logs for converting outlier columns to "PASS" or "FAIL"
    writelog(
      '\nConvert outlier columns to "PASS" or "FAIL"',
      logfile = logfile,
      code = "
        stats$outlier.05[stats$outlier.05 == 0] <- 'PASS'
        stats$outlier.05[stats$outlier.05 == 1] <- 'FAIL'
        stats$outlier.01[stats$outlier.01 == 0] <- 'PASS'
        stats$outlier.01[stats$outlier.01 == 1] <- 'FAIL'
      ",
      data = stats,
      verbose = verbose
    )
    create_download_link(data = stats, logfile = logfile, filename = 'stats-dataframe-pass-fail.csv', linktext = 'Download stats dataframe with PASS/FAIL outliers', verbose = verbose)


    writelog('\n## END: Calculate Scores function.\n', logfile = logfile, verbose = verbose)

    return(stats)

  }

  # Final Results
  results <- list(oe.table = CalculateScores(logfile = logfile, verbose = verbose),
                  observed = observed.data,
                  predicted = expected.data$predicted,
                  Pcutoff = Pcutoff,
                  region = "scb")
  # Write to the logs for final results
  writelog(
    '\nFinal Results',
    logfile = logfile,
    code = "
      results <- list(
        oe.table = CalculateScores(logfile = logfile, verbose = verbose),
        observed = observed.data,
        predicted = expected.data$predicted,
        Pcutoff = Pcutoff,
        region = 'scb'
      )
    ",
    data = results,
    verbose = verbose
  )
  create_download_link(data = results$oe.table, logfile = logfile, filename = 'final-results-oe-table.csv', linktext = 'Download final results OE table', verbose = verbose)
  create_download_link(data = results$observed, logfile = logfile, filename = 'final-results-observed.csv', linktext = 'Download final results observed', verbose = verbose)
  create_download_link(data = results$predicted, logfile = logfile, filename = 'final-results-predicted.csv', linktext = 'Download final results predicted', verbose = verbose)
  create_download_link(data = data.frame(results$Pcutoff), logfile = logfile, filename = 'final-results-pcutoff.csv', linktext = 'Download final results Pcutoff', verbose = verbose)
  create_download_link(data = data.frame(results$region), logfile = logfile, filename = 'final-results-region.csv', linktext = 'Download final results region', verbose = verbose)


  writelog('\n### END: So Cal RIVPACS function.\n', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  return(results)

}



