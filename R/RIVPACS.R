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

  writelog('\nBEGIN: RIVPACS function.\n', logfile = logfile, verbose = verbose)

  writelog('*** DATA *** Input to RIVPACS function - RIVPACS-step0.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(benthic_data, logfile = file.path(dirname(logfile), 'RIVPACS-step0.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # Split to SoCal and SFBay.
  ## We are only working with SoCal data so we don't need to do this!


  # SCB Predictors - needs to be logged
  scb.predictors <- data.frame(Latitude = benthic_data$Latitude,
                               Longitude = benthic_data$Longitude,
                               SampleDepth = benthic_data$SampleDepth) %>%
    dplyr::distinct()

  writelog('*** DATA *** RIVPACS-SCBpredictors-initial.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(scb.predictors, logfile = file.path(dirname(logfile), 'RIVPACS-SCBpredictors-initial.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # data prep step 1 - rename taxa to Taxon
  writelog("Renamed Taxa to Taxon in the input data (benthic_data)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  benthic_data <- benthic_data %>% dplyr::rename(Taxa = Taxon)

  # data prep step 2 - get distinct records on StationID, Latitude, Longitude, SampleDepth - also set row names to StationID
  writelog("get distinct records on StationID, Latitude, Longitude, SampleDepth - also set row names to StationID", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  scb.taxa <- benthic_data %>% dplyr::select(StationID, Latitude, Longitude, SampleDepth) %>%
    dplyr::distinct()
  writelog('*** DATA *** get distinct records on StationID, Latitude, Longitude, SampleDepth - also set row names to StationID - RIVPACS-SoCalPrep-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(scb.taxa, logfile = file.path(dirname(logfile), 'RIVPACS-SoCalPrep-step1.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # *Log scb.taxa

  # data prep step 3 - set up the scb predictors
  row.names(scb.predictors) <- scb.taxa$StationID
  scb.predictors <- as.matrix(scb.predictors)
  writelog('*** DATA *** SCB predictors with stationids - RIVPACS-SCBpredictors-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(scb.predictors, logfile = file.path(dirname(logfile), 'RIVPACS-SCBpredictors-step1.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix, include.row.names = T)


  # data prep step 4 - Filter to replicate one and get distinct values on StationID Taxa and Abundance
  writelog("Filter to replicate one and get distinct values on StationID Taxa and Abundance", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  scb.taxa <- benthic_data %>%
    dplyr::filter(Replicate == 1) %>%
    dplyr::select(StationID, Taxa, Abundance) %>%
    dplyr::distinct()
  writelog('*** DATA *** Filter to replicate one and get distinct values on StationID Taxa and Abundance - RIVPACS-SoCalPrep-step2.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(scb.taxa, logfile = file.path(dirname(logfile), 'RIVPACS-SoCalPrep-step2.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # Data prep step 5 - remove certain special characters from taxa name
  writelog("Filter to replicate one and get distinct values on StationID Taxa and Abundance", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  scb.taxa$Taxa <- gsub(" ", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub("(", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub(")", "_", scb.taxa$Taxa, fixed = TRUE)
  writelog('*** DATA *** Filter to replicate one and get distinct values on StationID Taxa and Abundance - RIVPACS-SoCalPrep-step3.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(scb.taxa, logfile = file.path(dirname(logfile), 'RIVPACS-SoCalPrep-step3.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # Data prep step 6 - pivot the data out wide and make it a data.frame
  writelog("pivot the data out wide and make it a data.frame", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  scb.taxa <- scb.taxa %>%
    tidyr::pivot_wider(id_cols = "StationID", names_from = "Taxa",
                       values_from = "Abundance", values_fn = list(Abundance = list))
  scb.taxa <- as.data.frame(scb.taxa)

  # Log the pivoting action
  writelog('*** DATA *** pivot the data out wide and make it a data.frame - RIVPACS.SoCalPrep.step4.RData', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  # Save scb.taxa as an RData file
  rdata_file_path <- file.path(dirname(logfile), 'RIVPACS.SoCalPrep.step4.RData')
  RIVPACS.SoCalPrep.step4 <- scb.taxa
  save(RIVPACS.SoCalPrep.step4, file = rdata_file_path)
  # Log the saving action
  writelog(paste('Saved scb.taxa as RData file:', rdata_file_path), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste('** NOTE ** I needed to save it as an RData file since there were objects in the dataframe that were not able to be written to a csv', rdata_file_path), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  # Data prep step 7 - ...
  writelog("Take the last column of scb taxa", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  scb.taxa <- scb.taxa[, -1]

  # data prep step 8 - remove Abundance. from column names
  colnames(scb.taxa) <- gsub("Abundance.", "", colnames(scb.taxa))
  writelog('RIVPACS-SoCalPrep-step6', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  # writelog(paste(c('RIVPACS-SoCalPrep-step6', scb.taxa), collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('** NOTE ** RIVPACS-SoCalPrep-step6 was a very long list that cluttered the log - omitting logging for now', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)


  # data prep step 9 - Replace NAs with zero.
  scb.taxa[scb.taxa == "NULL"] <- 0
  scb.taxa = as.data.frame(lapply(scb.taxa, as.numeric))
  row.names(scb.taxa) <- row.names(scb.predictors)
  writelog('*** DATA *** Replace NAs with zero. - RIVPACS-SoCalPrep-step7.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(scb.taxa, logfile = file.path(dirname(logfile), 'RIVPACS-SoCalPrep-step7.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # RIVPACS calculations. By default the functions use the example user data.
  writelog('About to call SoCal RIVPACS within RIVPACS function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  socal <- SoCalRivpacs(observed.predictors = scb.predictors, observed.taxa = scb.taxa, logfile = logfile, verbose = verbose)
  writelog('Done calling SoCal RIVPACS within RIVPACS function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  # the stations column of the oe table dataframe was being returned as a factor. Need to make that a character
  writelog('Transform factor columns into characters', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  socal$oe.table <- socal$oe.table %>%
    mutate_if(is.factor,as.character)
  writelog('*** DATA *** SoCal OE Table - after making factor columns into characters - RIVPACS-socal-oe-factors-to-characters.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(socal$oe.table, logfile = file.path(dirname(logfile), 'RIVPACS-socal-oe-factors-to-characters.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)




  # **Please Log socal$oe.table**

  # Get the distinct values in benthic data based on StationID Replicate SampleDate Stratum
  benthic_data <- benthic_data %>%
    dplyr::select(StationID, Replicate, SampleDate, Stratum) %>%
    dplyr::distinct()

  writelog('*** DATA *** benthic data distinct on StationID, Replicate, SampleDate, Stratum - RIVPACS-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(benthic_data, logfile = file.path(dirname(logfile), 'RIVPACS-step1.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # Calculate RIVPACS Scores

  # Riv step 0 - select stations and Observed/Expected
  riv0 <- socal$oe.table %>%
    dplyr::select(stations, O.over.E)

  writelog('*** DATA *** select stations and Observed/Expected - RIVPACS-step2.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(riv0, logfile = file.path(dirname(logfile), 'RIVPACS-step2.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # Riv step 1 - join with benthic data
  riv1 <- riv0 %>%
    dplyr::rename(StationID = stations, Score = O.over.E) %>%
    dplyr::full_join(benthic_data) %>%
    dplyr::mutate(Index = "RIVPACS")

  writelog('*** DATA *** join with benthic data - RIVPACS-step3.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(riv0, logfile = file.path(dirname(logfile), 'RIVPACS-step3.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  # Get the scores based on the thresholds
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


  writelog('*** DATA *** Get the scores based on the thresholds - FINAL RIVPACS OUTPUT - RIVPACS-final.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(riv0, logfile = file.path(dirname(logfile), 'RIVPACS-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  writelog('\nEND: RIVPACS function.\n', logfile = logfile, verbose = verbose)


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

  writelog('\nBEGIN: So Cal RIVPACS function.\n', logfile = logfile, verbose = verbose)

  # set up the hyphen log prefix - which hasnt yet worked as i want it to
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("Begin Logging input to So Cal Rivpacs", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog(paste(c('Reference Groups: ', reference.groups), collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog('So Cal RIVPACS observed predictors - observed.predictors.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(observed.predictors, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-observed.predictors.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix, include.row.names = T)

  writelog('So Cal RIVPACS observed predictors - reference.taxa.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(reference.taxa, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-reference.taxa.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix, include.row.names = T)


  writelog('So Cal RIVPACS observed predictors - group.means.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(group.means, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-group.means.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix, include.row.names = T)

  writelog('So Cal RIVPACS observed predictors - reference.cov.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(reference.cov, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-reference.cov.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix, include.row.names = T)

  writelog('So Cal RIVPACS observed predictors - observed.taxa.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(observed.taxa, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-observed.taxa.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix, include.row.names = T)



  # Pcutoff is the probability cutoff

  # Names of predictor variables.
  predictor.variables <- c("Latitude", "Longitude", "SampleDepth")

  writelog('These are the predictor variables in socal rivpacs', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste(predictor.variables, collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  # ----- Format Observed Data -----

  FormatObservedData <- function(
    logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
    verbose = F
  ) {

    # set up the hyphen log prefix - which hasnt yet worked as i want it to
    hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

    writelog('\nBEGIN: FormatObservedData function (within SoCal RIVPACS.\n', logfile = logfile, verbose = verbose)

    # Align observed (user) data columns with reference data columns. Columns in same
    # order. Observed data may have a different number of taxa (columns) than
    # reference data.

    # Convert observed.taxa to presence/absence (0/1)
    tmp.pa <- observed.taxa

    writelog('Observed Taxa - FormatObservedData-observed.taxa.csv', logfile = logfile, verbose = verbose)
    writelog(observed.taxa, logfile = file.path(dirname(logfile), 'FormatObservedData-observed.taxa.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    #tmp.pa <- lapply(tmp.pa, as.numeric)
    tmp.pa[tmp.pa > 0] <- 1

    writelog('Observed Taxa converted to Presence Absence - FormatObservedData-observed.taxa.PA-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(tmp.pa, logfile = file.path(dirname(logfile), 'FormatObservedData-observed.taxa.PA-step1.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    # Align rows using predictor variables.
    tmp.pa <- tmp.pa[row.names(observed.predictors), ]    # !!! is this required???


    # Container matrix.
    n.observed.sites <- dim(tmp.pa)[1]

    writelog(paste('Number of observed sites: ', n.observed.sites), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    n.reference.taxa <- dim(reference.taxa)[2]

    writelog(paste('Number of referenced taxa: ', n.reference.taxa), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa),
                               nrow = n.observed.sites, ncol = n.reference.taxa,
                               dimnames = list(rownames(tmp.pa), names(reference.taxa)))

    writelog('Observed Taxa Presence Absence - step 2 - FormatObservedData-observed.taxa.PA-step2.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(observed.taxa.pa, logfile = file.path(dirname(logfile), 'FormatObservedData-observed.taxa.PA-step2.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    writelog('match the observed taxa PA (step 2) with the taxa PA from step 1', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    # Fill container with observed data.
    col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])
    #tmp.pa <- as.data.frame(lapply(tmp.pa, as.numeric))


    writelog('replace columns of observed.taxa.pa with those of tmp.pa if the column names match\n', logfile = logfile, verbose = verbose)
    for(i in 1:n.reference.taxa) {
      if(!is.na(col.match[i])) observed.taxa.pa[, i] <- tmp.pa[, col.match[i]]
    }

    writelog('Final Observed Taxa PA - FormatObservedData-observed.taxa.PA-final.csv', logfile = logfile, verbose = verbose)
    writelog(observed.taxa.pa, logfile = file.path(dirname(logfile), 'FormatObservedData-observed.taxa.PA-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

    writelog('\nEND: FormatObservedData function (within SoCal RIVPACS.\n', logfile = logfile, verbose = verbose)

    # The matrix observed.taxa.pa contains the observed.scores used for O/E.
    return(observed.taxa.pa)

  }

  writelog('About to Call the Final Observed Data function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  observed.data <- FormatObservedData(logfile = logfile, verbose = verbose)

  writelog('Done Calling the Final Observed Data function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog('Output from the function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(observed.data, logfile = file.path(dirname(logfile), 'FormatObservedData-observed.taxa-function.output.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)




  # ----- Calculate Expected Data -----

  CalculateExpectedData <- function(
    logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
    verbose = F
  ) {

    # set up the hyphen log prefix - which hasnt yet worked as i want it to
    hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

    writelog('\nBEGIN: Calculate Expected Data function.\n', logfile = logfile, verbose = verbose)

    # Calculate probability of sites belonging to groups. Follow RIVPACS assumption
    # of weighting the group probabilities by reference group size. Flags outlier
    # sites, using the chi-squared statistic.


    # Definitions.
    n.predictor.variables <- length(predictor.variables)
    group.size <- table(reference.groups)
    n.groups <- length(group.size)

    writelog( paste( c('N predictor variables', n.predictor.variables), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog('*** DATA *** Calculate Expected Data - Reference groups - CalculateExpectedData-reference.groups.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(group.size, logfile = file.path(dirname(logfile), 'CalculateExpectedData-reference.groups.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)
    writelog( paste( c('N Groups', n.predictor.variables), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)


    # Chi-squared values for flagging outlier samples.
    degrees.freedom <- min(c(n.predictor.variables, (n.groups - 1)))
    crit.01 <- qchisq(0.99, df = degrees.freedom)
    crit.05 <- qchisq(0.95, df = degrees.freedom)

    writelog( paste( c('degrees.freedom', degrees.freedom), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog( paste( c('crit.01', crit.01), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog( paste( c('crit.05', crit.05), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)


    # Container for probabilities.
    n.observed.sites.filtered <- dim(observed.predictors)[[1]]
    writelog( paste( c('n.observed.sites.filtered', n.observed.sites.filtered), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)


    group.probabilities <- matrix(rep(0, n.observed.sites.filtered * n.groups),
                                  nrow = n.observed.sites.filtered,
                                  dimnames = list(dimnames(observed.predictors)[[1]],
                                                  dimnames(group.means)[[1]]))


    writelog('*** DATA *** Group Probabilities - CalculateExpectedData-group-prob.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(group.probabilities, logfile = file.path(dirname(logfile), 'CalculateExpectedData-group-prob.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    # Container for outlier flags and minimum distance.
    outlier.flag <- data.frame(outlier.05 = rep(0, n.observed.sites.filtered),
                               outlier.01 = rep(0, n.observed.sites.filtered),
                               min.distance = rep(0, n.observed.sites.filtered),
                               row.names = dimnames(observed.predictors)[[1]])

    writelog('*** DATA *** Ouitlier flags - CalculateExpectedData-outlier.flags.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(outlier.flag, logfile = file.path(dirname(logfile), 'CalculateExpectedData-outlier.flags.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    # calculate group membership probabilities for each sample and find outliers.
    writelog('Calculate group membership probabilities for each sample and find outliers.', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog('Get mahalanobis distance from each sample to each group mean', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
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

    writelog('*** DATA *** Ouitlier flags AFTER distance calculation - CalculateExpectedData-outlier.flags-after-dist-calc.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(outlier.flag, logfile = file.path(dirname(logfile), 'CalculateExpectedData-outlier.flags-after-dist-calc.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)



    # Occurrence frequencies of all taxa in the reference groups.
    writelog('Get Occurrence frequencies of all taxa in the reference groups.', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    freq.in.group <- apply(reference.taxa, 2,
                           function(x){tapply(x, reference.groups, function(y){sum(y) / length(y)})})

    # Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).
    writelog('Get Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    predicted.prob.all <- group.probabilities %*% freq.in.group

    writelog('*** DATA *** Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4). - CalculateExpectedData-predicted.prob.all.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(predicted.prob.all, logfile = file.path(dirname(logfile), 'CalculateExpectedData-predicted.prob.all.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    # predicted.prob.all are the predicted (expected) probabilites.
    expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)

    writelog("Final output is a list of items which have already been written out:", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog('-- CalculateExpectedData-predicted.prob.all.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog('-- CalculateExpectedData-outlier.flags-after-dist-calc.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog('-- Also n.observed.sites.filtered which was printed earlier, but here it is again', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog( paste( c('n.observed.sites.filtered', n.observed.sites.filtered), collapse = ': '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)


    writelog('\nEND: Calculate Expected Data function.\n', logfile = logfile, verbose = verbose)

    return(expected.data)

  }

  writelog('About to Call the Final Expected Data function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  expected.data <- CalculateExpectedData(logfile = logfile, verbose = verbose)

  writelog('Done Calling the Final Expected Data function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog('Output from the CalculateExpectedData function:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(expected.data, logfile = file.path(dirname(logfile), 'CalculateExpectedData-expected.data.output.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)





  # ----- Calculate Scores -----

  CalculateScores <- function(
    logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ),
    verbose = F
  ) {

    # set up the hyphen log prefix - which hasnt yet worked as i want it to
    hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

    writelog('\nBEGIN: Calculate Scores function.\n', logfile = logfile, verbose = verbose)

    observed.score <- vector(mode = "numeric", length = expected.data$n)
    expected.score <- vector(mode = "numeric", length = expected.data$n)
    BC <- vector(mode = "numeric", length = expected.data$n) # Bray-Curtis dissimilarity

    writelog("Calculating Observed, Expected, and Bray Curtis dissimilarity", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    for(i in 1:expected.data$n) {
      tryCatch(
        {
          predicted.prob <- expected.data$predicted[i, ] # predicted probabilities for current sample

          taxa.subset <- names(predicted.prob)[predicted.prob >= Pcutoff]  # subset of taxa with probabilities >= Pcutoff

          expected.prob <- predicted.prob[taxa.subset] # probabilites for subset of included taxa

          observed.pa <- observed.data[i, taxa.subset] # observed presence/absence for those taxa

          observed.score[i] <- sum(observed.pa) # observed richness (O)
          expected.score[i] <- sum(expected.prob) # expected richness (E)
          BC[i] <- sum(abs(observed.pa - expected.prob)) /
            (observed.score[i] + expected.score[i]) # BC value
        },
        error = function(e) {
          observed.score[i] <- NA_real_ # observed richness (O)
          expected.score[i] <- NA_real_ # expected richness (E)
          BC[i] <- NA_real_ # BC value
        }
      )


    }

    O.over.E <- observed.score/expected.score

    writelog('Observed', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(paste(observed.score, collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    writelog('Expected', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(paste(expected.score, collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    writelog('BC Dissimilarity', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(paste(BC, collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    writelog('O.over.E before any rounding', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(paste(O.over.E, collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    stats <- data.frame(stations = row.names(observed.predictors),
                        O = observed.score,
                        E = round(expected.score, digits = 4),
                        O.over.E = round(O.over.E, digits = 4))

    writelog('All the above in a dataframe - rounding applied to Expected and Observed/Expected - SOCAL-RIVPACS-O.over.E-initial.csv:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(stats, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-O.over.E-initial.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    #     stats <- data.frame(stations = row.names(observed.predictors),
    #                         O = observed.score,
    #                         E = expected.score,
    #                         O.over.E = O.over.E)


    writelog('Set outlier columns and use to determine pass or fail', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    stats$outlier.05 <- expected.data$outliers$outlier.05
    stats$outlier.01 <- expected.data$outliers$outlier.01

    # Convert to "PASS" or "FAIL"
    stats$outlier.05[stats$outlier.05 == 0] <- "PASS"
    stats$outlier.05[stats$outlier.05 == 1] <- "FAIL"

    stats$outlier.01[stats$outlier.01 == 0] <- "PASS"
    stats$outlier.01[stats$outlier.01 == 1] <- "FAIL"

    #   mean.O.over.E <- mean(OE.stats$O.over.E)
    #   stdev.O.over.E <- sqrt(var(OE.stats$O.over.E))

    writelog('Final So Cal RIVPACS df - SOCAL-RIVPACS-O.over.E-final.csv:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(stats, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-O.over.E-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


    writelog('\nEND: Calculate Scores function.\n', logfile = logfile, verbose = verbose)

    return(stats)

  }

  results <- list(oe.table = CalculateScores(logfile = logfile, verbose = verbose),
                  observed = observed.data,
                  predicted = expected.data$predicted,
                  Pcutoff = Pcutoff,
                  region = "scb")

  writelog("SoCal RIVPACS output - apart from the O/E Table, which was already written out upon calling Calculate Scores (SOCAL-RIVPACS-O.over.E-final.csv)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('observed - SOCAL-RIVPACS-observed.data-final.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(results$observed, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-observed.data-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  writelog('SOCAL RIVPACS predicted', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(results$predicted, logfile = file.path(dirname(logfile), 'SOCAL-RIVPACS-predicted-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  writelog( paste('P Cutoff', results$Pcutoff, collapse = ': ') , logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog( paste('Region', results$region, collapse = ': ') , logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog('\nEND: So Cal RIVPACS function.\n', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  return(results)

}



