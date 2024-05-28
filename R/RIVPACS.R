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
RIVPACS <- function(benthic_data, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T){

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  # Split to SoCal and SFBay.
  ## We are only working with SoCal data so we don't need to do this!


  # SCB Predictors - needs to be logged
  scb.predictors <- data.frame(Latitude = benthic_data$Latitude,
                               Longitude = benthic_data$Longitude,
                               SampleDepth = benthic_data$SampleDepth) %>%
    dplyr::distinct()

  # data prep step 1 - rename taxa to Taxon
  benthic_data <- benthic_data %>% dplyr::rename(Taxa = Taxon)

  # data prep step 2 - get distinct records on StationID, Latitude, Longitude, SampleDepth - also set row names to StationID
  scb.taxa <- benthic_data %>% dplyr::select(StationID, Latitude, Longitude, SampleDepth) %>%
    dplyr::distinct()

  # *Log scb.taxa

  # data prep step 3 - set up the scb predictors
  row.names(scb.predictors) <- scb.taxa$StationID
  scb.predictors <- as.matrix(scb.predictors)
  # *Log scb.predictors

  # data prep step 4 - Filter to replicate one and get distinct values on StationID Taxa and Abundance
  scb.taxa <- benthic_data %>%
    dplyr::filter(Replicate == 1) %>%
    dplyr::select(StationID, Taxa, Abundance) %>%
    dplyr::distinct()

  # Data prep step 5 - remove certain special characters from taxa name
  scb.taxa$Taxa <- gsub(" ", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub("(", "_", scb.taxa$Taxa, fixed = TRUE)
  scb.taxa$Taxa <- gsub(")", "_", scb.taxa$Taxa, fixed = TRUE)

  # Data prep step 6 - pivot the data out wide and make it a data.frame
  scb.taxa <- scb.taxa %>%
    tidyr::pivot_wider(id_cols = "StationID", names_from = "Taxa",
                       values_from = "Abundance", values_fn = list(Abundance = list))
  scb.taxa <- as.data.frame(scb.taxa)

  # Data prep step 7 - ...
  scb.taxa <- scb.taxa[, -1]

  # data prep step 8 - remove Abundance. from column names
  colnames(scb.taxa) <- gsub("Abundance.", "", colnames(scb.taxa))

  # data prep step 9 - Replace NAs with zero.
  scb.taxa[scb.taxa == "NULL"] <- 0
  scb.taxa = as.data.frame(lapply(scb.taxa, as.numeric))
  row.names(scb.taxa) <- row.names(scb.predictors)

  # RIVPACS calculations. By default the functions use the example user data.
  socal <- SoCalRivpacs(observed.predictors = scb.predictors, observed.taxa = scb.taxa, logfile = logfile, verbose = verbose)

  # the stations column of the oe table dataframe wwas being returned as a factor. Need to make that a character
  socal$oe.table <- socal$oe.table %>%
    mutate_if(is.factor,as.character)

  # **Please Log socal$oe.table**

  # Get the distinct values in benthic data based on StationID Replicate SampleDate Stratum
  benthic_data <- benthic_data %>%
    dplyr::select(StationID, Replicate, SampleDate, Stratum) %>%
    dplyr::distinct()

  # Calculate RIVPACS Scores

  # Riv step 0 - select stations and Observed/Expected
  riv0 <- socal$oe.table %>%
    dplyr::select(stations, O.over.E)

  # Riv step 1 - join with benthic data
  riv1 <- riv0 %>%
    dplyr::rename(StationID = stations, Score = O.over.E) %>%
    dplyr::full_join(benthic_data) %>%
    dplyr::mutate(Index = "RIVPACS")

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
                         logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ),
                         verbose = T
                         ) {

  # Pcutoff is the probability cutoff

  # Names of predictor variables.
  predictor.variables <- c("Latitude", "Longitude", "SampleDepth")

  # ----- Format Observed Data -----

  FormatObservedData <- function() {

    # Align observed (user) data columns with reference data columns. Columns in same
    # order. Observed data may have a different number of taxa (columns) than
    # reference data.

    # Convert observed.taxa to presence/absence (0/1)
    tmp.pa <- observed.taxa
    #tmp.pa <- lapply(tmp.pa, as.numeric)
    tmp.pa[tmp.pa > 0] <- 1

    # Align rows using predictor variables.
    tmp.pa <- tmp.pa[row.names(observed.predictors), ]    # !!! is this required???


    # Container matrix.
    n.observed.sites <- dim(tmp.pa)[1]
    n.reference.taxa <- dim(reference.taxa)[2]

    observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa),
                               nrow = n.observed.sites, ncol = n.reference.taxa,
                               dimnames = list(rownames(tmp.pa), names(reference.taxa)))

    # Fill container with observed data.
    col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])
    #tmp.pa <- as.data.frame(lapply(tmp.pa, as.numeric))

    for(i in 1:n.reference.taxa) {

      if(!is.na(col.match[i])) observed.taxa.pa[, i] <- tmp.pa[, col.match[i]]

    }

    # The matrix observed.taxa.pa contains the observed.scores used for O/E.
    return(observed.taxa.pa)

  }

  observed.data <- FormatObservedData()

  # ----- Calculate Expected Data -----

  CalculateExpectedData <- function() {

    # Calculate probability of sites belonging to groups. Follow RIVPACS assumption
    # of weighting the group probabilities by reference group size. Flags outlier
    # sites, using the chi-squared statistic.

    # Definitions.
    n.predictor.variables <- length(predictor.variables)
    group.size <- table(reference.groups)
    n.groups <- length(group.size)

    # Chi-squared values for flagging outlier samples.
    degrees.freedom <- min(c(n.predictor.variables, (n.groups - 1)))
    crit.01 <- qchisq(0.99, df = degrees.freedom)
    crit.05 <- qchisq(0.95, df = degrees.freedom)

    # Container for probabilities.
    n.observed.sites.filtered <- dim(observed.predictors)[[1]]

    group.probabilities <- matrix(rep(0, n.observed.sites.filtered * n.groups),
                                  nrow = n.observed.sites.filtered,
                                  dimnames = list(dimnames(observed.predictors)[[1]],
                                                  dimnames(group.means)[[1]]))

    # Container for outlier flags and minimum distance.
    outlier.flag <- data.frame(outlier.05 = rep(0, n.observed.sites.filtered),
                               outlier.01 = rep(0, n.observed.sites.filtered),
                               min.distance = rep(0, n.observed.sites.filtered),
                               row.names = dimnames(observed.predictors)[[1]])

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

    # Occurrence frequencies of all taxa in the reference groups.
    freq.in.group <- apply(reference.taxa, 2,
                           function(x){tapply(x, reference.groups, function(y){sum(y) / length(y)})})

    # Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).
    predicted.prob.all <- group.probabilities %*% freq.in.group

    # predicted.prob.all are the predicted (expected) probabilites.
    expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)

  }

  expected.data <- CalculateExpectedData()

  # ----- Calculate Scores -----

  CalculateScores <- function() {

    observed.score <- vector(mode = "numeric", length = expected.data$n)
    expected.score <- vector(mode = "numeric", length = expected.data$n)
    BC <- vector(mode = "numeric", length = expected.data$n) # Bray-Curtis dissimilarity



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

    stats <- data.frame(stations = row.names(observed.predictors),
                        O = observed.score,
                        E = round(expected.score, digits = 4),
                        O.over.E = round(O.over.E, digits = 4))

    #     stats <- data.frame(stations = row.names(observed.predictors),
    #                         O = observed.score,
    #                         E = expected.score,
    #                         O.over.E = O.over.E)

    stats$outlier.05 <- expected.data$outliers$outlier.05
    stats$outlier.01 <- expected.data$outliers$outlier.01

    # Convert to "PASS" or "FAIL"
    stats$outlier.05[stats$outlier.05 == 0] <- "PASS"
    stats$outlier.05[stats$outlier.05 == 1] <- "FAIL"

    stats$outlier.01[stats$outlier.01 == 0] <- "PASS"
    stats$outlier.01[stats$outlier.01 == 1] <- "FAIL"

    #   mean.O.over.E <- mean(OE.stats$O.over.E)
    #   stdev.O.over.E <- sqrt(var(OE.stats$O.over.E))


    return(stats)

  }

  results <- list(oe.table = CalculateScores(),
                  observed = observed.data,
                  predicted = expected.data$predicted,
                  Pcutoff = Pcutoff,
                  region = "scb")

}



