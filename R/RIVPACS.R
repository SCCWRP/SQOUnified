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
#' @param retrofit_taxonomy Logical. If TRUE (default), modern SCAMIT Edition 14 taxonomy in the
#'    submitted data is retrofitted back to SQO-compatible names via \code{\link{benthicdata_prep}}
#'    before the index is calculated. Set FALSE if the data has already been prepped/retrofitted.
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
RIVPACS <- function(BenthicData,
                                retrofit_taxonomy = TRUE,
                                logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'RIVPACS_generic_log.Rmd'),
                                verbose = FALSE,
                                knitlog = FALSE)
{
  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)

  writelog('\n## BEGIN: Generic RIVPACS function.\n', logfile = logfile, verbose = verbose)

  # Reference data (socal.reference.*, socal.example.*) are available as package datasets — see ?socal.reference.taxa
  # SoCalRivpacs.2() is available in the package namespace from R/SoCalRivpacs2.R

  # Standardize and (optionally) retrofit the submitted taxonomy to SQO-compatible names
  BenthicData <- benthicdata_prep(BenthicData, retrofit = retrofit_taxonomy, logfile = logfile, verbose = verbose)$benthic_data

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
  # Collapse any duplicate NoOrganismsPresent rows within a sample so we don't
  # emit multiple defaunated rows per (stationid, sampledate, replicate).
  defaunated <- BenthicData %>%
    filter(taxon == "NoOrganismsPresent") %>%
    distinct(stationid, sampledate, replicate, .keep_all = TRUE) %>%
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
    '### RIVPACS Step 1 - Predictors and taxa matrices prepared\n',
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
    distinct(stationid, sampledate, replicate, .keep_all = TRUE)

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
    '### RIVPACS Step 2 - O/E details\n',
    logfile = logfile,
    data = rivpacs.scores %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rivpacs.scores, logfile = logfile, filename = 'RIVPACS_generic-step1-OE_details.csv', linktext = 'Download RIVPACS O/E details', verbose = verbose)

  # Drop any samples from defaunated that already have a calculated score, so
  # the bind_rows + join below can never emit two rows per sample.
  defaunated <- defaunated %>%
    anti_join(rivpacs.scores, by = c("stationid", "sampledate", "replicate"))

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
    '### RIVPACS Final - RIVPACS Scores\n',
    logfile = logfile,
    data = rivpacs.scores.2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rivpacs.scores.2, logfile = logfile, filename = 'RIVPACS_generic-final_scores.csv', linktext = 'Download RIVPACS scores', verbose = verbose)

  writelog('\n## END: Generic RIVPACS function.\n', logfile = logfile, verbose = verbose)

  return(rivpacs.scores.2)
}



# This function is private, used by alt.RIVPACS.generic
# RIVPACS discriminant function model for Southern California
# Based on SoCalRivpacs v2.R from the SccwrpRivpacs package
# Reference data objects (socal.reference.*) are loaded from R/sysdata.rda
SoCalRivpacs.2 <- function(Pcutoff = 0.5,
                           reference.groups = socal.reference.groups,
                           observed.predictors = socal.example.habitat,
                           reference.taxa = socal.reference.taxa,
                           group.means = socal.reference.group.means,
                           reference.cov = socal.reference.covariance,
                           observed.taxa = socal.example.taxa) {

  # Pcutoff is the probability cutoff

  # Names of predictor variables.
  predictor.variables <- c("Latitude", "Longitude", "SampleDepth")

  # ----- Format Observed Data -----

  FormatObservedData <- function() {

    # Align observed (user) data columns with reference data columns.
    # Convert observed.taxa to presence/absence (0/1)
    tmp.pa <- observed.taxa
    tmp.pa[tmp.pa > 0] <- 1

    # Align rows using predictor variables.
    tmp.pa <- tmp.pa[row.names(observed.predictors), ]

    # Container matrix.
    n.observed.sites <- dim(tmp.pa)[[1]]
    n.reference.taxa <- dim(reference.taxa)[[2]]

    observed.taxa.pa <- matrix(rep(0, times = n.observed.sites * n.reference.taxa),
                               nrow = n.observed.sites,
                               dimnames = list(rownames(tmp.pa), names(reference.taxa)))

    # Fill container with observed data.
    col.match <- match(dimnames(observed.taxa.pa)[[2]], dimnames(tmp.pa)[[2]])

    for (i in 1:n.reference.taxa) {
      if (!is.na(col.match[i])) observed.taxa.pa[, i] <- tmp.pa[, col.match[i]]
    }

    return(observed.taxa.pa)
  }

  observed.data <- FormatObservedData()

  # ----- Calculate Expected Data -----

  CalculateExpectedData <- function() {

    # Calculate probability of sites belonging to groups.
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
    for (i in 1:n.observed.sites.filtered) {

      # Squared Mahalanobis distance from current sample to each group mean.
      distances <- mahalanobis(group.means,
                               observed.predictors[i, ],
                               reference.cov,
                               inverted = TRUE)

      group.probabilities[i, ] <- group.size * exp(-0.5 * distances)
      group.probabilities[i, ] <- group.probabilities[i, ] / sum(group.probabilities[i, ])

      # Outlier criteria is minimum distance.
      outlier.flag$min.distance[i] <- min(distances)

      # Check for outliers.
      if (outlier.flag$min.distance[i] > crit.05) outlier.flag[i, "outlier.05"] <- 1
      if (outlier.flag$min.distance[i] > crit.01) outlier.flag[i, "outlier.01"] <- 1
    }

    # Occurrence frequencies of all taxa in the reference groups.
    freq.in.group <- apply(reference.taxa, 2,
                           function(x) { tapply(x, reference.groups, function(y) { sum(y) / length(y) }) })

    # Matrix algebra form of the RIVPACS combining formula (Clarke et al. 2003, Eq. 4).
    predicted.prob.all <- group.probabilities %*% freq.in.group

    expected.data <- list(predicted = predicted.prob.all, outliers = outlier.flag, n = n.observed.sites.filtered)
  }

  expected.data <- CalculateExpectedData()

  # ----- Calculate Scores -----

  CalculateScores <- function() {

    observed.score <- vector(mode = "numeric", length = expected.data$n)
    expected.score <- vector(mode = "numeric", length = expected.data$n)
    BC <- vector(mode = "numeric", length = expected.data$n)

    for (i in 1:expected.data$n) {

      predicted.prob <- expected.data$predicted[i, ]
      taxa.subset <- names(predicted.prob)[predicted.prob >= Pcutoff]
      expected.prob <- predicted.prob[taxa.subset]
      observed.pa <- observed.data[i, taxa.subset]

      observed.score[i] <- sum(observed.pa)
      expected.score[i] <- sum(expected.prob)
      BC[i] <- sum(abs(observed.pa - expected.prob)) /
        (observed.score[i] + expected.score[i])
    }

    O.over.E <- observed.score / expected.score

    stats <- data.frame(sample_id = row.names(observed.predictors),
                        O = observed.score,
                        E = round(expected.score, digits = 4),
                        O.over.E = round(O.over.E, digits = 4))

    stats$outlier.05 <- expected.data$outliers$outlier.05
    stats$outlier.01 <- expected.data$outliers$outlier.01

    # Convert to "PASS" or "FAIL"
    stats$outlier.05[stats$outlier.05 == 0] <- "PASS"
    stats$outlier.05[stats$outlier.05 == 1] <- "FAIL"

    stats$outlier.01[stats$outlier.01 == 0] <- "PASS"
    stats$outlier.01[stats$outlier.01 == 1] <- "FAIL"

    return(stats)
  }

  results <- list(oe.table = CalculateScores(),
                  observed = observed.data,
                  predicted = expected.data$predicted,
                  Pcutoff = Pcutoff,
                  region = "scb")
}
