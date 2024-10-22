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
SQOUnified <- function(benthic = NULL, chem = NULL, tox = NULL, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.Rmd' ), verbose = F) {

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)


  #load("data/site_assessment_criteria.RData")

  if (all(is.null(c(benthic,chem,tox)))){
    stop("
      No data was provided.
      Please provide benthic, chemistry and toxicity data to get the integrated site assessments
    ")
  }
  # check the data coming in before anything
  checkdata(benthic, chem, tox, logfile = logfile, verbose = verbose)



  # Compute ALL SQO scores
  writelog('Compute ALL SQO scores', logfile = logfile, verbose = verbose)

  # ---- Benthic ----
  if (!is.null(benthic)) {

    benthiclogfile <- file.path( dirname(logfile), 'Benthic', 'benthiclog.Rmd' )
    benthiclibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(benthiclogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = benthiclibs)

    benthic <- benthic.sqo(benthic, logfile = benthiclogfile, verbose = verbose) %>%
      mutate(LOE = 'Benthic') %>%
      select(StationID, Replicate, SampleDate, LOE, Index, Score, Category, `Category Score`) %>%
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
      select(-c(Replicate,SampleDate))

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

  # ---- Chemistry ----
  if (!is.null(chem)) {

    chemlogfile <- file.path( dirname(logfile), 'Chemistry', 'chemlog.Rmd' )
    chemlibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(chemlogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = chemlibs)

    chem <- chem.sqo(chem, logfile = chemlogfile, verbose = verbose) %>%
      mutate(LOE = 'Chemistry') %>%
      select(StationID, LOE, Index, Score, Category, `Category Score`)
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

  # ---- Toxicity ----
  if (!is.null(tox)) {

    toxlogfile <- file.path( dirname(logfile), 'Toxicity', 'toxlog.Rmd' )
    toxlibs <- c('tidyverse', 'DT', 'knitr', 'rmarkdown', 'SQOUnified')

    init.log(toxlogfile, base.func.name = sys.call(), type = 'RMarkdown', current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, libraries = toxlibs)

    tox <- tox.sqo(tox, logfile = toxlogfile, verbose = verbose) %>%
      mutate(LOE = 'Toxicity') %>%
      select(StationID, LOE, Index, Score, Category, `Category Score`)
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

  integrated <- bind_rows(benthic, chem, tox) %>%
    filter(
      grepl("SQO",Index)
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
    left_join(
      site_assessment_criteria,
      by = c("Benthic","Chemistry","Toxicity")
    ) %>%
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
    benthic, chem, tox, integrated
  ) %>%
  arrange(
    StationID, LOE, Index
  )

  writelog('Done computing ALL SQO scores', logfile = logfile, verbose = verbose)

  return(out)


}
