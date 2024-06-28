#' Compute the SQO index scores.
#' Based on the Technical Manual - June 2021 edition
#' https://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf
#'
#' @param benthic_data A data file string name that we want to compute SQO scores for.
#' @param SQO A list of the type of SQO scores that we want to compute
#'     (e.g., \code{c("MAMBI", "RBI")}).
#'     The default is \code{"all"}, meaning that all scores will be computed.
#' @usage
#' data(benthic_data)
#' data(site_assessment_criteria)
#' @examples
#' SQOUnified(benthic_data = benthic_data, LOE = "all")
#' SQOUnified(benthic_data = benthic_data, LOW = 'tox')
#'
#' @importFrom dplyr case_when full_join select rename mutate arrange
#' @importFrom purrr map
#' @importFrom tidyr spread


#' @export
SQOUnified <- function(benthic = NULL, chem = NULL, tox = NULL, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.Rmd' ), verbose = F) {

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

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
