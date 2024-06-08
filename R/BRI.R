

BRI <- function(BenthicData, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T)
{

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\nBEGIN: BRI function.\n', logfile = logfile, verbose = verbose)

  writelog('*** DATA *** Input to BRI function - BRI-step0.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(BenthicData, logfile = file.path(dirname(logfile), 'BRI-step0.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  writelog('*** DATA *** sqo.list.new - which gets joined with BenthicData', logfile = logfile, verbose = verbose)
  writelog(sqo.list.new, logfile = file.path(dirname(logfile), 'BRI-sqo.list.new.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  out <- BenthicData %>%
    left_join(sqo.list.new, by = c('Taxon' = 'TaxonName'))

  writelog('*** DATA *** Benthic data joined with sqo.list.new', logfile = logfile, verbose = verbose)
  writelog(out, logfile = file.path(dirname(logfile), 'BRI-step1.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  out <- out %>%
    # I assume that the next line is something they had in there as a method of removing duplicates
    # for this reason, this next line will likely be eliminated.
    # They grouped by all the columns that were selected (In query BRI - 1)
    # Instead, if need be we can use something from dplyr that deals with duplicates
    # I actually found that it didn't appear to make a difference
    filter(!is.na(ToleranceScore)) %>%
    #rename(Stratum) %>%
    select(Stratum, StationID, SampleDate, Replicate, Taxon, Abundance, ToleranceScore)

  writelog('*** DATA *** Remove NA tolerance scores and select only the columns: Stratum, StationID, SampleDate, Replicate, Taxon, Abundance, ToleranceScore', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(out, logfile = file.path(dirname(logfile), 'BRI-step2.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  out <- out %>%
    # End of BRI - 2 query. Begin BRI - 3 query
    mutate(
      fourthroot_abun = Abundance ** 0.25,
      tolerance_score = fourthroot_abun * ToleranceScore
    )

  writelog('*** DATA *** Calculate fourthroot of abundance - also multiply by tolerancs score - BRI-step3.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(out, logfile = file.path(dirname(logfile), 'BRI-step3.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)

  writelog('Next get the Score - group by Stratum, StationID, SampleDate, Replicate and do: (sum of the tolerance scores)/(sum of fourthroot abundances)', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('Then Get Categories (CASQO Technical Manual 3rd Edition Page 72 - Table 4.24)', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- < 39.96 is Reference', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- >=39.96 and <49.15 is Low', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- >=49.15 and <73.27 is Reference', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('---- >=73.27 is High', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  out <- out %>%
    # End of BRI - 3. Begin BRI - 4
    group_by(
      Stratum, StationID, SampleDate, Replicate
    ) %>%
    summarize(
      Score = sum(tolerance_score, na.rm = T) / sum(fourthroot_abun, na.rm = T)
    ) %>%
    # Output the BRI category given the BRI score and the thresholds for Southern California Marine Bays
    # CASQO Technical Manual 3rd Edition Page 72 - Table 4.24
    mutate(
      Category = case_when( (Score < 39.96) ~ "Reference",
                            (Score >= 39.96 & Score < 49.15) ~ "Low Disturbance",
                            (Score >= 49.15 & Score < 73.27) ~ "Moderate Disturbance",
                            (Score >= 73.27) ~ "High Disturbance"
      )) %>%
    # Output the BRI category score given the category for thresholds for Southern CA Marine Bays
    mutate(
      `Category Score` = case_when( (Category == "Reference") ~ 1,
                                    (Category == "Low Disturbance") ~ 2,
                                    (Category == "Moderate Disturbance") ~ 3,
                                    (Category == "High Disturbance") ~ 4 )
    ) %>%
    dplyr::mutate(Index = "BRI")


  writelog('*** DATA *** Final BRI dataframe: BRI-final.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(out, logfile = file.path(dirname(logfile), 'BRI-final.csv'), filetype = 'csv', verbose = verbose, prefix = hyphen.log.prefix)


  writelog('\nEND: BRI function.\n', logfile = logfile, verbose = verbose)

  return(out)
}
