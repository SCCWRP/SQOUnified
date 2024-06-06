#---- tox summary ----
#' Get the Statistical Summary for Tox Results Data
#'
#' @description
#' This function will calculate the tox summary,
#'
#' given an argument which is a dataframe of tox results data
#'
#' @param toxresults a dataframe with the following columns: stationid, toxbatch, species, sampletypecode
#'    matrix, labrep, result. This data must also include the control samples
#'    (stationcode 0000, sampletypecode CNEG etc)
#'
#'    The input dataframe is structured as follows
#'
#'    \strong{\code{toxresults}} -  a dataframe that contains the toxicity results
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the location;
#'
#'    \strong{\code{toxbatch}} - the toxbatch id - used to join with the control sample
#'
#'    \strong{\code{species}} - The Genus and species of the animale that was tested
#'
#'    \strong{\code{sampletypecode}} - The sampletype used Grab, CNEG etc. Control samples must be included
#'
#'    \strong{\code{matrix}} - Whole Sediment, Sediment Water Interface, etc. Probably useless to include.
#'                    I Just have it to make sure they dont put Reference Toxicant
#'
#'    \strong{\code{labrep}} - There should be 5 per station, species pair
#'
#'    \strong{\code{result}} - the percentage that survived the test, or had normal development
#'
#'
#' @usage tox.summary(toxresults)
#'
#' @examples
#' data(tox_sampledata)
#' tox.summary(tox_sampledata)
#'
#' @importFrom plyr rbind.fill
#' @importFrom stats t.test
#' @importFrom tidyr separate
#' @import dplyr
#' @export

# Version 0.3.0 update - allow a user to select sampletypes to include - allows QA to be included if a user so chooses
# DEFAULT leave it out and do only grabs
tox.summary <- function(toxresults, results.sampletypes = c('Grab'), logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T) {

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\nBEGIN: Tox Summary function.\n', logfile = logfile, verbose = verbose)
  writelog('*** DATA *** Tox results upon entry into tox.summary function - tox.summary-initial-input-step0.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(toxresults, logfile = file.path(dirname(logfile), 'tox.summary-initial-input-step0.csv'), filetype = 'csv', verbose = verbose)

  writelog("These are the sampletypes which are to be analyzed as actual results (Typically Grab or Grab/QA):", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste0(results.sampletypes, collapse = ', '), logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  "tox_categories"

  toxresults$stationid <- replace(toxresults$stationid, toxresults$stationid %>% as.character() == '0', '0000')
  toxresults$result <- replace(toxresults$result, toxresults$result < 0, NA_real_)

  writelog('*** DATA *** Tox results after replacing 0 stationid with 0000 and negative results with NA_real_ - tox.summary-input-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(toxresults, logfile = file.path(dirname(logfile), 'tox.summary-input-step1.csv'), filetype = 'csv', verbose = verbose)


  if (!('lab' %in% tolower(names(toxresults)))) {
    toxresults <- toxresults %>%
      mutate(
        lab = 'Not Recorded'
      )
    writelog('Added lab column to input data - filled with Not Recorded', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  }

  if (!('fieldrep' %in% tolower(names(toxresults)) )) {
    toxresults <- toxresults %>%
      mutate(
        fieldrep = 1
      )
    writelog('Added fieldrep column to input data - filled with 1', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  }

  # here we separate the controls from the rest of the samples
  controls <- toxresults %>%
    filter(
      stationid == '0000',
      !(sampletypecode %in% results.sampletypes)
    )

  writelog('*** DATA *** Controls - tox.summary-controls.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(controls, logfile = file.path(dirname(logfile), 'tox.summary-controls.csv'), filetype = 'csv', verbose = verbose)

  # get rid of the controls. They will be merged later
  results <- toxresults %>%
    filter(
      ( (stationid != '0000') | (sampletypecode %in% results.sampletypes) )
    )

  writelog('*** DATA *** Results separated from controls - tox.summary-results-without-controls.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(results, logfile = file.path(dirname(logfile), 'tox.summary-results-without-controls.csv'), filetype = 'csv', verbose = verbose)

  writelog("Join results and controls on 'toxbatch', 'species', 'labrep', and 'lab' tox.summary", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  summary <- results %>%
    inner_join(
      controls,
      by = c('toxbatch','species','labrep','lab'),
      suffix = c('','_control')
    )

  writelog('*** DATA *** Results separated from controls - tox.summary-results-joined-with-controls-step2.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(summary, logfile = file.path(dirname(logfile), 'tox.summary-results-joined-with-controls-step2.csv'), filetype = 'csv', verbose = verbose)

  writelog('Group by lab, stationid, toxbatch, species, fieldrep, sampletypecode', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('calculate two sided t-test p value (t.test function in R) - the two samples being the result and the control result', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('t.test function in R by default will treat them as two samples from two populations (not assuming equal variances)', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('The function call to t.test is this:', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- t.test(x = result, y = result_control, mu = 0, var.equal = F, alternative = 'two.sided')$p.value / 2", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("Further explanation of the function call:", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- x = result and y = result_control: These are the two vectors (samples) being compared.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- mu = 0: The null hypothesis is that the difference between the means of the two samples is zero.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- var.equal = F: The variances of the two samples are not assumed to be equal (Welch's t-test).", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- alternative = 'two.sided': The alternative hypothesis is that the means are different (two-sided test).", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- $p.value / 2: The p-value of the t-test is divided by 2. This division suggests that the intention is to obtain a one-tailed p-value from a two-tailed test.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)


  summary <- summary %>%
    # filter(
    #   sampletypecode != 'QA'
    # ) %>%
    group_by(
      lab, stationid, toxbatch, species, fieldrep, sampletypecode
    ) %>%
    summarize(
      p = tryCatch({
        # it errors out when you pass two constant vectors as the first two arguments
        t.test(x = result, y = result_control, mu = 0, var.equal = F, alternative = 'two.sided')$p.value / 2
      },
      warning = function(w){
        # If there was a warning, we still want the output of the function
        return(
          t.test(x = result, y = result_control, mu = 0, var.equal = F, alternative = 'two.sided')$p.value / 2
        )
      },
      error = function(err){
        if(all(result == result_control)){
          # if the two vectors were exactly the same, they are not significantly different for obvious reasons
          return(1)
        } else {
          # if they were two completely different constant vectors, then they are significantly different
          return(0)
        }
      }
      ),
      pct_result = mean(result, na.rm = T),
      pct_control = mean(result_control, na.rm = T),
      pct_result_adj = (pct_result / pct_control) * 100,
      stddev = sd(result, na.rm = T),
      cv = stddev / pct_result,
      n = n()
    ) %>%
    ungroup() %>%
    mutate(
      sigdiff = if_else(p < 0.05, TRUE, FALSE)
    ) %>%
    left_join(
      tox_categories, by = 'species'
    ) %>%
    mutate(
      # CASQO Technical Manual page 45 - 47
      sqo_category_value_initial = case_when(
        # if the endpoint method is not Growth, we look at the non control adjusted mean percentage to determine nontoxicity
        (endpoint_method != 'Growth') & (pct_result >= nontox) ~ 1,
        # for Growth, we consider the control adjusted mean percentage
        (endpoint_method == 'Growth') & (pct_result_adj >= nontox) ~ 1,
        # For all the other toxicity SQO categories, we always look at the control adjusted mean percentage
        # if lowtox <= pct_result_adj < nontox, put it in the low toxicity category - always
        pct_result_adj >= lowtox ~ 2,
        # if modtox <= pct_result_adj < lowtox, put it in the moderate toxicity category - always
        pct_result_adj >= modtox ~ 3,
        # below lower bound of moderate toxicity renders it in the category of high toxicity - always
        pct_result_adj < modtox ~ 4,
        TRUE ~ NA_real_
      ),
      sqo_category_value = if_else(
        # I noticed by reading the manual that the low and moderate scores got improved by one category if the means were not significantly different
        sqo_category_value_initial %in% c(2,3) & !sigdiff,
        sqo_category_value_initial - 1,
        sqo_category_value_initial
      ),
      sqo_category = case_when(
        sqo_category_value == 1 ~ "Nontoxic",
        sqo_category_value == 2 ~ "Low Toxicity",
        sqo_category_value == 3 ~ "Moderate Toxicity",
        sqo_category_value == 4 ~ "High Toxicity",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-c(nontox, lowtox, modtox, hightox)) %>%
    rename(
      lab = lab,
      stationid = stationid,
      species = species,
      toxbatch = toxbatch,
      sampletypecode = sampletypecode,
      `P Value` = p,
      `Mean` = pct_result,
      `Control Adjusted Mean` = pct_result_adj,
      `Endpoint Method` = endpoint_method,
      `Standard Deviation` = stddev,
      `Coefficient of Variance` = cv,
      Score = sqo_category_value,
      Category = sqo_category
    )
  writelog('*** DATA *** Tox summary dataframe final - tox.summary-final.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(summary, logfile = file.path(dirname(logfile), 'tox.summary-final.csv'), filetype = 'csv', verbose = verbose)

  writelog('\nEND: Tox Summary function.\n', logfile = logfile, verbose = verbose)

  return(summary)
}


# ---- Tox SQO ----
#' Get Tox SQO Scores and Categories
#'
#' @description
#' This function will calculate the tox SQO scores for stations given a dataframe structured as described
#' in the details or the Arguments section. The funtion will get SQO scores for each individual test
#' conducted for a station, as well as the integrated Toxicity LOE SQO score and category
#'
#' @param toxresults a dataframe with the following columns: stationid, toxbatch, species, sampletypecode
#'    matrix, labrep, result. This data must also include the control samples
#'    (stationcode 0000, sampletypecode CNEG etc)
#'
#'    The input dataframe is structured as follows:
#'
#'    \strong{\code{stationid}} - an alpha-numeric identifier of the location;
#'
#'    \strong{\code{toxbatch}} - the toxbatch id - used to join with the control sample
#'
#'    \strong{\code{species}} - The Genus and species of the animale that was tested
#'
#'    \strong{\code{sampletypecode}} - The sampletype used Grab, CNEG etc. Control samples must be included
#'
#'    \strong{\code{matrix}} - Whole Sediment, Sediment Water Interface, etc. Probably useless to include.
#'                    I Just have it to make sure they dont put Reference Toxicant
#'
#'    \strong{\code{labrep}} - There should be 5 per station, species pair
#'
#'    \strong{\code{result}} - the percentage that survived the test, or had normal development
#'
#' @usage tox.sqo(toxresults)
#'
#' @examples
#' data(tox_sampledata)
#' tox.sqo(tox_sampledata)
#'
#' @export
tox.sqo <- function(toxresults, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T) {

  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\nBEGIN: Tox SQO function.\n', logfile = logfile, verbose = verbose)

  writelog('*** DATA *** Tox results upon entry into tox.sqo function - tox.sqo-initial-input-step0.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(toxresults, logfile = file.path(dirname(logfile), 'tox.sqo-initial-input-step0.csv'), filetype = 'csv', verbose = verbose)

  # It is common that there is a station called "0" which is actually supposed to be a string of four zeros ('0000')
  toxresults$stationid <- replace(toxresults$stationid, toxresults$stationid %>% as.character() == '0', '0000')

  writelog('*** DATA *** Tox results after replacing 0 with 0000 in stationid column - tox.sqo-initial-input-step0.1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(toxresults, logfile = file.path(dirname(logfile), 'tox.sqo-initial-input-step0.1.csv'), filetype = 'csv', verbose = verbose)

  writelog('Prep the tox data by taking the summary - then categorizing the endpoint method to acute or sublethal', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  tox_nonintegrated <- tox.summary(toxresults, logfile = logfile, verbose = verbose) %>%
    select(
      stationid,
      species,
      `Endpoint Method`,
      Score,
      Category
    ) %>%
    separate(
      species,
      c("Genus","Species")
    ) %>%
    select(-Species) %>%
    mutate(
      Index = paste(Genus, `Endpoint Method`),
      `Test Type` = case_when(
        `Endpoint Method` %in% c('Growth','NormDev') ~ 'Sublethal',
        `Endpoint Method` == 'Survival' ~ 'Acute',
        TRUE ~ NA_character_
      )
    ) %>%
    select(stationid, Index, `Test Type`, Score, Category) %>%
    mutate(
      `Category Score` = Score # just for purposes of the very final unified output, all three LOE's in one table
    )

  writelog('*** DATA *** After taking summary and categorizing endpoint method - tox.sqo-step1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(tox_nonintegrated, logfile = file.path(dirname(logfile), 'tox.sqo-step1.csv'), filetype = 'csv', verbose = verbose)

  if ('stratum' %in% names(toxresults)) {

    writelog('Include the stratum if it was provided in the initial toxresults', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    tox_nonintegrated <- tox_nonintegrated %>%
      left_join(
        toxresults %>% select(stationid, stratum) %>% distinct() %>% filter(!is.na(stratum)),
        by = 'stationid'
      )
    writelog('*** DATA *** After tacking on the stratum - tox.sqo-step1.1.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(tox_nonintegrated, logfile = file.path(dirname(logfile), 'tox.sqo-step1.1.csv'), filetype = 'csv', verbose = verbose)
  }

  writelog("# For Toxicity, we take the mean of SQO scores among the tox tests at the site", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# CASQO manual (3rd edition) page 109", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog("# The score is the mean of the two, however, if one is missing, the score should come out to NA", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# This means that in this particular case, it should be na.rm = F", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# (NOTE: Not for offshore sites - offhsore sites only require the Amphipod test)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog("# This is not explicitly mentioned on page 109 of the manual, where the instructions are given for integrating the results of the tox tests to determine the Tox LOE", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# However, we can understand that it must be na.rm = F because there are essentially two types of test methods, Acute and Sublethal (Page 85) and in the intro", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# to the instructions for tox assessment on page 84, first paragraph, it says the following:", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog("#   Multiple toxicity tests are needed to assess toxicity because no single method exists that can capture the", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("#   full spectrum of potential contaminant effects. Toxicity assessment under the California", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("#   Sediment Quality Objectives (CASQO) framework requires information from two types of tests:", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("#     1) short-term amphipod survival and 2) a sublethal test.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog("# For this reason, if one of the test results is missing, we will declare the Tox LOE to be indeterminate or unknown", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# So here, we will actually need to check if both tests were present", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("# Acute tests mean the endpoint method is 'Survival' and the Sublethal tests mean the endpoint is 'Growth' or 'NormDev' (Normal Development)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)



  tox_integrated <- tox_nonintegrated %>%
    group_by(stationid) %>%
    summarize(
      # For Toxicity, we take the mean
      # CASQO manual (3rd edition) page 109

      # PER DARRIN GREENSTEIN, in July 2022 - yes, the score is the mean of the two, however, if one is missing, the score should come out to NA
      # This means that in this particular case, it should be na.rm = F

      # This is not explicitly mentioned on page 109 of the manual, where the instructions are given for integrating the results of the tox tests to determine the Tox LOE
      # However, we can understand that it must be na.rm = F because there are essentially two types of test methods, Acute and Sublethal (Page 85) and in the intro
      #  to the instructions for tox assessment on page 84, first paragraph, it says the following:

      #   Multiple toxicity tests are needed to assess toxicity because no single method exists that can capture the
      #   full spectrum of potential contaminant effects. Toxicity assessment under the California
      #   Sediment Quality Objectives (CASQO) framework requires information from two types of tests:
      #     1) short-term amphipod survival and 2) a sublethal test.

      # For this reason, if one of the test results is missing, we will declare the Tox LOE to be indeterminate or unknown
      # So here, we will actually need to check if both tests were present
      # Acute tests mean the endpoint method is "Survival" and the Sublethal tests mean the endpoint is "Growth" or "NormDev" (Normal Development)

      #Score = ceiling(mean(Score, na.rm = F))
      has.all.tests = all(c('Acute','Sublethal') %in% `Test Type`),
      offshore = 'stratum' %in% names(.) && stratum %in% c('Outer Shelf','Mid Shelf','Inner Shelf','Channel Islands'),
      Score = case_when(
        # offhsore sites only require the Amphipod test - per Ken Schiff August 1, 2022
        offshore ~ ceiling(mean(Score, na.rm = T)),
        has.all.tests ~ ceiling(mean(Score, na.rm = F)),
        TRUE ~ NA_real_
      )
    ) %>%
    mutate(
      Category = case_when(
        Score == 1 ~ "Nontoxic",
        Score == 2 ~ "Low Toxicity",
        Score == 3 ~ "Moderate Toxicity",
        Score == 4 ~ "High Toxicity",
        TRUE ~ NA_character_
      ),
      `Category Score` = Score,
      Index = "Integrated SQO"
    ) %>%
    mutate(
      `Category Score` = Score # just for purposes of the very final unified output, all three LOE's in one table
    ) %>%
    select(stationid, Index, Score, Category, `Category Score`)

  tox_final <- rbind.fill(tox_integrated, tox_nonintegrated) %>%
    rename(StationID = stationid) %>%
    arrange(
      StationID, Index, Category
    )

  writelog('*** DATA *** Final Tox dataframe: tox.sqo-final.csv', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(tox_final, logfile = file.path(dirname(logfile), 'tox.sqo-final.csv'), filetype = 'csv', verbose = verbose)

  writelog('\nEND: Tox SQO function.\n', logfile = logfile, verbose = verbose)


  return(tox_final)
}
