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

  "tox_categories"

  toxresults$stationid <- replace(toxresults$stationid, toxresults$stationid %>% as.character() == '0', '0000')
  toxresults$result <- replace(toxresults$result, toxresults$result == -88, NA_real_)


  if (!('lab' %in% tolower(names(toxresults)))) {
    toxresults <- toxresults %>%
      mutate(
        lab = 'Not Recorded'
      )
  }

  if (!('fieldrep' %in% tolower(names(toxresults)) )) {
    toxresults <- toxresults %>%
      mutate(
        fieldrep = 1
      )
  }

  # here we separate the controls from the rest of the samples
  controls <- toxresults %>%
    filter(
      stationid == '0000',
      !(sampletypecode %in% results.sampletypes)
    )

  # get rid of the controls. They will be merged later
  results <- toxresults %>%
    filter(
      ( (stationid != '0000') | (sampletypecode %in% results.sampletypes) )
    )

  summary <- results %>%
    inner_join(
      controls,
      by = c('toxbatch','species','labrep','lab'),
      suffix = c('','_control')
    ) %>%
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
tox.sqo <- function(toxresults, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ), verbose = T) {

  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  # It is common that there is a station called "0" which is actually supposed to be a string of four zeros ('0000')
  toxresults$stationid <- replace(toxresults$stationid, toxresults$stationid %>% as.character() == '0', '0000')

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

  if ('stratum' %in% names(toxresults)) {
    tox_nonintegrated <- tox_nonintegrated %>%
      left_join(
        toxresults %>% select(stationid, stratum) %>% distinct() %>% filter(!is.na(stratum)),
        by = 'stationid'
      )
  }



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
  return(tox_final)
}
