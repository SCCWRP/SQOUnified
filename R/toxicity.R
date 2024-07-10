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
tox.summary <- function(tox.summary.input, results.sampletypes = c('Grab'), logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = F) {

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\n## BEGIN: Tox Summary function.\n', logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'tox.summary.input.RData'
  if (verbose) {
    save(tox.summary.input, file = file.path( dirname(logfile), rawinput.filename ))
  }

  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n### Raw input to tox.summary:\n  ",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') # This will load a dataframe called 'tox.summary.input' into your environment"),
    data = tox.summary.input %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox.summary.input, logfile = logfile, filename = 'tox.summary-input-step0.csv', linktext = 'Download Raw Input to Tox Summary', verbose = verbose)

  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n### Ensure the results.sampletypes argument is defined for the RMarkdown:\n  ",
    logfile = logfile,
    code = paste0("results.sampletypes <- c('", paste0( results.sampletypes, collapse = "','" ), "')"  ),
    verbose = verbose
  )


  # write to the logs the sampletypes used for results
  writelog("These are the sampletypes which are to be analyzed as actual results (Typically Grab or Grab/QA):\n  ", logfile = logfile, verbose = verbose)
  writelog(paste0(results.sampletypes, collapse = ', '), logfile = logfile, verbose = verbose)

  # I think this needs to be here to load to the local environment
  "tox_categories"

  # replace stationid of 0 with the Bight QA StationID of 0000
  tox.summary.input$stationid <- replace(tox.summary.input$stationid, tox.summary.input$stationid %>% as.character() == '0', '0000')
  writelog(
    "\n### replace stationid of 0 with the Bight QA StationID of 0000:\n  ",
    logfile = logfile,
    code = "tox.summary.input$stationid <- replace(tox.summary.input$stationid, tox.summary.input$stationid %>% as.character() == '0', '0000')",
    data = tox.summary.input %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox.summary.input, logfile = logfile, filename = 'tox.summary-input-step0.1.csv', linktext = 'Download Input to Tox Summary with stationid of 0 replaced with 0000', verbose = verbose)

  # replace negative result values with NA_real_ (-88's and -99's shouldnt be included, nor should the result value ever be negative)
  tox.summary.input$result <- replace(tox.summary.input$result, tox.summary.input$result < 0, NA_real_)
  writelog(
    "\n### Replace negative result values with NA's:\n  ",
    logfile = logfile,
    code = "tox.summary.input$result <- replace(tox.summary.input$result, tox.summary.input$result < 0, NA_real_)",
    data = tox.summary.input %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox.summary.input, logfile = logfile, filename = 'tox.summary-input-step0.2.csv', linktext = 'Download Input to Tox Summary with negative result values replaced with NA', verbose = verbose)


  # Add lab column if it is not there already
  if (!('lab' %in% tolower(names(tox.summary.input)))) {

    tox.summary.input <- tox.summary.input %>%
      mutate(
        lab = 'Not Recorded'
      )

    writelog(
      "\n### Added lab column to input data - filled with Not Recorded:\n  ",
      logfile = logfile,
      code = "tox.summary.input <- tox.summary.input %>% mutate(lab = 'Not Recorded')",
      data = tox.summary.input %>% head(25),
      verbose = verbose
    )
    create_download_link(data = tox.summary.input, logfile = logfile, filename = 'tox.summary-input-step0.2.1.csv', linktext = 'Download Input to Tox Summary with lab column added', verbose = verbose)

  }

  # Add fieldrep column if it is not there already
  if (!('fieldrep' %in% tolower(names(tox.summary.input)) )) {
    tox.summary.input <- tox.summary.input %>%
      mutate(
        fieldrep = 1
      )
    writelog(
      "\n### Added fieldrep column to input data - filled with 1:\n  ",
      logfile = logfile,
      code = "tox.summary.input <- tox.summary.input %>% mutate(fieldrep = 1)",
      data = tox.summary.input %>% head(25),
      verbose = verbose
    )
    create_download_link(data = tox.summary.input, logfile = logfile, filename = 'tox.summary-input-step0.2.2.csv', linktext = 'Download Input to Tox Summary with fieldrep column added', verbose = verbose)
  }

  # here we separate the controls from the rest of the samples
  controls <- tox.summary.input %>%
    filter(
      stationid == '0000',
      !(sampletypecode %in% results.sampletypes)
    )

  writelog(
    "\n### Separate the controls from the rest of the samples:\n  ",
    logfile = logfile,
    code = "controls <- tox.summary.input %>% filter( stationid == '0000', !(sampletypecode %in% results.sampletypes) )",
    data = controls,
    verbose = verbose
  )
  create_download_link(data = controls, logfile = logfile, filename = 'tox.summary-controls.csv', linktext = 'Download Tox Summary Controls', verbose = verbose)


  # get rid of the controls. They will be merged later
  results <- tox.summary.input %>%
    filter(
      ( (stationid != '0000') | (sampletypecode %in% results.sampletypes) )
    )

  writelog(
    "\n### Get the results dataframe without the controls\n  ",
    logfile = logfile,
    code = "results <- tox.summary.input %>% filter( ( (stationid != '0000') | (sampletypecode %in% results.sampletypes) ) )",
    data = results,
    verbose = verbose
  )
  create_download_link(data = results, logfile = logfile, filename = 'tox.summary-toxresults.csv', linktext = 'Download Tox Results (No controls)', verbose = verbose)


  # Join results and controls on 'toxbatch', 'species', 'labrep', and 'lab'
  summary <- results %>%
    inner_join(
      controls,
      by = c('toxbatch','species','labrep','lab'),
      suffix = c('','_control')
    )
  writelog(
    "\n### Join results and controls on 'toxbatch', 'species', 'labrep', and 'lab'\n  ",
    logfile = logfile,
    code = "
      summary <- results %>%
        inner_join(
          controls,
          by = c('toxbatch','species','labrep','lab'),
          suffix = c('','_control')
        )
    ",
    data = results %>% head(25),
    verbose = verbose
  )
  create_download_link(data = results, logfile = logfile, filename = 'tox.summary-results-w-controls.csv', linktext = 'Download Tox Results joined with controls', verbose = verbose)

  writelog("\n#### Explanation of what the next steps are going to be",logfile = logfile, verbose = verbose)
  writelog('Group by lab, stationid, toxbatch, species, fieldrep, sampletypecode\n', logfile = logfile, verbose = verbose)
  writelog('Calculate two sided t-test p value (t.test function in R) - the two samples being the result and the control result\n', logfile = logfile, verbose = verbose)
  writelog('t.test function in R by default will treat them as two samples from two populations (not assuming equal variances)\n', logfile = logfile, verbose = verbose)
  writelog('The function call to t.test is this:\n', logfile = logfile, verbose = verbose)
  writelog("> t.test(x = result, y = result_control, mu = 0, var.equal = F, alternative = 'two.sided')$p.value / 2\n  ", logfile = logfile, verbose = verbose)
  writelog("Further explanation of the function call:\n", logfile = logfile, verbose = verbose)
  writelog("> x = result and y = result_control: These are the two vectors (samples) being compared.", logfile = logfile, verbose = verbose)
  writelog("> mu = 0: The null hypothesis is that the difference between the means of the two samples is zero.", logfile = logfile, verbose = verbose)
  writelog("> var.equal = F: The variances of the two samples are not assumed to be equal (Welch's t-test).", logfile = logfile, verbose = verbose)
  writelog("> alternative = 'two.sided': The alternative hypothesis is that the means are different (two-sided test).", logfile = logfile, verbose = verbose)
  writelog("> $p.value / 2: The p-value of the t-test is divided by 2. This division suggests that the intention is to obtain a one-tailed p-value from a two-tailed test.", logfile = logfile, verbose = verbose)

  # Get the stats
  summary <- summary %>%
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
    )
  writelog(
    "\n### Get the tox summary statistics\n  ",
    logfile = logfile,
    code = "
      # Get the stats
      summary <- summary %>%
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
        )
    ",
    data = summary %>% head(25),
    verbose = verbose
  )
  create_download_link(data = summary, logfile = logfile, filename = 'tox.summary-stats.csv', linktext = 'Download Tox Summary (statistics)', verbose = verbose)

  # Write tox categories to the logs
  writelog("#### These are the tox category thresholds that we use for binning into Tox SQO Categories (CASQO Manual June 2021 ed. pages 106-108):", data = tox_categories, logfile = logfile, verbose = verbose)

  # join with tox categories
  summary <- summary %>%
    left_join(
      tox_categories, by = 'species'
    )
  writelog(
    "\n### Join with the Tox Categories dataframe\n  ",
    code = "
      # join with tox categories
      summary <- summary %>%
        left_join(
          tox_categories, by = 'species'
        )
    ",
    data = summary %>% head(25),
    logfile = logfile,
    verbose = verbose
  )
  create_download_link(data = summary, logfile = logfile, filename = 'tox.summary-with-categories.csv', linktext = 'Download Tox Summary joined with the categories dataframe', verbose = verbose)


  # Assign SQO categories
  summary <- summary %>%
    mutate(
      # CASQO Technical Manual page 106-108 (June 2021 Edition)
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
        sqo_category_value == 1 ~ 'Nontoxic',
        sqo_category_value == 2 ~ 'Low Toxicity',
        sqo_category_value == 3 ~ 'Moderate Toxicity',
        sqo_category_value == 4 ~ 'High Toxicity',
        TRUE ~ NA_character_
      )
    )
  writelog(
    "\n### Assign SQO categories\n  ",
    code = "
      summary <- summary %>%
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
            sqo_category_value == 1 ~ 'Nontoxic',
            sqo_category_value == 2 ~ 'Low Toxicity',
            sqo_category_value == 3 ~ 'Moderate Toxicity',
            sqo_category_value == 4 ~ 'High Toxicity',
            TRUE ~ NA_character_
          )
        )
    ",
    data = summary %>% head(25),
    logfile = logfile,
    verbose = verbose
  )
  create_download_link(data = summary, logfile = logfile, filename = 'tox.summary-with-sqo-categories.csv', linktext = 'Download Tox Summary with SQO Categories', verbose = verbose)


  # Finalize selection and naming of columns
  summary <- summary %>%
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

  writelog(
    "\n### Finalize selection and naming of columns\n  ",
    code = "
      # Finalize selection and naming of columns
      summary <- summary %>%
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
    ",
    data = summary %>% head(25),
    logfile = logfile,
    verbose = verbose
  )
  create_download_link(data = summary, logfile = logfile, filename = 'tox.summary-final.csv', linktext = 'Download Final Tox Summary Dataframe', verbose = verbose)


  writelog('\n## END: Tox Summary function.\n', logfile = logfile, verbose = verbose)

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
tox.sqo <- function(toxresults, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = F) {

  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\n# BEGIN: Tox SQO function.\n  ', logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'toxicity.sqo.input.RData'
  if (verbose) {
    save(toxresults, file = file.path( dirname(logfile), rawinput.filename ))
  }


  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n## Raw input to tox.sqo:\n  ",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') # This will load a dataframe called 'toxresults' into your environment"),
    data = toxresults %>% head(25),
    verbose = verbose
  )
  create_download_link(data = toxresults, logfile = logfile, filename = 'tox.sqo-initial-input-step0.csv', linktext = 'Download Raw Input to Tox SQO Function', verbose = verbose)

  # It is common that there is a station called "0" which is actually supposed to be a string of four zeros ('0000')
  toxresults$stationid <- replace(toxresults$stationid, toxresults$stationid %>% as.character() == '0', '0000')

  writelog(
    "\n## Replace stations that may say '0' with the bight 'QA station' '0000'\n  ",
    logfile = logfile,
    code = "
      # It is common that there is a station called '0' which is actually supposed to be a string of four zeros ('0000')
      toxresults$stationid <- replace(toxresults$stationid, toxresults$stationid %>% as.character() == '0', '0000')
    ",
    verbose = verbose
  )
  create_download_link(
    data = toxresults %>% head(25),
    logfile = logfile,
    filename = 'tox.sqo-initial-input-step0.1.csv',
    linktext = 'Download ToxResults after replacing 0 with 0000 in stationid column',
    verbose = verbose
  )

  # Call tox summary function
  tox_nonintegrated1 <- tox.summary(toxresults, logfile = logfile, verbose = verbose)
  writelog(
    "\n## Calling Tox Summary within Tox.SQO function\n  ",
    logfile = logfile,
    code = "
      # Call tox summary function
      tox_nonintegrated1 <- tox.summary(toxresults)
    ",
    data = tox_nonintegrated1 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_nonintegrated1, logfile = logfile, filename = 'tox.summary-output.csv', linktext = 'Download Tox Summary function initial output', verbose = verbose)

  # Select only certain columns
  tox_nonintegrated2 <- tox_nonintegrated1 %>%
    select(
      stationid,
      species,
      `Endpoint Method`,
      Score,
      Category
    )
  writelog(
    "\n## Selecting certain columns stationid, species, `Endpoint Method`, Score, Category\n  ",
    logfile = logfile,
    code = "
      tox_nonintegrated2 <- tox_nonintegrated1 %>%
      select(
        stationid,
        species,
        `Endpoint Method`,
        Score,
        Category
      )
    ",
    data = tox_nonintegrated2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_nonintegrated2, logfile = logfile, filename = 'tox.summary-output-selectcolumns.csv', linktext = 'Download Tox Summary output after selecting certain columns', verbose = verbose)

  # Isolate Genus
  tox_nonintegrated3 <- tox_nonintegrated2 %>%
    separate(
      species,
      c("Genus","Species")
    ) %>%
    select(-Species)
  writelog(
    "\n## Isolate Genus - separate Genus from the Species\n  ",
    logfile = logfile,
    code = "
      # Isolate the Genus
      tox_nonintegrated3 <- tox_nonintegrated2 %>%
        separate(
          species,
          c('Genus','Species')
        ) %>%
        select(-Species)
    ",
    data = tox_nonintegrated3 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_nonintegrated3, logfile = logfile, filename = 'tox.summary-output-isolate-genus.csv', linktext = 'Download Tox Summary output after Genus was isolated', verbose = verbose)

  # Define the Type of test based on Genus and Endpoint Method
  tox_nonintegrated4 <- tox_nonintegrated3 %>%
    mutate(
      Index = paste(Genus, `Endpoint Method`),
      `Test Type` = case_when(
        `Endpoint Method` %in% c('Growth','NormDev') ~ 'Sublethal',
        `Endpoint Method` == 'Survival' ~ 'Acute',
        TRUE ~ NA_character_
      )
    )
  writelog(
    "\n## Isolate Genus - separate Genus from the Species\n  ",
    logfile = logfile,
    code = "
      # Define the Type of test based on Genus and Endpoint Method
      tox_nonintegrated4 <- tox_nonintegrated3 %>%
        mutate(
          Index = paste(Genus, `Endpoint Method`),
          `Test Type` = case_when(
            `Endpoint Method` %in% c('Growth','NormDev') ~ 'Sublethal',
            `Endpoint Method` == 'Survival' ~ 'Acute',
            TRUE ~ NA_character_
          )
        )
    ",
    data = tox_nonintegrated4 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_nonintegrated4, logfile = logfile, filename = 'tox.summary-output-define-testtype.csv', linktext = 'Download Tox Summary output after defining test type', verbose = verbose)


  # Create category score column for purpose of final unified output
  tox_nonintegrated <- tox_nonintegrated4 %>%
    select(stationid, Index, `Test Type`, Score, Category) %>%
    mutate(
      `Category Score` = Score # just for purposes of the very final unified output, all three LOE's in one table
    )
  writelog(
    "\n## Create Category Score column\n  ",
    logfile = logfile,
    code = "
      # Create category score column for purpose of final unified output
      tox_nonintegrated <- tox_nonintegrated4 %>%
        select(stationid, Index, `Test Type`, Score, Category) %>%
        mutate(
          `Category Score` = Score # just for purposes of the very final unified output, all three LOE's in one table
        )
    ",
    data = tox_nonintegrated %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_nonintegrated, logfile = logfile, filename = 'tox.summary-withcatscore-column.csv', linktext = 'Download Tox Summary output with Category Score Column initialized', verbose = verbose)


  if ('stratum' %in% names(toxresults)) {

    # Tack back on the stratum column if it was initially included
    tox_nonintegrated <- tox_nonintegrated %>%
      left_join(
        toxresults %>% select(stationid, stratum) %>% distinct() %>% filter(!is.na(stratum)),
        by = 'stationid'
      )
    writelog(
      "\n## Include the stratum if it was provided in the initial toxresults\n  ",
      logfile = logfile,
      code = "
      # Include the stratum if it was provided in the initial toxresults
      tox_nonintegrated <- tox_nonintegrated %>%
        left_join(
          toxresults %>% select(stationid, stratum) %>% distinct() %>% filter(!is.na(stratum)),
          by = 'stationid'
        )
    ",
      data = tox_nonintegrated %>% head(25),
      verbose = verbose
    )
    create_download_link(data = tox_nonintegrated, logfile = logfile, filename = 'tox.summary-withstratum.csv', linktext = 'Download Tox Summary output with stratum column', verbose = verbose)


  }

  writelog("### Explanation of assigning the Tox SQO Score based on the Manual", logfile = logfile, verbose = verbose)

  writelog("For Toxicity, we take the mean of SQO scores among the tox tests at the site (CASQO manual, June 2021 edition, page 109).\n--\n", logfile = logfile, verbose = verbose)

  writelog("The score is the mean of the two, however, if one is missing, the score should come out to NA", logfile = logfile, verbose = verbose)
  writelog("This means that in this particular case, it should be na.rm = F", logfile = logfile, verbose = verbose)
  writelog("**(NOTE: Not for offshore sites - offhsore sites only require the Amphipod test)**\n--\n", logfile = logfile, verbose = verbose)

  writelog("This is not explicitly mentioned on page 109 of the manual, where the instructions are given for integrating the results of the tox tests to determine the Tox LOE.\n", logfile = logfile, verbose = verbose)
  writelog("However, we can understand that it must be na.rm = F because there are essentially two types of test methods, Acute and Sublethal (Page 85) and in the intro ", logfile = logfile, verbose = verbose)
  writelog("to the instructions for tox assessment on page 84, first paragraph, it says the following:\n", logfile = logfile, verbose = verbose)

  writelog("> Multiple toxicity tests are needed to assess toxicity because no single method exists that can capture the ", logfile = logfile, verbose = verbose)
  writelog("> full spectrum of potential contaminant effects. ", logfile = logfile, verbose = verbose)
  writelog("> Toxicity assessment under the California ", logfile = logfile, verbose = verbose)
  writelog("> Sediment Quality Objectives (CASQO) framework requires information from two types of tests: ", logfile = logfile, verbose = verbose)
  writelog("> 1) short-term amphipod survival and \n2) a sublethal test.\n\n  ", logfile = logfile, verbose = verbose)

  writelog("For this reason, if one of the test results is missing, we will declare the Tox LOE to be indeterminate or unknown. ", logfile = logfile, verbose = verbose)
  writelog("So here, we will actually need to check if both tests were present. ", logfile = logfile, verbose = verbose)
  writelog("Acute tests mean the endpoint method is 'Survival' and the Sublethal tests mean the endpoint is 'Growth' or 'NormDev' (Normal Development)\n--\n", logfile = logfile, verbose = verbose)



  tox_integrated0 <- tox_nonintegrated %>%
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
    )
  writelog(
    "\n## Determine if all necessary tests are present and assign a score for the site\n  ",
    logfile = logfile,
    code = "
      tox_integrated0 <- tox_nonintegrated %>%
        group_by(stationid) %>%
        summarize(
          has.all.tests = all(c('Acute','Sublethal') %in% `Test Type`),
          offshore = 'stratum' %in% names(.) && stratum %in% c('Outer Shelf','Mid Shelf','Inner Shelf','Channel Islands'),
          Score = case_when(
            # offhsore sites only require the Amphipod test - per Ken Schiff August 1, 2022
            offshore ~ ceiling(mean(Score, na.rm = T)),
            has.all.tests ~ ceiling(mean(Score, na.rm = F)),
            TRUE ~ NA_real_
          )
        )
    ",
    data = tox_integrated0 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_integrated0, logfile = logfile, filename = 'tox.integrated.score.initial.csv', linktext = 'Download Initial step of determining the integrated tox score', verbose = verbose)


  # Assign Category and select final columns
  tox_integrated <- tox_integrated0 %>%
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

  writelog(
    "\n## Assign Category and select final columns\n  ",
    logfile = logfile,
    code = '
      # Assign Category and select final columns
      tox_integrated <- tox_integrated0 %>%
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
          `Category Score` = Score # just for purposes of the very final unified output, all three LOEs in one table
        ) %>%
        select(stationid, Index, Score, Category, `Category Score`)
    ',
    data = tox_integrated %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_integrated, logfile = logfile, filename = 'tox.integrated.scores.csv', linktext = 'Download Integrated Tox Scores for each site', verbose = verbose)

  # final output - integrated scores and scores of individual tests (non integrated scores)
  tox_final <- rbind.fill(tox_integrated, tox_nonintegrated) %>%
    rename(StationID = stationid) %>%
    arrange(
      StationID, Index, Category
    )
  writelog(
    "\n## Concat Integrated scores with the regular non integrated scores (scores for each test)\n  ",
    logfile = logfile,
    code = '
      # final output - integrated scores and scores of individual tests (non integrated scores)
      tox_final <- plyr::rbind.fill(tox_integrated, tox_nonintegrated) %>%
        rename(StationID = stationid) %>%
        arrange(
          StationID, Index, Category
        )
    ',
    data = tox_final %>% head(25),
    verbose = verbose
  )
  create_download_link(data = tox_final, logfile = logfile, filename = 'tox.sqo.final.csv', linktext = 'Download Final Tox SQO Data', verbose = verbose)


  writelog('\n# END: Tox SQO function.\n', logfile = logfile, verbose = verbose)


  return(tox_final)
}
