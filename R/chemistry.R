# ------------ LRM (Logistic Regression Model)--------------
#' Logistic Regression Model Score
#'
#' This function calculates the Logistic Regression Model score
#' in order to calculate the Chemistry SQO scores for the given sites.
#' The ultimate guide for this function is the CASQO Technical Manual Chapter 3
#' (end of page 31 to beginning of page 34)
#'
#' @usage LRM(chemdata)
#'
#' @param chemdata a dataframe with the following headings:
#'
#'    \code{StationID},
#'
#'    \code{AnalyteName},
#'
#'    \code{Result},
#'
#'    \code{RL},
#'
#'    \code{MDL}
#'
#' @examples
#' data(chem_sampledata) # load sample data to your environment
#' LRM(chem_sampledata) # get scores and see output
#'
#' @import dplyr
#' @export
LRM <- function(chemdata.lrm.input, preprocessed = F, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-chemlog.Rmd') ), verbose = T)  {
  "lrm_table"

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\n\n## Chem CA LRM Function\n", logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'lrm.input.RData'
  if (verbose) {
    save(chemdata.lrm.input, file = file.path( dirname(logfile), rawinput.filename ))
  }

  writelog(
    "\n### Initial input to LRM:",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "')"),
    data = chemdata.lrm.input,
    verbose = verbose
  )
  create_download_link(data = chemdata.lrm.input, logfile = logfile, filename = 'LRM_InitialInputData.csv', linktext = 'Download Initial Input to Chem LRM Function', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  #  ---- Make the call to the Preprocessing function (If not already preprocessed) ----
  if (!preprocessed) {
    # Write to the log file
    writelog(
      "\n#### Preprocessing chemistry data (Details within chemdata_prep to be shown later as well):\n",
      logfile = logfile,
      verbose = verbose
    )

    # Log the separation space between steps
    writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)


    # Actually call the function
    chemdata.lrm.input <- chemdata_prep(chemdata.lrm.input, logfile = logfile, verbose = verbose)


    # Log the separation space between steps
    writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)

    # Create code block and download link to the preprocessed data
    writelog(
      "\n#### Chemdata Pre processing function is finished executing - Here is its final output along with a code block (for R Studio users):",
      logfile = logfile,
      code = 'chemdata.lrm.input <- chemdata_prep(chemdata.lrm.input, verbose = FALSE)',
      verbose = verbose
    )
    create_download_link(data = chemdata.lrm.input, logfile = logfile, filename = 'LRM_PreProcessedInput.csv', linktext = 'Download Preprocessed Input to LRM Function (after calling chemdata_prep)', verbose = verbose)

  }


  # ---- Display the LRM Table ----
  writelog(
    "### LRM table by AnalyteName (Table 3.5 on page 39 of Technical Manual)",
    data = lrm_table,
    logfile = logfile,
    verbose = verbose,
    pageLength = 12
  )


  # ---- CA LRM Step 0 - Take Log10 of the results ----

  # Write to log file
  writelog("\n### Step 0: Take the Log10 of the chemistry concentration.", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata.log10 <- chemdata.lrm.input %>%
    mutate(
      logResult = log10(Result)
    )

  # Write to code portion of the logs
  writelog(
    "",
    code = '
      chemdata.log10 <- chemdata.lrm.input %>%
        mutate(
          logResult = log10(Result)
        )
    ',
    data = chemdata.log10,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata.log10, logfile = logfile, filename = 'LRM_Step0 (Log Transform).csv', linktext = 'Download Log Transformed Input to LRM Function (Step 0)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  writelog("\n### Manipulate the LRM chem data per Technical Manual page 37-40", logfile = logfile, verbose = verbose)


  # -------- Step 1 - Join LRM Table with the log transformed chemistry results data --------
  # Write to log
  writelog("\n#### Step 1", logfile = logfile, verbose = verbose)
  writelog("\n##### Join with the LRM table by AnalyteName (Table 3.5 on page 39 of Technical Manual)", logfile = logfile, verbose = verbose)

  # Write the table 3.6 (Since it is called table 3.6 in the SQO Technical Manual 3rd Edition) to the log file
  table3.6 <- "
    | Category          | Range           | Category Score |
    |-------------------|-----------------|----------------|
    | Minimal Exposure  | < 0.33          | 1              |
    | Low Exposure      | ≥ 0.33 - ≤ 0.49 | 2              |
    | Moderate Exposure | > 0.49 - ≤ 0.66 | 3              |
    | High Exposure     | > 0.66          | 4              |
  "
  writelog("##### Here is the table:", logfile = logfile, verbose = verbose)
  writelog(table3.6, logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_lrm1 <- lrm_table %>%
    left_join(chemdata.log10, by = 'AnalyteName')

  # Write to code portion of the logs
  writelog(
    "",
    code = '
      chemdata_lrm1 <- lrm_table %>%
        left_join(chemdata.log10, by = \'AnalyteName\')
    ',
    data = chemdata_lrm1,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_lrm1, logfile = logfile, filename = 'LRM_Step1.csv', linktext = 'Download Step 1 LRM Data', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)



  # -------- Step 2 - Get P Values for each analyte of each station --------
  # Write to log
  writelog("\n#### Step 2", logfile = logfile, verbose = verbose)
  writelog("\n##### p = (exp(B0 + B1 * logResult) / (1 + exp(B0 + B1 * logResult))) ---> Round to 2 decimal places", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_lrm2 <- chemdata_lrm1 %>%
    mutate(
      # page 38 of Technical Manual 3rd Edition June 2021
      p = (exp(B0 + B1 * logResult) / (1 + exp(B0 + B1 * logResult))) %>% round(2)
    )

  # Write to code portion of the logs
  writelog(
    "",
    code = '
      chemdata_lrm2 <- chemdata_lrm1 %>%
        mutate(
          # page 38 of Technical Manual 3rd Edition June 2021
          p = (exp(B0 + B1 * logResult) / (1 + exp(B0 + B1 * logResult))) %>% round(2)
        )
    ',
    data = chemdata_lrm2,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_lrm2, logfile = logfile, filename = 'LRM_Step2.csv', linktext = 'Download Step 2 LRM Data', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n____" , logfile = logfile, verbose = verbose)




  # -------- Step 3 - Get Max P Value --------
  writelog("\n#### Step 3", logfile = logfile, verbose = verbose)
  writelog("\n##### Group by StationID (this function and R package overall assumes you are feeding it data for one sediment sample per stationid)", logfile = logfile, verbose = verbose)
  writelog("\n##### take max value of the variable p within each grouping (grouped by StationID)", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_lrm3 <- chemdata_lrm2 %>%
    # Page 39 of Technical Manual
    # Get the max of the "p" values
    group_by(StationID) %>%
    summarize(
      Score = if_else(
        all(is.na(p)), NA_real_, suppressWarnings(max(p, na.rm = T))
      )
    ) %>%
    ungroup()


  # Write to code portion of the logs
  writelog(
    "",
    code = '
    chemdata_lrm3 <- chemdata_lrm2 %>%
      # Page 39 of Technical Manual
      # Get the max of the "p" values
      group_by(StationID) %>%
      summarize(
        Score = if_else(
          all(is.na(p)), NA_real_, suppressWarnings(max(p, na.rm = T))
        )
      ) %>%
      ungroup()
    ',
    data = chemdata_lrm3,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_lrm3, logfile = logfile, filename = 'LRM_Step3.csv', linktext = 'Download Step 3 LRM Data', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)




  # ---- Main Intermediate CA LRM Calculation QA Step ----

  writelog("\n#### *Main Intermediate CA LRM Calculation QA Step*\n", logfile = logfile, verbose = verbose)
  writelog("\n- Just basically selecting certain columns, and 'R Binding' the pmax value as 'SAMPLE_PMAX' for the whole station\n", logfile = logfile, verbose = verbose)
  writelog("\n**The ouput from this intermediate calculation step is designed to make an easier comparison with output from the Excel Tool**\n", logfile = logfile, verbose = verbose)

  # Actually execute the code
  lrm.main.intermediate.qa.calc.step <- chemdata_lrm2 %>%
    select(
      StationID,
      CAP_parameter = AnalyteName,
      Value = p
    ) %>%
    mutate(
      CAP_parameter = paste0(
        'p_',
        case_when(
          CAP_parameter == 'alpha-Chlordane' ~ 'CHLORDAN_A',
          CAP_parameter == 'trans-Nonachlor' ~ 'NONACHL_TR',
          CAP_parameter == 'PCBs_total' ~ 'PCB_SUM',
          CAP_parameter == "4,4'-DDT" ~ "4,4'_DDT",
          TRUE ~ toupper(CAP_parameter)
        )
      )
    ) %>%
    rbind(
      chemdata_lrm3 %>%
        mutate(
          CAP_parameter = 'SAMPLE_PMAX'
        ) %>%
        select(
          StationID,
          CAP_parameter,
          Value = Score
        )
    ) %>%
    arrange(StationID, CAP_parameter)

  # Write code portion of the logs
  writelog(
    "",
    code = "
      lrm.main.intermediate.qa.calc.step <- chemdata_lrm2 %>%
        select(
          StationID,
          CAP_parameter = AnalyteName,
          Value = p
        ) %>%
        mutate(
          CAP_parameter = paste0(
            'p_',
            case_when(
              CAP_parameter == 'alpha-Chlordane' ~ 'CHLORDAN_A',
              CAP_parameter == 'trans-Nonachlor' ~ 'NONACHL_TR',
              CAP_parameter == 'PCBs_total' ~ 'PCB_SUM',
              CAP_parameter == \"4,4'-DDT\" ~ \"4,4'_DDT\",
              TRUE ~ toupper(CAP_parameter)
            )
          )
        ) %>%
        rbind(
          chemdata_lrm3 %>%
            mutate(
              CAP_parameter = 'SAMPLE_PMAX'
            ) %>%
            select(
              StationID,
              CAP_parameter,
              Value = Score
            )
        ) %>%
        arrange(StationID, CAP_parameter)
    ",
    data = lrm.main.intermediate.qa.calc.step,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = lrm.main.intermediate.qa.calc.step, logfile = logfile, filename = 'LRM_IntermediateCalcQA.csv', linktext = 'Download Main LRM Intermediate Calculation QA Data', verbose = verbose)





  # Log the separation space between steps
  writelog("\n___\n___\n___\n____" , logfile = logfile, verbose = verbose)






  # -------- Step 4 -  Assign the category and category score --------

  # Write to the logs
  writelog("\n#### Step 4", logfile = logfile, verbose = verbose)
  writelog("\n##### Assign LRM Categories based on thresholds defined in table 3.6 of Technical Manual page 40", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_lrm.final <- chemdata_lrm3 %>%
    mutate(
      # Page 40 of Technical Manual (3rd edition June 2021)
      `Category Score` = case_when(
        Score < 0.33 ~ 1,
        Score >= 0.33 & Score <= 0.49 ~ 2,
        Score > 0.49 & Score <= 0.66 ~ 3,
        Score > 0.66 ~ 4,
        TRUE ~ NA_real_
      ),
      Category = case_when(
        `Category Score` == 1 ~ "Minimal Exposure",
        `Category Score` == 2 ~ "Low Exposure",
        `Category Score` == 3 ~ "Moderate Exposure",
        `Category Score` == 4 ~ "High Exposure",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      Index = 'LRM'
    ) %>%
    select(
      StationID, Index, Score, Category, `Category Score`
    )


  # Write to code portion of the logs
  writelog(
    "",
    code = '
      chemdata_lrm.final <- chemdata_lrm3 %>%
        mutate(
          # Page 40 of Technical Manual (3rd edition June 2021)
          `Category Score` = case_when(
            Score < 0.33 ~ 1,
            Score >= 0.33 & Score <= 0.49 ~ 2,
            Score > 0.49 & Score <= 0.66 ~ 3,
            Score > 0.66 ~ 4,
            TRUE ~ NA_real_
          ),
          Category = case_when(
            `Category Score` == 1 ~ "Minimal Exposure",
            `Category Score` == 2 ~ "Low Exposure",
            `Category Score` == 3 ~ "Moderate Exposure",
            `Category Score` == 4 ~ "High Exposure",
            TRUE ~ NA_character_
          )
        ) %>%
        mutate(
          Index = \'LRM\'
        ) %>%
        select(
          StationID, Index, Score, Category, `Category Score`
        )
    ',
    data = chemdata_lrm.final,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_lrm.final, logfile = logfile, filename = 'LRM_Final.csv', linktext = 'Download Final LRM Data (Within LRM Function)', verbose = verbose)

  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)

  writelog("##### END Chem LRM Function\n", logfile = logfile, verbose = verbose)

  return(chemdata_lrm.final)

# uncomment to make the above block a function again
}


# ------------ CSI (Chemical Score Index) --------------
#' Chemical Score Index
#'
#' This function calculates the Chemical Score Index for the given sites.
#' The ultimate guide for this function is the CASQO Technical Manual Chapter 3
#' (beginning of page 34 to the beginning of page 36)
#'
#' @usage CSI(chemdata)
#'
#' @param chemdata a dataframe with the following headings:
#'
#'    \code{StationID},
#'
#'    \code{AnalyteName},
#'
#'    \code{Result},
#'
#'    \code{RL},
#'
#'    \code{MDL}
#'
#' @examples
#' data(chem_sampledata) # load sample data to your environment
#' CSI(chem_sampledata) # get scores and see output
#'
#' @export

# Conversation with Darrin on April 16th 2020
# We need to fix the rounding in this thing.
# The fixes may also need to be implemented in the chemdata_prep function
# Certain analytes result values need to be rounded to different numbers of decimal places
# going off memory here. Darrin is the one to consult, and he can show the latest greatest version of the excel tool
# Copper - 1
# Lead  - 1
# Mercury - 2
# Zinc - 1
# HPAH - 1
# LPAH - 1
# All the rest - 2
CSI <- function(chemdata.csi.input, preprocessed = F, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-chemlog.Rmd') ), verbose = T) {
  "csi_weight"

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\n\n## Chem CSI Function\n", logfile = logfile, verbose = verbose)



  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'csi.input.RData'
  if (verbose) {
    save(chemdata.csi.input, file = file.path( dirname(logfile), rawinput.filename ))
  }

  writelog(
    "\n### Initial input to CSI:",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "')"),
    data = chemdata.csi.input,
    verbose = verbose
  )
  create_download_link(data = chemdata.csi.input, logfile = logfile, filename = 'CSI_InitialInputData.csv', linktext = 'Download Initial Input to Chem CSI Function', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  #  ---- Make the call to the Preprocessing function (If not already preprocessed) ----
  if (!preprocessed) {
    # Write to the log file
    writelog(
      "\n#### Preprocessing chemistry data (Details within chemdata_prep to be shown later as well):\n",
      logfile = logfile,
      verbose = verbose
    )

    # Log the separation space between steps
    writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)

    # Actually call the function
    chemdata.csi.input <- chemdata_prep(chemdata.csi.input, logfile = logfile, verbose = verbose)


    # Log the separation space between steps
    writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)


    # Create code block and download link to the preprocessed data
    writelog(
      "\n#### Chemdata Pre processing function is finished executing - Here is its final output along with a code block (for R Studio users):",
      logfile = logfile,
      code = 'chemdata.csi.input <- chemdata_prep(chemdata.csi.input, verbose = FALSE)',
      verbose = verbose
    )
    create_download_link(data = chemdata.csi.input, logfile = logfile, filename = 'CSI_PreProcessedInput.csv', linktext = 'Download Preprocessed Input to CSI Function (after calling chemdata_prep)', verbose = verbose)

  }

  # ---- CSI Weight DataFrame - Do not include a download option ----
  writelog(
    "\n#### This is the CSI Weight DataFrame (Based on table 3.7 and 3.8 on page 40-41 of the manual - 3rd edition June 2021):",
    logfile = logfile,
    verbose = verbose
  )
  writelog("\n##### (Based on table 3.7 and 3.8 on page 40-41 of the manual - 3rd edition June 2021):", logfile = logfile, verbose = verbose)
  writelog("\n###### (the \"BDC\" columns represent the cutoff points to separate the exposure score bins):", logfile = logfile, verbose = verbose)
  writelog(
    "",
    data = csi_weight,
    logfile = logfile,
    verbose = verbose,
    pageLength = 12
  )

  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)

  # Here we begin to manipulate the CSI chem data per Technical Manual page 40-42
  writelog("\n### Manipulate the CSI chem data per Technical Manual page 40-42\n", logfile = logfile, verbose = verbose)



  # ---- Step 0. Combine CSI Weight values with data based on the compound. Exclude compounds not in CSI calculation. ----

  # Write to the log file
  writelog(
    "\n#### Combine CSI Weight values with data based on the compound. Exclude compounds not in CSI calculation.",
    logfile = logfile,
    verbose = verbose
  )

  # Actually execute the code
  chemdata_csi <- csi_weight %>% left_join(chemdata.csi.input, by = "AnalyteName")

  # Write code portion to the log file
  writelog(
    "",
    code = 'chemdata_csi <- csi_weight %>% left_join(chemdata.csi.input, by = "AnalyteName")',
    data = chemdata_csi,
    logfile = logfile,
    verbose = verbose
  )

  # Serve it up for download
  create_download_link(data = chemdata_csi, logfile = logfile, filename = 'CSI_Step0.csv', linktext = 'Download CSI Step0', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)



  # ---- Step 0.5: Round the result value according to the convention given by Darrin. ----

  # Write to the log file
  writelog("\n#### Round the result value according to the convention.", logfile = logfile, verbose = verbose)
  writelog("\nIf the result value is 0.005 or less, we round to 4 decimal places.  \n", logfile = logfile, verbose = verbose)
  writelog("\nIf the result value is under 10, we round to 2 decimal places.  \n", logfile = logfile, verbose = verbose)
  writelog("\nIf the result value is under 100, we round to 1 decimal place.  \n", logfile = logfile, verbose = verbose)
  writelog("\nIf the result value is 100 or more, we round to the nearest integer  \n", logfile = logfile, verbose = verbose)
  writelog("\nThe weights which we are comparing the results values to were rounded this way, so we are rounding the result values this way as well  \n", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_csi <- chemdata_csi %>%
    mutate(
      Result = case_when(
        # June 10, 2024 - It was decided we will round to 4 decimal places if the result value is less than 1
        Result <= 1 ~ round(Result, 4),
        Result < 10 ~ round(Result, 2),
        Result < 100 ~ round(Result, 1),
        TRUE ~ round(Result)
      )
    )

  # Write code portion to the log file
  writelog(
    "",
    code = '
    chemdata_csi <- chemdata_csi %>%
      mutate(
        Result = case_when(
          # June 10, 2024 - It was decided we will round to 4 decimal places if the result value is less than 1
          Result <= 1 ~ round(Result, 4),
          Result < 10 ~ round(Result, 2),
          Result < 100 ~ round(Result, 1),
          TRUE ~ round(Result)
        )
      )
    ',
    data = chemdata_csi,
    logfile = logfile,
    verbose = verbose
  )

  # Serve it up for download
  create_download_link(data = chemdata_csi, logfile = logfile, filename = 'CSI_Step0-rounded.csv', linktext = 'Download CSI Step0 (rounded result values)', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)





  # ---- CSI Calclulation Step Description ----
  writelog("\n\n#### Calculate CSI\n", logfile = logfile, verbose = verbose)
  writelog("\n##### Step 1:\n", logfile = logfile, verbose = verbose)
  writelog("- Step 1a: Make the CSI weights NA where the Result value is NA\n", logfile = logfile, verbose = verbose)
  writelog("- Step 1b: Assign exposure score according to appropriate category (see Table 3.7 on page 40 of the Technical Manual 3rd edition)\n", logfile = logfile, verbose = verbose)
  writelog("- Step 1c: Assign weighted exposure score (exposure_score x weight)\n", logfile = logfile, verbose = verbose)
  writelog("\n##### Step 2:\n", logfile = logfile, verbose = verbose)
  writelog("- Step 2a: Group by StationID (this function and R package overall assumes you are feeding it data for one sediment sample per stationid)\n", logfile = logfile, verbose = verbose)
  writelog("- Step 2b: Sum the weighted scores and the weights - then take the sum of the weighted scores divided by the sum of the weights. - This is the CSI Score\n", logfile = logfile, verbose = verbose)
  writelog("- Step 2c: Assign CSI Categories based on thresholds defined in table 3.9 of Technical Manual page 42 (3rd edition)\n", logfile = logfile, verbose = verbose)
  table3.9 <- "
| Category          | Range           | Category Score |
|-------------------|-----------------|----------------|
| Minimal Exposure  | < 1.69          | 1              |
| Low Exposure      | ≥ 1.69 - ≤ 2.33 | 2              |
| Moderate Exposure | > 2.33 - ≤ 2.99 | 3              |
| High Exposure     | > 2.99          | 4              |
  "
  writelog(table3.9, logfile = logfile, verbose = verbose)
  writelog("\n##### Step 3:\n", logfile = logfile, verbose = verbose)
  writelog("- Step 3a: Assign Numerical Category Score and SQO Category Name", logfile = logfile, verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)



  # ---- Calculate CSI Step 1 ----

  # Write to the log file
  writelog("\n##### Execute Step 1:\n", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_csi1 <- chemdata_csi %>%
    mutate(
      weight = case_when(
        # We needed this because the weights are only summed for mean CSI if the Result value exists
        # So we have to make it NA_real_ where the Result is NA
        !is.na(Result) ~ weight,
        TRUE ~ NA_real_
      ),
      # Page 40 of Technical Manual (3rd edition June 2021)
      exposure_score = case_when(
        Result <= BDC1 ~ 1,
        (Result > BDC1) & (Result <= BDC2) ~ 2,
        (Result > BDC2) & (Result <= BDC3) ~ 3,
        Result > BDC3 ~ 4,
        TRUE ~ NA_real_
      ),
      # Page 41 of Technical Manual (3rd edition June 2021)
      weighted_score = exposure_score * weight # manual calls this the weighted score (page 41) (3rd edition June 2021)
    )

  # Write code portion of the logs
  writelog(
    "",
    code = '
      chemdata_csi1 <- chemdata_csi %>%
        mutate(
          weight = case_when(
            # We needed this because the weights are only summed for mean CSI if the Result value exists
            # So we have to make it NA_real_ where the Result is NA
            !is.na(Result) ~ weight,
            TRUE ~ NA_real_
          ),
          # Page 40 of Technical Manual (3rd edition June 2021)
          exposure_score = case_when(
            Result <= BDC1 ~ 1,
            (Result > BDC1) & (Result <= BDC2) ~ 2,
            (Result > BDC2) & (Result <= BDC3) ~ 3,
            Result > BDC3 ~ 4,
            TRUE ~ NA_real_
          ),
          # Page 41 of Technical Manual (3rd edition June 2021)
          weighted_score = exposure_score * weight # manual calls this the weighted score (page 41) (3rd edition June 2021)
        )
    ',
    data = chemdata_csi1,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_csi1, logfile = logfile, filename = 'CSI_Step1.csv', linktext = 'Download CSI Step1', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  # ---- Calculate CSI Step 2 ----

  # Write to the log file
  writelog("\n##### Execute Step 2:\n", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_csi2 <- chemdata_csi1 %>%
    group_by(
      # We would need to also group by SampleDate - although typically calculating SQO for Bight data, we dont need a sampledate
      StationID
    ) %>%
    summarize(
      # Page 41 of Technical Manual (Chart)
      weighted_score_sum = sum(weighted_score, na.rm = T),
      weight = sum(weight, na.rm = T),
      Score = round(weighted_score_sum / weight, 2)
    ) %>%
    ungroup()

  # Write code portion of the logs
  writelog(
    "",
    code = '
      chemdata_csi2 <- chemdata_csi1 %>%
        group_by(
          # We would need to also group by SampleDate - although typically calculating SQO for Bight data, we dont need a sampledate
          StationID
        ) %>%
        summarize(
          # Page 41 of Technical Manual 3rd edition (Chart)
          weighted_score_sum = sum(weighted_score, na.rm = T),
          weight = sum(weight, na.rm = T),
          Score = round(weighted_score_sum / weight, 2)
        ) %>%
        ungroup()
    ',
    data = chemdata_csi2,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_csi2, logfile = logfile, filename = 'CSI_Step2.csv', linktext = 'Download CSI Step2', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___\n____" , logfile = logfile, verbose = verbose)




  # ---- Main Intermediate CSI Calculation QA Step ----

  writelog("\n#### *Main Intermediate CSI Calculation QA Step*\n", logfile = logfile, verbose = verbose)
  writelog("\n- Join output from Step 1 and 2\n", logfile = logfile, verbose = verbose)
  writelog("\n- (step1 %>% left_join step2)\n", logfile = logfile, verbose = verbose)
  writelog("\n**The ouput from this intermediate calculation step is designed to make an easier comparison with output from the Excel Tool**\n", logfile = logfile, verbose = verbose)

  # Actually execute the code
  csi.main.intermediate.qa.calc.step <- chemdata_csi1 %>%
    select(
      StationID,
      Chemical = AnalyteName,
      Weight = weight,
      Category = exposure_score,
      CCS = weighted_score
    ) %>%
    left_join(
      chemdata_csi2 %>%
        select(
          StationID,
          `Sum CCS` = weighted_score_sum,
          `Sum Weight` = weight,
          `Weighted Mean` = Score
        )
      ,
      by = 'StationID'
    ) %>%
    mutate(
      Chemical = case_when(
        Chemical == 'alpha-Chlordane' ~ 'Alpha Chlordane',
        Chemical == 'gamma-Chlordane' ~ 'Gamma Chlordane',
        Chemical == 'DDDs_total' ~ 'DDDs, total',
        Chemical == 'DDEs_total' ~ 'DDEs, total',
        Chemical == 'DDTs_total' ~ 'DDTs, total',
        Chemical == 'PCBs_total' ~ 'PCBs, total',
        TRUE ~ Chemical
      )
    ) %>%
    arrange(StationID, Chemical)

  # Write code portion of the logs
  writelog(
    "",
    code = '
      csi.main.intermediate.qa.calc.step <- chemdata_csi1 %>%
        select(
          StationID,
          Chemical = AnalyteName,
          Weight = weight,
          Category = exposure_score,
          CCS = weighted_score
        ) %>%
        left_join(
          chemdata_csi2 %>%
            select(
              StationID,
              `Sum CCS` = weighted_score_sum,
              `Sum Weight` = weight,
              `Weighted Mean` = Score
            )
          ,
          by = \'StationID\'
        ) %>%
        mutate(
          Chemical = case_when(
            Chemical == \'alpha-Chlordane\' ~ \'Alpha Chlordane\',
            Chemical == \'gamma-Chlordane\' ~ \'Gamma Chlordane\',
            Chemical == \'DDDs_total\' ~ \'DDDs, total\',
            Chemical == \'DDEs_total\' ~ \'DDEs, total\',
            Chemical == \'DDTs_total\' ~ \'DDTs, total\',
            Chemical == \'PCBs_total\' ~ \'PCBs, total\',
            TRUE ~ Chemical
          )
        ) %>%
        arrange(StationID, Chemical)
    ',
    data = csi.main.intermediate.qa.calc.step,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = csi.main.intermediate.qa.calc.step, logfile = logfile, filename = 'CSI_IntermediateCalcQA.csv', linktext = 'Download Main CSI Intermediate Calculation QA Data', verbose = verbose)




  # Log the separation space between steps
  writelog("\n___\n___\n___\n____" , logfile = logfile, verbose = verbose)




  # ---- Calculate CSI Step 3 ----
  # write to the logs
  writelog("\n##### Execute Step 3:\n", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_csi.final <- chemdata_csi2 %>%
    # Remove unnecessary columns - weighted_score_sum and weight
    select(-c(weighted_score_sum, weight)) %>%

    # Assign category based on the bins specified in page 42 of Technical Manual (3rd edition from June 2021)
    mutate(
      `Category Score` = case_when(
        Score < 1.69 ~ 1,
        (Score >= 1.69) & (Score <= 2.33) ~ 2,
        (Score >= 2.33) & (Score <= 2.99) ~ 3,
        Score > 2.99 ~ 4,
        TRUE ~ NA_real_
      ),
      Category = case_when(
        `Category Score` == 1 ~ "Minimal Exposure",
        `Category Score` == 2 ~ "Low Exposure",
        `Category Score` == 3 ~ "Moderate Exposure",
        `Category Score` == 4 ~ "High Exposure",
        TRUE ~ NA_character_
      )
    ) %>%

    # To identify that the score is for the Chemical Score Index
    mutate(
      Index = "CSI"
    ) %>%

    # Ordering the columns
    select(
      StationID, Index, Score, Category, `Category Score`
    )

  # write to the code portion of the logs - provide final data for download
  writelog(
    "",
    code = '
      chemdata_csi.final <- chemdata_csi2 %>%
        # Remove unnecessary columns - weighted_score_sum and weight
        select(-c(weighted_score_sum, weight)) %>%

        # Assign category based on the bins specified in page 42 of Technical Manual (3rd edition from June 2021)
        mutate(
          `Category Score` = case_when(
            Score < 1.69 ~ 1,
            (Score >= 1.69) & (Score <= 2.33) ~ 2,
            (Score >= 2.33) & (Score <= 2.99) ~ 3,
            Score > 2.99 ~ 4,
            TRUE ~ NA_real_
          ),
          Category = case_when(
            `Category Score` == 1 ~ "Minimal Exposure",
            `Category Score` == 2 ~ "Low Exposure",
            `Category Score` == 3 ~ "Moderate Exposure",
            `Category Score` == 4 ~ "High Exposure",
            TRUE ~ NA_character_
          )
        ) %>%

        # To identify that the score is for the Chemical Score Index
        mutate(
          Index = "CSI"
        ) %>%

        # Ordering the columns
        select(
          StationID, Index, Score, Category, `Category Score`
        )
    ',
    data = chemdata_csi.final,
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(data = chemdata_csi.final, logfile = logfile, filename = 'CSI_Final.csv', linktext = 'Download CSI Final (Final DataFrame within CSI Function)', verbose = verbose)

  return(chemdata_csi.final)
}

# ------------ Integrated Score --------------
#' Get Integrated Chem SQO Scores and Categories
#'
#' This function will not only calculate CSI and LRM, but it
#' will also get the overall integrated SQO Category and Score
#' for the stations that are given to it.
#'
#' @usage chem.sqo(chemdata)
#'
#' @param chemdata a dataframe with the following columns:
#'
#'    \code{StationID},
#'
#'    \code{AnalyteName},
#'
#'    \code{Result},
#'
#'    \code{RL},
#'
#'    \code{MDL}
#'
#' @examples
#' data(chem_sampledata) # load sample data to your environment
#' chem.sqo(chem_sampledata) # get scores and see output
#'
#' @export
chem.sqo <- function(chemdata, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'chemlog.Rmd' ), verbose = T, logtitle = 'Chemistry SQO Logs') {

  # ---- Initialize Logging ----
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose, title = logtitle)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\n# Chemistry SQO Main Function\n", logfile = logfile, verbose = verbose)


  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'chem.sqo.input.RData'
  if (verbose) {
    save(chemdata, file = file.path( dirname(logfile), rawinput.filename ))
  }

  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n## Raw input to chem.sqo:",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') ### This will load a dataframe called 'chemdata' into your environment"),
    verbose = verbose
  )
  create_download_link(data = chemdata, logfile = logfile, filename = 'ChemSQO-RawInput.csv', linktext = 'Download Raw Input to Chem SQO Function', verbose = verbose)



  #  ---- Make the call to the Preprocessing function ----
  # Write to the log file
  writelog(
    "\n## Preprocessing chemistry data (Details within chemdata_prep to be shown later as well):\n",
    logfile = logfile,
    verbose = verbose
  )

  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)


  # Actually call the function
  chemdata <- chemdata_prep(chemdata, logfile = logfile, verbose = verbose)

  # Create code block and download link to the preprocessed data
  writelog(
    "\n## Chemdata Pre processing function is finished executing - Here is its final output along with a code block (for R Studio users):",
    logfile = logfile,
    code = 'chemdata <- chemdata_prep(chemdata, verbose = FALSE)',
    data = chemdata,
    verbose = verbose
  )
  create_download_link(data = chemdata, logfile = logfile, filename = 'ChemSQO-PreProcessedInput.csv', linktext = 'Download Preprocessed Input to Chem SQO Function (after calling chemdata_prep)', verbose = verbose)

  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)


  #  ---- Make the call to the LRM function ----
  # Write to the log file
  writelog(
    "\n## Call LRM function within chem.sqo....\n",
    logfile = logfile,
    verbose = verbose
  )
  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)

  # Actually call it
  chemdata_lrm <- LRM(chemdata, preprocessed = T, logfile = logfile, verbose = verbose)

  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)

  # Create code block and download link
  writelog(
    "\n## LRM Function finished executing",
    logfile = logfile,
    verbose = verbose
  )
  writelog(
    "#### Here is its final output along with a code block (for R Studio users):",
    logfile = logfile,
    code = 'chemdata_lrm <- LRM(chemdata, preprocessed = TRUE,  verbose = FALSE)',
    data = chemdata_lrm,
    verbose = verbose
  )
  create_download_link(data = chemdata_lrm, logfile = logfile, filename = 'ChemSQO-LRM-Output.csv', linktext = 'Download Final LRM Output', verbose = verbose)



  #  ---- Make the call to the CSI function ----
  # Write to the log file
  writelog(
    "\n## Call CSI function within chem.sqo\n",
    logfile = logfile,
    verbose = verbose
  )

  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)

  # Actually call it
  chemdata_csi <- CSI(chemdata, preprocessed = T, logfile = logfile, verbose = verbose)

  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)

  # Create code block and download link
  writelog(
    "\n## CSI Function finished executing",
    logfile = logfile,
    verbose = verbose
  )
  writelog(
    "\n### Here is its final output along with a code block (for R Studio users):",
    logfile = logfile,
    code = 'chemdata_csi <- CSI(chemdata, preprocessed = TRUE, verbose = FALSE)',
    data = chemdata_csi,
    verbose = verbose
  )
  create_download_link(data = chemdata_csi, logfile = logfile, filename = 'ChemSQO-CSI-Output.csv', linktext = 'Download Final CSI Output', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)


  # ---- Subsequent steps after CSI and LRM ----
  writelog(
    "\n## Subsequent Steps After CSI and LRM",
    logfile = logfile,
    verbose = verbose
  )


  # ---- First Subsequent Step: Rbind the CSI and LRM Dataframes ----

  # Actually execute the code
  combined1 <- rbind(chemdata_lrm, chemdata_csi) %>%
    arrange(StationID)

  # Write to the logs and make the code block and datatable display
  writelog(
    "\n### First Subsequent Step: RBind the CSI and LRM Dataframes",
    logfile = logfile,
    code = '
      combined1 <- rbind(chemdata_lrm, chemdata_csi) %>%
        arrange(StationID)
    ',
    data = combined1,
    verbose = verbose
  )

  # Make the download link
  create_download_link(data = combined1, logfile = logfile, filename = 'CSI_LRM_Combine_step1.csv', linktext = 'Download CSI_LRM_Combine_step1.csv', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  # ---- Second Subsequent Step: Group by StationID and take the average of CSI and LRM ----
  # --- If the average is not an integer, take the ceiling of it (round it up)

  # Here, we actually execute the code
  combined2 <- combined1 %>%
    group_by(StationID) %>%
    summarize(
      # we call ceiling because we are always wanting to round up
      # err on the side of determining that a site is more impacted, rather than not
      # na.rm is FALSE by default - so I'm not sure why it is specified as false in the Category Score one - It is the same exact thing
      Score = ceiling(mean(`Category Score`)),

      # na.rm = F, since we need both CSI and LRM to determine the Chemistry LOE score
      #   although if you can calculate one, you should be able to calculate the other as well
      # na.rm is FALSE by default - so I'm not sure why it is specified as false here - It is the same exact thing
      `Category Score` = ceiling(mean(`Category Score`, na.rm = F)) # Category Score is the same as the Score in this case
    ) %>%
    ungroup()

  # Write to the log
  writelog(
    "\n### Second Subsequent Step: Group by StationID and take the average of CSI and LRM",
    logfile = logfile,
    code = '
      combined2 <- combined1 %>%
        group_by(StationID) %>%
        summarize(
          # we call ceiling because we are always wanting to round up
          # err on the side of determining that a site is more impacted, rather than not
          # na.rm is FALSE by default - so I\'m not sure why it is specified as false in the Category Score one - It is the same exact thing
          Score = ceiling(mean(`Category Score`)),

          # na.rm = F, since we need both CSI and LRM to determine the Chemistry LOE score
          #   although if you can calculate one, you should be able to calculate the other as well
          # na.rm is FALSE by default - so I\'m not sure why it is specified as false here - It is the same exact thing
          `Category Score` = ceiling(mean(`Category Score`, na.rm = F)) # Category Score is the same as the Score in this case

        ) %>%
        ungroup()
    ',
    data = combined2,
    verbose = verbose
  )

  # Make the download link
  create_download_link(data = combined2, logfile = logfile, filename = 'CSI_LRM_Combine_step2.csv', linktext = 'Download CSI_LRM_Combine_step2.csv', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  # ---- Third Subsequent Step: Convert numeric score to the human readable category ----

  # Actually execute the code
  combined3 <- combined2 %>%
    mutate(
      Index = "Integrated SQO",
      Category = case_when(
        Score == 1 ~ "Minimal Exposure",
        Score == 2 ~ "Low Exposure",
        Score == 3 ~ "Moderate Exposure",
        Score == 4 ~ "High Exposure",
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      StationID, Index, Score, Category, `Category Score`
    )

  # Write to the log
  writelog(
    "\n### Third Subsequent Step: Convert numeric score to the human readable category",
    logfile = logfile,
    code = '
      combined3 <- combined2 %>%
        mutate(
          Index = "Integrated SQO",
          Category = case_when(
            Score == 1 ~ "Minimal Exposure",
            Score == 2 ~ "Low Exposure",
            Score == 3 ~ "Moderate Exposure",
            Score == 4 ~ "High Exposure",
            TRUE ~ NA_character_
          )
        ) %>%
        select(
          StationID, Index, Score, Category, `Category Score`
        )
    ',
    data = combined3,
    verbose = verbose
  )

  # Make the download link
  create_download_link(data = combined3, logfile = logfile, filename = 'CSI_LRM_Combine_step3.csv', linktext = 'Download CSI_LRM_Combine_step3.csv', verbose = verbose)

  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)

  # ---- Fourth Subsequent Step: Combine CSI, LRM, and Integrated Score dataframes, convert factors to characters, and order it by StationID ----
  combined.final <- combined3 %>%
    rbind(chemdata_lrm, chemdata_csi) %>%
    mutate_if(is.factor,as.character) %>%
    arrange(
      StationID, Index
    )

  # Write to the log
  writelog(
    "\n### Fourth Subsequent Step: Combine CSI, LRM, and Integrated Score dataframes, convert factors to characters, and order it by StationID",
    logfile = logfile,
    code = '
      combined.final <- combined3 %>%
        rbind(chemdata_lrm, chemdata_csi) %>%
        mutate_if(is.factor,as.character) %>%
        arrange(
          StationID, Index
        )
    ',
    data = combined.final,
    verbose = verbose
  )

  # Make the download link
  create_download_link(data = combined.final, logfile = logfile, filename = 'FinalChemSQOScores.csv', linktext = 'Download FinalChemSQOScores.csv', verbose = verbose)


  # Log the separation space between steps
  writelog("\n___\n___\n___" , logfile = logfile, verbose = verbose)


  writelog("\nEND Chem SQO Function\n", logfile = logfile, verbose = verbose)

  return(combined.final)

}


# ---- Utility function ----
#' Prepare raw chemistry data for analysis, per specifications in the
#'
#' This function will not only calculate CSI and LRM, but it
#' will also get the overall integrated SQO Category and Score
#' for the stations that are given to it.
#'
#' We make the assumption that Non-detects are marked with -88 in the result column. We also assume that
#' all constituents are expressed on a sediment dry-weight basis.
#' Specifically, all metals should be in mg/dry kg and all organic constituents should be in ug/dry kg.
#'
#' @usage chemdata_prep(chemdata)
#'
#' @param chemdata_prep.input a dataframe with the following columns:
#'
#'    \code{StationID},
#'
#'    \code{AnalyteName},
#'
#'    \code{Result},
#'
#'    \code{RL},
#'
#'    \code{MDL}
#'
#' @examples
#' data(chem_sampledata) # load sample data to your environment
#' chemdata_prep(chem_sampledata) # get scores and see output
#'
#' @export
chemdata_prep <- function(chemdata_prep.input, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'chemlog.Rmd' ), verbose = T){

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\n\n## Chemistry Preprocessing function (chemdata_prep)\n", logfile = logfile, verbose = verbose)


  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'chemdata_prep.input.RData'
  if (verbose) {
    save(chemdata_prep.input, file = file.path( dirname(logfile), rawinput.filename ))
  }

  writelog(
    "\n### Initial input to chemdata_prep:",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "')"),
    data = chemdata_prep.input %>% head(15),
    verbose = verbose
  )
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.csv', linktext = 'Download Initial Input to Chem Preprocessing Function', verbose = verbose)




  # Log the separation space between steps
  writelog("\n___\n___\n___\n___" , logfile = logfile, verbose = verbose)



  # ---- Step 1: lowercase column names ----
  writelog(
    "\n#### Step 1: Lowercase column names of input data:",
    logfile = logfile,
    verbose = verbose
  )
  names(chemdata_prep.input) <- names(chemdata_prep.input) %>% tolower()
  writelog(
    "",
    logfile = logfile,
    code = 'names(chemdata_prep.input) <- names(chemdata_prep.input) %>% tolower()',
    data = chemdata_prep.input %>% head(15),
    verbose = verbose,
    pageLength = 5
  )
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step1.csv', linktext = 'Download Preprocessing Data (Step 1)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 2: Rename fieldduplicate columns to "fieldrep" ----
  writelog("\n#### Step 2: Renaming fielddup colnames to 'fieldrep'", logfile = logfile, verbose = verbose)
  if ( 'fieldduplicate' %in% names(chemdata_prep.input) ) {
    chemdata_prep.input <- chemdata_prep.input %>% rename(fieldrep = fieldduplicate)
  }
  if ( 'fieldreplicate' %in% names(chemdata_prep.input) ) {
    chemdata_prep.input <- chemdata_prep.input %>% rename(fieldrep = fieldreplicate)
  }
  if ( 'fielddup' %in% names(chemdata_prep.input) ) {
    chemdata_prep.input <- chemdata_prep.input %>% rename(fieldrep = fielddup)
  }
  # Write code to the logs
  writelog(
    "",
    logfile = logfile,
    code = "
      if ( 'fieldduplicate' %in% names(chemdata_prep.input) ) {
        chemdata_prep.input <- chemdata_prep.input %>% rename(fieldrep = fieldduplicate)
      }
      if ( 'fieldreplicate' %in% names(chemdata_prep.input) ) {
        chemdata_prep.input <- chemdata_prep.input %>% rename(fieldrep = fieldreplicate)
      }
      if ( 'fielddup' %in% names(chemdata_prep.input) ) {
        chemdata_prep.input <- chemdata_prep.input %>% rename(fieldrep = fielddup)
      }
    ",
    data = chemdata_prep.input %>% head(15),
    verbose = verbose,
    pageLength = 5
  )
  # Serve it up for download
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step2.csv', linktext = 'Download Preprocessing Data (Step 2)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 3: Rename LabReplicate columns to "labrep" ----

  writelog("\n#### Step 3: Rename labreplicate to labrep", logfile = logfile, verbose = verbose)

  # Actually execute
  if ('labreplicate' %in% names(chemdata_prep.input)) {
    chemdata_prep.input <- chemdata_prep.input %>% rename(labrep = labreplicate)
  }

  # Write code to logs
  writelog(
    "",
    code = "
      if ('labreplicate' %in% names(chemdata_prep.input)) {
        chemdata_prep.input <- chemdata_prep.input %>% rename(labrep = labreplicate)
      }
    ",
    data = chemdata_prep.input %>% head(15),
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )

  # Serve it up for download
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step3.csv', linktext = 'Download Preprocessing Data (Step 3)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 4: Remove the Field/Lab Replicate #2 (Keep only fieldrep/labrep #1) ----
  writelog(
    "\n#### Step 4: Remove the Field/Lab Replicate #2 (Keep only fieldrep/labrep #1)",
    logfile = logfile,
    verbose = verbose
  )

  # Actually Execute the Step
  # Check for 'labrep' column
  if ('labrep' %in% names(chemdata_prep.input)) {
    chemdata_prep.input <- chemdata_prep.input %>%
      filter(as.numeric(labrep) == 1)
  } else {
    msg <- "Warning: Column 'labrep' was not provided - this may affect results if there are duplicate records for certain analytes"
    warning(msg)
    writelog(msg, logfile = logfile, verbose = verbose)
  }

  # Check for 'fieldrep' column
  if ('fieldrep' %in% names(chemdata_prep.input)) {
    chemdata_prep.input <- chemdata_prep.input %>%
      filter(as.numeric(fieldrep) == 1)
  } else {
    msg <- "Warning: Column 'fieldrep' was not provided - this may affect results if there are duplicate records for certain analytes"
    warning(msg)
    writelog(msg, logfile = logfile, verbose = verbose)
  }

  # Write to code portion of the log
  writelog(
    "",
    code = "
      # Check for 'labrep' column
      if ('labrep' %in% names(chemdata_prep.input)) {
        chemdata_prep.input <- chemdata_prep.input %>%
          filter(as.numeric(labrep) == 1)
      } else {
        msg <- \"Warning: Column 'labrep' was not provided - this may affect results if there are duplicate records for certain analytes\"
        warning(msg)
        writelog(msg, logfile = logfile, verbose = verbose)
      }

      # Check for 'fieldrep' column
      if ('fieldrep' %in% names(chemdata_prep.input)) {
        chemdata_prep.input <- chemdata_prep.input %>%
          filter(as.numeric(fieldrep) == 1)
      } else {
        msg <- \"Warning: Column 'fieldrep' was not provided - this may affect results if there are duplicate records for certain analytes\"
        warning(msg)
        writelog(msg, logfile = logfile, verbose = verbose)
      }
    ",
    data = chemdata_prep.input %>% head(15),
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )
  # Serve it up for download
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step4.csv', linktext = 'Download Preprocessing Data (Step 4)', verbose = verbose)




  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 5: Filter for only "Result" sampletypes (i.e. no matrix spikes, lab blanks, etc - those should definitely not be in here) ----
  writelog(
    '\n#### Step 5: Filter for only "Result" sampletypes',
    logfile = logfile,
    verbose = verbose
  )
  writelog(
    '\n##### (i.e. no matrix spikes, lab blanks, etc - those should definitely not be in here)',
    logfile = logfile,
    verbose = verbose
  )

  # SampleTypeCode should be "Result"
  # This is not explicitly from the SQO Manual, but rather just based on "Common sense" and guidance from Darrin Greenstein
  # We take labreplicate 1 to avoid bias.
  # Darrin said on 7/20/2022 "I'd just pick the first rep. That will randomly even out any bias."
  # Check for 'sampletypecode' column
  if ('sampletypecode' %in% names(chemdata_prep.input)) {
    chemdata_prep.input <- chemdata_prep.input %>%
      filter(sampletypecode == 'Result')
  } else {
    msg <- "Warning: Column 'sampletypecode' was not provided - this may affect results if there are duplicate records for certain analytes"
    warning(msg)
    writelog(msg, logfile = logfile, verbose = verbose)
  }

  # Write to code portion of the logs
  writelog(
    "",
    code = '
    if (\'sampletypecode\' %in% names(chemdata_prep.input)) {
      chemdata_prep.input <- chemdata_prep.input %>%
        filter(sampletypecode == \'Result\')
    } else {
      msg <- "Warning: Column \'sampletypecode\' was not provided - this may affect results if there are duplicate records for certain analytes"
      warning(msg)
    }
    ',
    data = chemdata_prep.input %>% head(15),
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )
  # Serve it up for download
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step5.csv', linktext = 'Download Preprocessing Data (Step 5)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 6: Check if the resulting chemdata_prep.input is empty after filtering ----
  writelog(
    "\n#### Step 6: Check if the resulting chemdata_prep.input is empty after filtering",
    logfile = logfile,
    verbose = verbose
  )
  if (nrow(chemdata_prep.input) == 0) {
    stop("In chemdata_prep - chemistry input data is empty after filtering sampletypecode == Result and labrep == 1 and fieldrep == 1")
  }

  # Write to code portion of the logs
  writelog(
    "",
    code = '
      if (nrow(chemdata_prep.input) == 0) {
        stop("In chemdata_prep - chemistry input data is empty after filtering sampletypecode == Result and labrep == 1 and fieldrep == 1")
      }
    ',
    data = chemdata_prep.input %>% head(15),
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )
  # Serve it up for download
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step6.csv', linktext = 'Download Preprocessing Data (Step 6)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 7: Drop duplicates ----
  writelog(
    "\n#### Step 7: Drop Duplicates",
    logfile = logfile,
    verbose = verbose
  )

  # define columns for which to check for duplicates
  dupcols <- intersect(names(chemdata_prep.input), c('stationid', 'analytename', 'sampletypecode', 'labrep', 'fieldrep'))
  writelog(
    "\n##### define columns for which to check for duplicates:",
    code = "dupcols <- intersect(names(chemdata_prep.input), c('stationid', 'analytename', 'sampletypecode', 'labrep', 'fieldrep'))",
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )


  # Isolate duplicates for display purposes
  chemdupes <- chemdata_prep.input[(chemdata_prep.input %>% select(all_of(dupcols)) %>% duplicated), ]
  writelog(
    "\n##### Duplicate records:",
    code = 'chemdupes <- chemdata_prep.input[(chemdata_prep.input %>% select(all_of(dupcols)) %>% duplicated), ]',
    data = chemdupes,
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )


  # Purge duplicates from actual input data
  chemdata_prep.input <- chemdata_prep.input[!(chemdata_prep.input %>% select(all_of(dupcols)) %>% duplicated), ]
  writelog(
    "\n##### Chemdata Prep Input Data with Duplicates Removed:",
    code = 'chemdata_prep.input <- chemdata_prep.input[!(chemdata_prep.input %>% select(all_of(dupcols)) %>% duplicated), ]',
    data = chemdata_prep.input %>% head(15),
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )
  # Serve them up for download
  create_download_link(data = chemdupes, logfile = logfile, filename = 'chemdata_prep.dupes.csv', linktext = 'Download Duplicated Records', verbose = verbose)
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step7.csv', linktext = 'Download Preprocessing Data (Step 7 - removed duplicates)', verbose = verbose)



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 8: Convert Result, RL, MDL to numeric fields ----
  writelog("\n#### Convert Result, RL, MDL to numeric fields", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata_prep.input <- chemdata_prep.input %>%
     mutate(
      result = as.numeric(result),
      rl = as.numeric(rl),
      mdl = as.numeric(mdl),
    ) %>%
    filter(!is.na(stationid))

  # Write to the code portion of the logs
  writelog(
    "",
    code = '
      chemdata_prep.input <- chemdata_prep.input %>%
         mutate(
          result = as.numeric(result),
          rl = as.numeric(rl),
          mdl = as.numeric(mdl),
        ) %>%
        filter(!is.na(stationid))
    ',
    data = chemdata_prep.input %>% head(15),
    logfile = logfile,
    verbose = verbose,
    pageLength = 5
  )
  # Serve it up for download
  create_download_link(data = chemdata_prep.input, logfile = logfile, filename = 'chemdata_prep.input.step8.csv', linktext = 'Download Preprocessing Data (Step 8)', verbose = verbose)




  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)




  # ------------------------------------------------------------------- Relevant Compounds ---------------------------------------------------------------------------


  writelog("\n#### Analytes used in Chem SQO Calc (Based on Table 3.1 in SQO Manual - pages 19 and 20):", logfile = logfile, verbose = verbose)

  # Analytes that are not grouped in any particular category
  single_analytes <- c('Cadmium','Copper','Lead','Mercury','Zinc',
                       'alpha-Chlordane','gamma-Chlordane','trans-Nonachlor',"4,4'-DDT")
  writelog(
    "\n##### Single Compounds",
    code =  '
      # Analytes that are not grouped in any particular category
      single_analytes <- c(\'Cadmium\',\'Copper\',\'Lead\',\'Mercury\',\'Zinc\',
        \'alpha-Chlordane\',\'gamma-Chlordane\',\'trans-Nonachlor\',"4,4\'-DDT")
    ',
    logfile = logfile,
    verbose = verbose
  )


  # High PAH
  # Table 3.1 in the SQO Manual (page 19)
  hpah <- c('Benz(a)anthracene', 'Benzo(a)pyrene', 'Benzo(e)pyrene',
            'Chrysene', 'Dibenz(a,h)anthracene', 'Fluoranthene', 'Perylene','Pyrene')
  writelog(
    "\n##### High PAHs",
    code =  "
      # High PAH
      # Table 3.1 in the SQO Manual (page 19)
      hpah <- c('Benz(a)anthracene', 'Benzo(a)pyrene', 'Benzo(e)pyrene',
                'Chrysene', 'Dibenz(a,h)anthracene', 'Fluoranthene', 'Perylene','Pyrene')
    ",
    logfile = logfile,
    verbose = verbose
  )

  # Low PAH
  # Table 3.1 in the SQO Manual (page 19)
  lpah <- c('1-Methylnaphthalene', '1-Methylphenanthrene', '2,6-Dimethylnaphthalene',
            '2-Methylnaphthalene', 'Acenaphthene', 'Anthracene',
            'Biphenyl', 'Fluorene', 'Naphthalene', 'Phenanthrene')
  writelog(
    "\n##### Low PAHs",
    code =  "
      # Low PAH
      # Table 3.1 in the SQO Manual (page 19)
      lpah <- c('1-Methylnaphthalene', '1-Methylphenanthrene', '2,6-Dimethylnaphthalene',
                '2-Methylnaphthalene', 'Acenaphthene', 'Anthracene',
                'Biphenyl', 'Fluorene', 'Naphthalene', 'Phenanthrene')
    ",
    logfile = logfile,
    verbose = verbose
  )

  # The PCB's that we care about
  # Table 3.1 in the SQO Manual (page 20)
  relevant_pcbs <- paste(
    'PCB', c('008', '018', '028', '044', '052', '066', '101', '105', '110', '118', '128', '138', '153', '180', '187', '195'), sep = '-'
  )
  writelog(
    "\n##### PCB Compounds Which are used in Chemistry SQO",
    code =  "
      # The PCB's that we care about
      # Table 3.1 in the SQO Manual (page 20)
      relevant_pcbs <- paste(
        'PCB', c('008', '018', '028', '044', '052', '066', '101', '105', '110', '118', '128', '138', '153', '180', '187', '195'), sep = '-'
      )
    ",
    logfile = logfile,
    verbose = verbose
  )

  # (For the sake of Logging)
  # Find the maximum length of the columns
  max_length <- max(length(single_analytes), length(hpah), length(lpah), length(relevant_pcbs))

  # Create the data frame
  analytes_df <- data.frame(
    `Single Compounds` = c(single_analytes, rep(NA, max_length - length(single_analytes)) ),
    `High PAHs` = c(hpah, rep(NA, max_length - length(hpah))),
    `Low PAHs` = c(lpah, rep(NA, max_length - length(lpah))),
    `PCBs` = c(relevant_pcbs, rep(NA, max_length - length(relevant_pcbs)))
  )
  writelog(
    "\n##### Table of All compounds used in Chemistry SQO calculation (Except DDDs, DDEs and 2,4'-DDT)",
    logfile = logfile,
    verbose = verbose
  )
  writelog(
    "\n##### These are handled separately since they are aggregated.",
    logfile = logfile,
    verbose = verbose
  )
  writelog(
    "\n##### Furthermore - the final preprocessed data needs 4,4'-DDT by itself as well as the total DDTs, so those must be taken care of separately",
    code = '
      # Find the maximum length of the columns
      max_length <- max(length(single_analytes), length(hpah), length(lpah), length(relevant_pcbs))

      # Create the data frame
      analytes_df <- data.frame(
        `Single Compounds` = c(single_analytes, rep(NA, max_length - length(single_analytes)) ),
        `High PAHs` = c(hpah, rep(NA, max_length - length(hpah))),
        `Low PAHs` = c(lpah, rep(NA, max_length - length(lpah))),
        `PCBs` = c(relevant_pcbs, rep(NA, max_length - length(relevant_pcbs)))
      )
    ',
    data = analytes_df,
    logfile = logfile,
    verbose = verbose,
    pageLength = 16
  )
  # ----------------------------------------------------------------------------------------------------------------------------------------------------------------- #



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 9: Filter for relevant compounds (Excluding DDTs which will be handled separately later) ----
  writelog(
    "\n#### Step 9: Filter for relevant compounds (Excluding DDTs which will be handled separately later)",
    logfile = logfile,
    verbose = verbose
  )

  # Actually execute the code
  chemdata.filtered.no.DDT <- chemdata_prep.input %>%
    # create a new column called compound. This is what we will group by,
    mutate(
      compound = case_when(
        analytename %in% single_analytes ~ analytename,
        analytename %in% relevant_pcbs ~ "PCBs_total",
        analytename %in% hpah ~ "HPAH",
        analytename %in% lpah ~ "LPAH",

        # THe Manual (3rd Edition - page 35) states that the Total DDT's DDE's and DDD's are the sum of 2,4' and 4,4' - because those are the only ones measured in Bight
        # We decided to leave the script this way
        # If for whatever odd reason, some weird compound that shouldnt have even gotten measured (2,2') made its way into the dataset - we may as well include it
        # June 3, 2024
        grepl("DDD",analytename) ~ "DDDs_total",
        grepl("DDE",analytename) ~ "DDEs_total",
        # The DDT total will be handled separately
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(compound)
    )

  # Write to code portion of the logs
  writelog(
    "",
    code = '
      chemdata.filtered.no.DDT <- chemdata_prep.input %>%
        # create a new column called compound. This is what we will group by,
        mutate(
          compound = case_when(
            analytename %in% single_analytes ~ analytename,
            analytename %in% relevant_pcbs ~ "PCBs_total",
            analytename %in% hpah ~ "HPAH",
            analytename %in% lpah ~ "LPAH",

            # The Manual (3rd Edition - page 35) states that the Total DDT\'s DDE\'s and DDD\'s are the sum of 2,4\' and 4,4\' - because those are the only ones measured in Bight
            # We decided to leave the script this way
            # If for whatever odd reason, some weird compound that shouldnt have even gotten measured (2,2\') made its way into the dataset - we may as well include it
            # June 3, 2024
            grepl("DDD",analytename) ~ "DDDs_total",
            grepl("DDE",analytename) ~ "DDEs_total",
            # The DDT total will be handled separately
            TRUE ~ NA_character_
          )
        ) %>%
        filter(
          !is.na(compound)
        )
    ',
    data = chemdata.filtered.no.DDT %>% head(15),
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(chemdata.filtered.no.DDT, logfile = logfile, filename = 'chemdata_prep.input.step9.no.ddt.csv', verbose = verbose, linktext = 'Chem Preprocessed Data Filtered (No DDT)')



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ----------------------------------------------------------------- Check Units of Data -----------------------------------------------------------------------
  writelog('\n#### Check the units of the input data', logfile = logfile, verbose = verbose)

  # v0.3.5 update
  # June 3, 2024 - check units
  # 'Cadmium','Copper','Lead','Mercury','Zinc' - should be in Parts Per Million (ug/g, ug/g dw, or ppm)
  # Everything else should be in Parts Per Billion (ng/g, ng/g dw, or ppb)

  # Define the analytes that should be in ppm
  ppm_analytes <- c('Cadmium', 'Copper', 'Lead', 'Mercury', 'Zinc')

  # Define the acceptable units for ppm and ppb
  acceptable_units_ppm <- c('ug/g', 'ug/g dw', 'ppm')
  acceptable_units_ppb <- c('ng/g', 'ng/g dw', 'ppb')

  # Check if the units column exists
  if ('units' %in% names(chemdata.filtered.no.DDT)) {
    # Check for rows that do not meet the criteria
    incorrect_rows <- chemdata.filtered.no.DDT[
      (chemdata.filtered.no.DDT$analytename %in% ppm_analytes & !(chemdata.filtered.no.DDT$units %in% acceptable_units_ppm)) |
        (!chemdata.filtered.no.DDT$analytename %in% ppm_analytes & !(chemdata.filtered.no.DDT$units %in% acceptable_units_ppb)), ]

    # If there are incorrect rows, stop and display a message
    if (nrow(incorrect_rows) > 0) {
      stop("There are rows with incorrect units. Please check the following rows:\n",
           paste(capture.output(print(incorrect_rows)), collapse = "\n"))
    } else {
      writelog("INFO: All units for chemistry input data are correct.",logfile = logfile,verbose = verbose)
    }
  } else {
    warning("The 'units' column is missing from the chemistry input data")
    writelog("WARNING: The 'units' column is missing from the chemistry input data",logfile = logfile,verbose = verbose)
  }

  # write code to the Log R Markdown
  writelog(
    "",
    code = '
      # June 3, 2024 - check units
      # \'Cadmium\',\'Copper\',\'Lead\',\'Mercury\',\'Zinc\' - should be in Parts Per Million (ug/g, ug/g dw, or ppm)
      # Everything else should be in Parts Per Billion (ng/g, ng/g dw, or ppb)

      # Define the analytes that should be in ppm
      ppm_analytes <- c(\'Cadmium\', \'Copper\', \'Lead\', \'Mercury\', \'Zinc\')

      # Define the acceptable units for ppm and ppb
      acceptable_units_ppm <- c(\'ug/g\', \'ug/g dw\', \'ppm\')
      acceptable_units_ppb <- c(\'ng/g\', \'ng/g dw\', \'ppb\')

      # Check if the units column exists
      if (\'units\' %in% names(chemdata.filtered.no.DDT)) {
        # Check for rows that do not meet the criteria
        incorrect_rows <- chemdata.filtered.no.DDT[
          (chemdata.filtered.no.DDT$analytename %in% ppm_analytes & !(chemdata.filtered.no.DDT$units %in% acceptable_units_ppm)) |
            (!chemdata.filtered.no.DDT$analytename %in% ppm_analytes & !(chemdata.filtered.no.DDT$units %in% acceptable_units_ppb)), ]

        # If there are incorrect rows, stop and display a message
        if (nrow(incorrect_rows) > 0) {
          stop("There are rows with incorrect units. Please check the following rows:\n",
               paste(capture.output(print(incorrect_rows)), collapse = "\n"))
        } else {
          print("INFO: All units for chemistry input data are correct.")
        }
      } else {
        warning("The \'units\' column is missing from the chemistry input data")
      }

    ',
    logfile = logfile,
    verbose = verbose
  )

  # ----------------------------------------------------------------------------------------------------------------------------------------------------- #



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 10: Deal with missing values/Non Detects ----
  writelog("\n#### Step 10: Deal with missing values/Non Detects",logfile = logfile,verbose = verbose)
  writelog("\n##### Dealing with missing result values (non detects)",logfile = logfile,verbose = verbose)
  writelog("\n##### NA or negative result values are treated as missing values (covers -88, -99 or actual null values)",logfile = logfile,verbose = verbose)
  chemdata.step10 <- chemdata.filtered.no.DDT %>%
    mutate(
      result = case_when(
        # This is for dealing with Non detects (Page 37 of SQO Manual, Paragraph titled Data Preparation)
        # Note that if calculations using non-detected(ND) analytes are necessary, an estimated value must be used.
        # One estimation approach is to use 50% of the MDL for any samples with ND results for that analyte;
        # however, the previous section should be consulted for addressing ND values within summed groups of constituents.
        # NA or negative result values are treated as missing values (covers -88, -99 or actual null values)
        (coalesce(result, -88) < 0) & (compound %in% single_analytes) ~ as.numeric(1/2*mdl),

        # For the summed group of constituents, we get the directions of how to deal with them in page 36 of the SQO Manual
        # First paragraph below table 3.4
        # "Compounds qualified as non-detected are treated as having a concentration of zero for the purpose of summing"
        # NA or negative result values are treated as missing values (covers -88, -99 or actual null values)
        (coalesce(result, -88) < 0) & !(compound %in% single_analytes) ~ 0,
        TRUE ~ result
      )
    )

  writelog(
    "",
    code = '
      chemdata.step10 <- chemdata.filtered.no.DDT %>%
        mutate(
          result = case_when(
            # This is for dealing with Non detects (Page 37 of SQO Manual, Paragraph titled Data Preparation)
            # Note that if calculations using non-detected(ND) analytes are necessary, an estimated value must be used.
            # One estimation approach is to use 50% of the MDL for any samples with ND results for that analyte;
            # however, the previous section should be consulted for addressing ND values within summed groups of constituents.
            # NA or negative result values are treated as missing values (covers -88, -99 or actual null values)
            (coalesce(result, -88) < 0) & (compound %in% single_analytes) ~ as.numeric(1/2*mdl),

            # For the summed group of constituents, we get the directions of how to deal with them in page 36 of the SQO Manual
            # First paragraph below table 3.4
            # "Compounds qualified as non-detected are treated as having a concentration of zero for the purpose of summing"
            # NA or negative result values are treated as missing values (covers -88, -99 or actual null values)
            (coalesce(result, -88) < 0) & !(compound %in% single_analytes) ~ 0,
            TRUE ~ result
          )
        )
    ',
    data = chemdata.step10 %>% head(15),
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(chemdata.step10, logfile = logfile, filename = 'chemdata_prep.input.step10.no.ddt.csv', verbose = verbose, linktext = 'Chem Preprocessed Data Step 10 (Dealing with non detects)')



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 11: Multiply PCB Values by 1.72 ----
  writelog("\n#### Step 11: Multiply PCB Values by 1.72", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata.step11 <- chemdata.step10 %>%
    mutate(
      result = if_else(
        # CASQO manual page 36 (3rd edition June 2021), below table 3.4 - PCB result value gets multiplied by 1.72
        # Modified - check where compound is PCBs_total rather than analytename (2024-05-13 in version 0.3.0 update) - Robert Butler
        compound == "PCBs_total", 1.72 * result, result
      )
    )

  # Write code to R Markdown
  writelog(
    "",
    code = '
      chemdata.step11 <- chemdata.step10 %>%
        mutate(
          result = if_else(
            # CASQO manual page 36 (3rd edition June 2021), below table 3.4 - PCB result value gets multiplied by 1.72
            # Modified - check where compound is PCBs_total rather than analytename (2024-05-13 in version 0.3.0 update) - Robert Butler
            compound == "PCBs_total", 1.72 * result, result
          )
        )
    ',
    data = chemdata.step11 %>% head(15),
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(chemdata.step11, logfile = logfile, filename = 'chemdata_prep.input.step11.pcb.mult.no.ddt.csv', verbose = verbose, linktext = 'Chem Preprocessed Data Step 11 (Multiply PCBs by 1.72)')



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 12: If all are non detects in the group, assign it the max of the RLs ----
  writelog("\n#### Step 12: If all are non detects in the group, assign it the max of the RLs", logfile = logfile, verbose = verbose)

  # Actually execute the code
  chemdata.step12 <- chemdata.step11 %>%
    group_by(
      stationid, compound
    ) %>%
    summarize(
      # if the sum of the results is zero, assign it the max of the RL's
      # Page 36 of SQO Plan, first paragraph below table 3.4
      # "If all components of a sum are non-detected, then the highest reporting limit of any one compound in the group should be used to represent the sum value."
      result = if_else(
        sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
      )
    ) %>%
    ungroup()

  # Write the code to the logs
  writelog(
    "\n##### Group by stationid and compound and sum the result values - if all are missing take the Max RL (Page 36 from Technical Manual)",
    code = '
      chemdata.step12 <- chemdata.step11 %>%
        group_by(
          stationid, compound
        ) %>%
        summarize(
          # if the sum of the results is zero, assign it the max of the RLs
          # Page 36 of SQO Plan, first paragraph below table 3.4
          # "If all components of a sum are non-detected, then the highest reporting limit of any one compound in the group should be used to represent the sum value."
          result = if_else(
            sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
          )
        ) %>%
    ungroup()
    ',
    data = chemdata.step12 %>% head(15),
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(chemdata.step12, logfile = logfile, filename = 'chemdata_prep.input.step12.no.ddt.csv', verbose = verbose, linktext = 'Chem Preprocessed Data Step 12 (Handle case where all in analytegroup are Non detects)')



  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)



  # ---- Step 13: Handle DDTs ----

  writelog("\n#### Step 13:Handle DDTs (according to page 36 of SQO technical manual)", logfile = logfile, verbose = verbose)
  writelog("\n- 13a. Replace missing values with zero", logfile = logfile, verbose = verbose)
  writelog("\n- 13b. Sum them - if all are non-detects then take the max RL", logfile = logfile, verbose = verbose)

  # Code
  ddts_total.13a <- chemdata_prep.input %>%
    filter(grepl("DDT",analytename)) %>%
    mutate(
      compound = "DDTs_total",
      result = if_else(
        # For the summed group of constituents, we get the directions of how to deal with them in page 30 of the SQO Manual
        # First paragraph below table 3.4
        # "Compounds qualified as non-detected are treated as having a concentration of zero for the purpose of summing"
        coalesce(result, -88) < 0, 0, result
      )
    )
  writelog(
    "\n##### 13a: Replace non detects with zero before grouping and summing",
    code = '
      ddts_total.13a <- chemdata_prep.input %>%
        filter(grepl("DDT",analytename)) %>%
        mutate(
          compound = "DDTs_total",
          result = if_else(
            # For the summed group of constituents, we get the directions of how to deal with them in page 30 of the SQO Manual
            # First paragraph below table 3.4
            # "Compounds qualified as non-detected are treated as having a concentration of zero for the purpose of summing"
            coalesce(result, -88) < 0, 0, result
          )
        )
    ',
    data = ddts_total.13a %>% head(15),
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(ddts_total.13a, logfile = logfile, filename = 'chemdata_prep.input.ddt.13a.csv', verbose = verbose, linktext = 'Chem Preprocessed Data (Filtered to DDTs - before grouping and summing)')



  # Code
  ddts_total <- ddts_total.13a %>%
    group_by(stationid,compound) %>%
    summarize(
      result = if_else(
        # if the sum of the results is zero, assign it the max of the RLs
        # Page 30 of SQO Plan, first paragraph below table 3.4
        # "If all components of a sum are non-detected, then the highest reporting limit of any one compound in the group should be used to represent the sum value."
        sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
      )
    ) %>% ungroup()

  writelog(
    "\n##### Group by stationid and compound (only one compound DDTs_total) and take the sum",
    code = '
      ddts_total <- ddts_total.13a %>%
        group_by(stationid,compound) %>%
        summarize(
          result = if_else(
            # if the sum of the results is zero, assign it the max of the RLs
            # Page 30 of SQO Plan, first paragraph below table 3.4
            # "If all components of a sum are non-detected, then the highest reporting limit of any one compound in the group should be used to represent the sum value."
            sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
          )
        ) %>%
        ungroup()
    ',
    data = ddts_total %>% head(15),
    logfile = logfile,
    verbose = verbose
  )
  # Serve it up for download
  create_download_link(ddts_total, logfile = logfile, filename = 'chemdata_prep.input.ddt.csv', verbose = verbose, linktext = 'Chem Preprocessed Data (Filtered to DDTs)')


  # Log the separation space between steps
  writelog("\n___\n___\n___\n___  " , logfile = logfile, verbose = verbose)




  # ---- Step 14: Combine data ----
  writelog("\n#### Step 14: Combine all analytes", logfile = logfile, verbose = verbose)

  # Code
  # Step 13 dealt only with DDTs - essentially ddts_total is the step 13 dataframe

  # write to the logs
  writelog("\n##### Combine data with DDTs", logfile = logfile, verbose = verbose)

  # Actually perform the code
  chemdata.preprocessed.final <- chemdata.step12 %>%
    bind_rows(ddts_total) %>%
    arrange(stationid, compound) %>%
    rename(
      StationID = stationid,
      AnalyteName = compound,
      Result = result
    )

  # Write code to the R Markdown file
  writelog(
    "Final Chemdata Prep Output",
    code = '
      chemdata.preprocessed.final <- chemdata.step12 %>%
        bind_rows(ddts_total) %>%
        arrange(stationid, compound) %>%
        rename(
          StationID = stationid,
          AnalyteName = compound,
          Result = result
        )
    ',
    data = chemdata.preprocessed.final %>% head(15),
    logfile = logfile,
    verbose = verbose
  )


  writelog("\n#### END Function: chemdata_prep\n", logfile = logfile, verbose = verbose)


  return(chemdata.preprocessed.final)

}



