# ------------ LRM --------------
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
LRM <- function(chemdata, preprocessed = F, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ), verbose = T)  {
  "lrm_table"

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\nBEGIN Chem LRM Function\n", logfile = logfile, verbose = verbose)

  writelog("*** Chem data (Initial input to LRM function) may be found in LRM_initial_input_data.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chemdata, logfile = file.path(dirname(logfile), 'LRM_initial_input_data-step0.csv'), filetype = 'csv', verbose = verbose)

  # get it in a usable format
  if (!preprocessed) {
    writelog("Input to CSI function did not come preprocessed.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog("\nPrepping/Preprocessing chem data...", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    chemdata <- chemdata_prep(chemdata, logfile = logfile, verbose = verbose)
    writelog("\nDone prepping chem data...", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    writelog("*** Chem data (Input to LRM function after preprocessing) may be found in LRM_preprocessed_input_data-step0.1.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(chemdata, logfile = file.path(dirname(logfile), 'LRM_preprocessed_input_data-step0.1.csv'), filetype = 'csv', verbose = verbose)
  }


  # Take the Log10 of the chemistry concentration.
  writelog("Take the Log10 of the chemistry concentration.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  chemdata <- chemdata %>%
    mutate(
      logResult = log10(Result)
    )
  writelog("*** Chem Data for LRM function AFTER taking Log base 10 of Result Column found in LRM_chemdata_log10-step1.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chemdata, logfile = file.path(dirname(logfile), 'LRM_chemdata_log10-step1.csv'), filetype = 'csv', verbose = verbose)

  # Combine LRM values with data based on the compound. Exclude compounds not in lrm.

  writelog("Manipulate the LRM chem data per Technial Manual page 37-40", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Join with the LRM table by AnalyteName (Table 3.5 on page 39 of Technical Manual)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog("*** LRM table may be found in LRM_Table.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(lrm_table, logfile = file.path(dirname(logfile), 'LRM_Table.csv'), filetype = 'csv', verbose = verbose)

  writelog("-- p = (exp(B0 + B1 * logResult) / (1 + exp(B0 + B1 * logResult))) ---> Round to 2 decimal places", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Group by StationID (this function and R package overall assumes you are feeding it data for one sediment sample per stationid)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- take max value of the variable p within each grouping (grouped by StationID)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Sum the weighted scores and the weights - then take the sum of the weighted scores divided by the sum of the weights. - This is the CSI Score", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Assign LRM Categories based on thresholds defined in table 3.6 of Technical Manual page 40", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- Below 0.33 is 'Minimal'", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- From 0.33 to 0.49 is 'Low' (upper and lower bounds are inclusive)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- From 0.49 to 0.66 is 'Moderate' (lower bound is exclusive, upper bound is inclusive)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- Above 0.66 is 'High'", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  # Page 37 of Technical Manual
  chemdata_lrm <- lrm_table %>% left_join(chemdata, by = 'AnalyteName') %>%
    mutate(
      # page 38 of Technical Manual
      p = (exp(B0 + B1 * logResult) / (1 + exp(B0 + B1 * logResult))) %>% round(2)
    ) %>%
    # Page 39 of Technical Manual
    group_by(StationID) %>%
    summarize(
      Score = if_else(
        all(is.na(p)), NA_real_, suppressWarnings(max(p, na.rm = T))
      )
    ) %>%
    ungroup() %>%
    mutate(
      # Page 34 of Technical Manual
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

  writelog("*** Final LRM scores dataframe in LRM_scores-final.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chemdata_lrm, logfile = file.path(dirname(logfile), 'LRM_scores-final.csv'), filetype = 'csv', verbose = verbose)

  writelog("\nEND Chem LRM Function\n", logfile = logfile, verbose = verbose)
  writelog("----------------------------------------------------------------------------------------------------\n", logfile = logfile, verbose = verbose)

  return(chemdata_lrm)

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
CSI <- function(chemdata, preprocessed = F, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ), verbose = T) {
  "csi_weight"

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\nBEGIN Chem CSI Function\n", logfile = logfile, verbose = verbose)

  writelog("*** Chem data (Initial input to CSI function) may be found in csi_initial_input_data-step0.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chemdata, logfile = file.path(dirname(logfile), 'CSI_initial_input_data-step0.csv'), filetype = 'csv', verbose = verbose)

  # get it in a usable format
  if (!preprocessed) {
    writelog("Input to CSI function did not come preprocessed.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog("\nPrepping/Preprocessing chem data...", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    chemdata <- chemdata_prep(chemdata, logfile = logfile, verbose = verbose)
    writelog("\nDone prepping chem data...", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

    writelog("*** Chem data (Input to CSI function after preprocessing) may be found in LRM_preprocessed_input_data-step0.1.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    writelog(chemdata, logfile = file.path(dirname(logfile), 'CSI_preprocessed_input_data-step0.1.csv'), filetype = 'csv', verbose = verbose)

  }

  # Combine CSI Weight values with data based on the compound. Exclude compounds not in CSI calculation.
  print("Combine CSI Weight values with data based on the compound. Exclude compounds not in CSI calculation.", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  print("*** CSI Weight values may be found in CSI_weights.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(csi_weight, logfile = file.path(dirname(logfile), 'CSI_weights.csv'), filetype = 'csv', verbose = verbose)

  chemdata_csi <- csi_weight %>% left_join(chemdata, by = "AnalyteName")

  writelog("*** Chem data (Input to CSI function) combined with CSI Weight values may be found in CSI_input_data_with_weights-step1.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chemdata_csi, logfile = file.path(dirname(logfile), 'CSI_input_data_with_weights-step1.csv'), filetype = 'csv', verbose = verbose)

  writelog("Manipulate the CSI chem data per Technial Manual page 40-42", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Make the CSI weights NA where the Result value is NA", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Assign exposure score according to appropriate category (see Table 3.7 on page 40 of the Technical Manual)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Assign weighted exposure score (exposure_score x weight)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Group by StationID (this function and R package overall assumes you are feeding it data for one sediment sample per stationid)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Sum the weighted scores and the weights - then take the sum of the weighted scores divided by the sum of the weights. - This is the CSI Score", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Assign CSI Categories based on thresholds defined in table 3.9 of Technical Manual page 42", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- Below 1.69 is 'Minimal'", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- From 1.69 to 2.33 is 'Low' (upper and lower bounds are inclusive)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- From 2.33 to 2.99 is 'Moderate' (lower bound is exclusive, upper bound is inclusive)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("---- Above 2.99 is 'High'", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  chemdata_csi <- chemdata_csi %>%
    mutate(
      weight = case_when(
        # we needed this because the weights are only summed for mean CSI if the Result value exists
        # So we have to make it NA_real_ where the Result is NA
        !is.na(Result) ~ weight,
        TRUE ~ NA_real_
      ),
      # Page 40 of Technical Manual
      exposure_score = case_when(
        Result <= BDC1 ~ 1,
        (Result > BDC1) & (Result <= BDC2) ~ 2,
        (Result > BDC2) & (Result <= BDC3) ~ 3,
        Result > BDC3 ~ 4,
        TRUE ~ NA_real_
      ),
      # Page 41 of Technical Manual
      weighted_score = exposure_score * weight # manual calls this the weighted score (page 41)
    ) %>%
    group_by(
      # We need to also group by SampleDate.
      # Current SampleData does not have a date in it
      StationID #,SampleDate ?
    ) %>%
    summarize(
      # Page 41 of Technical Manual (Chart)
      weighted_score_sum = sum(weighted_score, na.rm = T),
      weight = sum(weight, na.rm = T),
      Score = round(weighted_score_sum / weight, 2)
    ) %>%
    ungroup() %>%
    select(-c(weighted_score_sum, weight)) %>%
    mutate(
      # Page 42 of Technical Manual
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
    mutate(
      Index = "CSI"
    ) %>%
    select(
      StationID, Index, Score, Category, `Category Score`
    )

  writelog("*** Final CSI scores dataframe in CSI_scores-final.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chemdata_csi, logfile = file.path(dirname(logfile), 'CSI_scores-final.csv'), filetype = 'csv', verbose = verbose)

  writelog("\nEND Chem CSI Function\n", logfile = logfile, verbose = verbose)
  writelog("----------------------------------------------------------------------------------------------------\n", logfile = logfile, verbose = verbose)

  return(chemdata_csi)
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
chem.sqo <- function(chemdata, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T) {

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\nBEGIN Chem SQO Function\n", logfile = logfile, verbose = verbose)

  writelog("About to preprocess chemistry data", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  chemdata <- chemdata_prep(chemdata, logfile = logfile, verbose = verbose)


  writelog("About to call LRM function within chem.sqo", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  chemdata_lrm <- LRM(chemdata, preprocessed = T, logfile = logfile, verbose = verbose)

  writelog("About to call CSI function within chem.sqo", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  chemdata_csi <- CSI(chemdata, preprocessed = T, logfile = logfile, verbose = verbose)

  # We should probably put some checks on the input data here
  # if it doesn't meet requirements, call stop function

  writelog("RBind LRM and CSI dataframes", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("Assign the Score according to the mean of the CSI and LRM score", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("If CSI and LRM are one step apart (for example, the CSI is 4 and the LRM is 3), take the ceiling of the mean (which basically would equate to the max of the two)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("The column called 'Category Score' will be the same as the score", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("The column called 'Category' is the 'human readable' SQO Assessment (Minimal, Low, Moderate, High)", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("The Final output dataframe will have the overall Chem Assessment score, along with the individual CSI and LRM values", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  combined <- rbind(chemdata_lrm, chemdata_csi) %>%
    arrange(StationID) %>%
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
    ungroup() %>%
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
    ) %>%
    rbind(chemdata_lrm, chemdata_csi) %>%
    mutate_if(is.factor,as.character) %>%
    arrange(
      StationID, Index
    )

  writelog("*** Final Chem SQO Dataframe in CHEMISTRY-SQO-FINAL.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(combined, logfile = file.path(dirname(logfile), 'CHEMISTRY-SQO-FINAL.csv'), filetype = 'csv', verbose = verbose)

  writelog("\nEND Chem SQO Function\n", logfile = logfile, verbose = verbose)
  writelog("----------------------------------------------------------------------------------------------------\n", logfile = logfile, verbose = verbose)

  return(combined)

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
#' chemdata_prep(chem_sampledata) # get scores and see output
#'
#' @export
chemdata_prep <- function(chem, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T){

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog("\nBEGIN Function: chemdata_prep\n", logfile = logfile, verbose = verbose)

  # Here chemdata consists of data in the same format as our database, with the columns
  # stationid, analytename, result, rl, mdl

  writelog("*** Chemistry data as fed into the function chemdata_prep may be found in chemdata-preprocessing-step0.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step0.csv'), filetype = 'csv', verbose = verbose)




  # lowercase column names
  names(chem) <- names(chem) %>% tolower()
  writelog("Lowercase column names of input data", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(
    "*** Chemistry data after lowercase of column names may be found in the function chemdata_prep may be found in chemdata-preprocessing-step1.1.csv",
    logfile = logfile,
    verbose = verbose,
    prefix = hyphen.log.prefix
  )
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step1.1.csv'), filetype = 'csv', verbose = verbose)

  writelog("Chemistry data after renaming fielddup colnames to 'fieldrep' may be found in chemdata-preprocessing-step1.2.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  if ( 'fieldduplicate' %in% names(chem) ) {
    chem <- chem %>% rename(fieldrep = fieldduplicate)
  }
  if ( 'fieldreplicate' %in% names(chem) ) {
    chem <- chem %>% rename(fieldrep = fieldreplicate)
  }
  if ( 'fielddup' %in% names(chem) ) {
    chem <- chem %>% rename(fieldrep = fielddup)
  }
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step1.2.csv'), filetype = 'csv', verbose = verbose)

  if ('labreplicate' %in% names(chem)) {
    writelog("Rename labreplicate to labrep", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
    chem <- chem %>% rename(labrep = labreplicate)
  }
  writelog("*** Chemistry data after renaming labrep colnames to 'labrep' may be found in chemdata-preprocessing-step1.3.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step1.3.csv'), filetype = 'csv', verbose = verbose)



  writelog(
    "Remove 2nd field/lab replicates from chem data. - chem data before removal of 2nd lab/field reps may be found in chemdata-preprocessing-step2.csv",
    logfile = logfile,
    verbose = verbose,
    prefix = hyphen.log.prefix
  )
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step2.csv'), filetype = 'csv', verbose = verbose)
  # SampleTypeCode should be "Result" and LabReplicate should be 1
  # This is not explicitly from the SQO Manual, but rather just based on "Common sense" and guidance from Darrin Greenstein
  # We take labreplicate 1 to avoid bias.
  # Darrin said on 7/20/2022 "I'd just pick the first rep. That will randomly even out any bias."
  if (all( c('sampletypecode','labrep', 'fieldrep') %in% names(chem)  )) {
    chem <- chem %>%
      filter(
        sampletypecode == 'Result'
      ) %>%
      filter(
        as.numeric(labrep) == 1
      ) %>%
      filter(
        as.numeric(labrep) == 1
      )

    if (nrow(chem_sampledata) == 0) {
      stop("In chemdata_prep - chemistry input data is empty after filtering sampletypecode == Result and labrep == 1")
    }
  } else {
    warning("Columns sampletypecode and labrep were not provided - this may affect results if there are duplicate records for certain analytes")
  }
  writelog(
    "Dropping duplicates from chem - chem data AFTER removal of duplicates may be found in chemdata-preprocessing-step3.csv",
    logfile = logfile,
    verbose = verbose,
    prefix = hyphen.log.prefix
  )
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step3.csv'), filetype = 'csv', verbose = verbose)


  writelog("Converting Result, RL, MDL to numeric fields", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  # result, rl, mdl should be numeric fields
  chem <- chem %>%
     mutate(
      result = as.numeric(result),
      rl = as.numeric(rl),
      mdl = as.numeric(mdl),
    ) %>%
    filter(!is.na(stationid))

  writelog(
    "*** Chem data after conversion of Result, RL, MDL may be found in chemdata-preprocessing-step4.csv",
    logfile = logfile,
    verbose = verbose,
    prefix = hyphen.log.prefix
  )
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step4.csv'), filetype = 'csv', verbose = verbose)


  # Analytes that are not grouped in any particular category
  single_analytes <- c('Cadmium','Copper','Lead','Mercury','Zinc',
                       'alpha-Chlordane','gamma-Chlordane','trans-Nonachlor',"4,4'-DDT")


  # High PAH
  # Table 3.1 in the SQO Manual (page 19)
  hpah <- c('Benz(a)anthracene', 'Benzo(a)pyrene', 'Benzo(e)pyrene',
            'Chrysene', 'Dibenz(a,h)anthracene', 'Fluoranthene', 'Perylene','Pyrene')

  # Low PAH
  # Table 3.1 in the SQO Manual (page 19)
  lpah <- c('1-Methylnaphthalene', '1-Methylphenanthrene', '2,6-Dimethylnaphthalene',
            '2-Methylnaphthalene', 'Acenaphthene', 'Anthracene',
            'Biphenyl', 'Fluorene', 'Naphthalene', 'Phenanthrene')

  # The PCB's that we care about
  # Table 3.1 in the SQO Manual (page 20)
  relevant_pcbs <- paste(
    'PCB', c('008', '018', '028', '044', '052', '066', '101', '105', '110', '118', '128', '138', '153', '180', '187', '195'), sep = "-"
  )

  writelog("Analytes used in Chem SQO Calc (Based on Table 3.1 in SQO Manual - pages 19 and 20):", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(
    paste(
      c(single_analytes, hpah, lpah, relevant_pcbs),
      collapse = ', '
    ),
    logfile = logfile,
    verbose = verbose,
    prefix = hyphen.log.prefix
  )
  writelog("To be clear - the script is filtering the chemistry data for the above analytes", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog('\n', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  writelog("These analytes have no grouping and are handled individually", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste(single_analytes, collapse = ', '),logfile = logfile,verbose = verbose,prefix = hyphen.log.prefix)

  writelog("\nThese analytes are the High PAHs", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste(hpah, collapse = ', '),logfile = logfile,verbose = verbose,prefix = hyphen.log.prefix)
  writelog("\nThese analytes are the Low PAHs", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste(lpah, collapse = ', '),logfile = logfile,verbose = verbose,prefix = hyphen.log.prefix)
  writelog("\nThese analytes are the the PCBs used in SQO analysis", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(paste(relevant_pcbs, collapse = ', '),logfile = logfile,verbose = verbose,prefix = hyphen.log.prefix)
  writelog('\n', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  # The other group is the "DD's" DDT, DDE, DDD
  # Can be done with grepl


  # Version 0.3.0 update (2024-05-13) - remove duplicates but give a warning to the user - Robert Butler
  # check for duplicates and issue a warning:
  if (any(chem %>% duplicated())) {
    warning("Duplicates detected in chemistry data. They will be removed in the calculation. Please check your input data.")
    writelog(
      "Duplicates detected in chemistry data. They will be removed in the calculation. Please check your input data."
      ,
      logfile = logfile,
      verbose = verbose,
      prefix = paste(hyphen.log.prefix, 'WARNING: ')
    )
    chem <- chem %>% distinct()
    writelog(
      "*** Duplicates removed - chemdata-preprocessing-step4.1.csv"
      ,
      logfile = logfile,
      verbose = verbose,
      prefix = paste(hyphen.log.prefix, 'WARNING: ')
    )
    writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step4.1.csv'), filetype = 'csv', verbose = verbose)
  }


  writelog(
    "*** Before filtering, grouping, converting missing result values - chemdata-preprocessing-step5.csv"
    ,
    logfile = logfile,
    verbose = verbose,
    prefix = paste(hyphen.log.prefix, 'WARNING: ')
  )
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step5.csv'), filetype = 'csv', verbose = verbose)

  writelog("filtering chemistry data to necessary analytes - etc", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("from step 5 to 6 this is what is performed:", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("Create a new column called compound which represents the group of the analyte (High/Low PAH, Total PCB, DDD, DDE. (DDTs not handled in this step)", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  writelog("remove records where compound is NA", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  writelog("Handle missing result values according to guidance of pages 30 and 31 of Technical manual", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  writelog("-- (half MDL for the analytes that are not part of a group, otherwise set to zero in the grouped analytes which get summed)", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  writelog("Multiply PCBs by 1.72", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  writelog("Sum results of grouped constituents/analytes", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  writelog("If all values in the group are non detects - take the max RL", logfile = logfile, verbose = verbose, prefix = paste0('-',hyphen.log.prefix))
  # its hard to think of useful, good variable names
  chemdata <- chem %>%
    # create a new column called compound. This is what we will group by,
    mutate(
      compound = case_when(
        analytename %in% single_analytes ~ analytename,
        analytename %in% relevant_pcbs ~ "PCBs_total",
        analytename %in% hpah ~ "HPAH",
        analytename %in% lpah ~ "LPAH",
        grepl("DDD",analytename) ~ "DDDs_total",
        grepl("DDE",analytename) ~ "DDEs_total",
        # The DDT total will be handled separately
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(compound)
    ) %>%
    mutate(
      result = case_when(
        # This is for dealing with Non detects (Page 31 of SQO Manual, Paragraph titled Data Preparation)
        # Note that if calculations using non-detected(ND) analytes are necessary, an estimated value must be used.
        # One estimation approach is to use 50% of the MDL for any samples with ND results for that analyte;
        # however, the previous section should be consulted for addressing ND values within summed groups of constituents.
        result == -88 & compound %in% single_analytes ~ as.numeric(1/2*mdl),

        # For the summed group of constituents, we get the directions of how to deal with them in page 30 of the SQO Manual
        # First paragraph below table 3.4
        # "Compounds qualified as non-detected are treated as having a concentration of zero for the purpose of summing"
        result == -88 & !(compound %in% single_analytes) ~ 0,
        TRUE ~ result
      )
    ) %>%
    mutate(
      result = if_else(
        # CASQO manual page 30, below table 3.4 - PCB result value gets multiplied by 1.72
        # Modified - check where compound is PCBs_total rather than analytename (2024-05-13 in version 0.3.0 update) - Robert Butler
        compound == "PCBs_total", 1.72 * result, result
      )
    ) %>%
    group_by(
      stationid, compound
    ) %>%
    summarize(
      # if the sum of the results is zero, assign it the max of the RL's
      # Page 30 of SQO Plan, first paragraph below table 3.4
      # "If all components of a sum are non-detected, then the highest reporting limit of any one compound in the group should be used to represent the sum value."
      result = if_else(
        sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
      )
    ) %>%
    ungroup()

  writelog("*** chemdata after performing these steps may be found in chemdata-preprocessing-step6.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step6.csv'), filetype = 'csv', verbose = verbose)


  writelog("Handle DDTs (according to page 30 of SQO technical manual", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- replace missing values with zero", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("-- Sum them - if all are non-detects then take the max RL", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  ddts_total <- chem %>%
    filter(grepl("DDT",analytename)) %>%
    mutate(
      compound = "DDTs_total",
      result = if_else(
        # For the summed group of constituents, we get the directions of how to deal with them in page 30 of the SQO Manual
        # First paragraph below table 3.4
        # "Compounds qualified as non-detected are treated as having a concentration of zero for the purpose of summing"
        result == -88, 0, result
      )
    ) %>%
    group_by(stationid,compound) %>%
    summarize(
      result = if_else(
        # if the sum of the results is zero, assign it the max of the RL's
        # Page 30 of SQO Plan, first paragraph below table 3.4
        # "If all components of a sum are non-detected, then the highest reporting limit of any one compound in the group should be used to represent the sum value."
        sum(result, na.rm = T) != 0, sum(result, na.rm = T), max(rl)
      )
    ) %>% ungroup()

  writelog("*** chemdata after performing these steps may be found in chemdata-preprocessing-step7 (DDTs).csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessing-step7 (DDTs).csv'), filetype = 'csv', verbose = verbose)

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
  writelog("Combine DDT data with rest of chem data - perform necessary rounding", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog("Copper, Lead, Zinc, High/Low PAHs will be rounded to 1 decimal place; All others will be rounded to 2", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  chemdata <- chemdata %>%
    bind_rows(ddts_total) %>%
    arrange(stationid, compound) %>%
    rename(
      StationID = stationid,
      AnalyteName = compound,
      Result = result
    ) %>%
    mutate(
      Result = case_when(
        AnalyteName %in% c('Copper','Lead','Zinc','HPAH','LPAH') ~ round(Result, digits = 1),
        TRUE ~ round(Result, digits = 2)
      )
    )

  writelog("*** final preprocessed chemistry data in chemdata-preprocessed-final.csv", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(chem, file.path(dirname(logfile), 'chemdata-preprocessed-final.csv'), filetype = 'csv', verbose = verbose)

  writelog("\nEND Function: chemdata_prep\n", logfile = logfile, verbose = verbose)
  writelog("----------------------------------------------------------------------------------------------------\n", logfile = logfile, verbose = verbose)

  return(chemdata)

}



