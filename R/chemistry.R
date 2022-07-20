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
LRM <- function(chemdata) {
  "lrm_table"

  # get it in a usable format
  chemdata <- chemdata_prep(chemdata)


  # Take the Log10 of the chemistry concentration.
  chemdata <- chemdata %>%
    mutate(
      logResult = log10(Result)
    )

  # Combine LRM values with data based on the compound. Exclude compounds not in lrm.

  # Page 32 of Technical Manual
  chemdata_lrm <- lrm_table %>% left_join(chemdata, by = 'AnalyteName') %>%
    mutate(
      p = (exp(B0 + B1 * logResult) / (1 + exp(B0 + B1 * logResult))) %>% round(2)
    ) %>%
    # Page 33 of Technical Manual
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
CSI <- function(chemdata) {
  "csi_weight"

  # get it in a usable format
  chemdata <- chemdata_prep(chemdata)

  # Combine CSI Weight values with data based on the compound. Exclued compounds not in CSI calculation.
  chemdata_csi <- csi_weight %>% left_join(chemdata, by = "AnalyteName")

  chemdata_csi <- csi_weight %>%
    left_join(chemdata, by = "AnalyteName") %>%
    mutate(
      weight = case_when(
        # we needed this because the weights are only summed for mean CSI if the Result value exists
        # So we have to make it NA_real_ where the Result is NA
        !is.na(Result) ~ weight,
        TRUE ~ NA_real_
      ),
      # Page 34 of Technical Manual
      exposure_score = case_when(
        Result <= BDC1 ~ 1,
        (Result > BDC1) & (Result <= BDC2) ~ 2,
        (Result > BDC2) & (Result <= BDC3) ~ 3,
        Result > BDC3 ~ 4,
        TRUE ~ NA_real_
      ),
      # Page 35 of Technical Manual
      weighted_score = exposure_score * weight # manual calls this the weighted score (page 35)
    ) %>%
    group_by(
      # We need to also group by SampleDate.
      # Current SampleData does not have a date in it
      StationID #,SampleDate ?
    ) %>%
    summarize(
      # Page 35 of Technical Manual (Chart)
      weighted_score_sum = sum(weighted_score, na.rm = T),
      weight = sum(weight, na.rm = T),
      Score = round(weighted_score_sum / weight, 2)
    ) %>%
    ungroup() %>%
    select(-c(weighted_score_sum, weight)) %>%
    mutate(
      # Page 36 of Technical Manual
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
chem.sqo <- function(chemdata) {

  chemdata_lrm <- LRM(chemdata)
  chemdata_csi <- CSI(chemdata)

  # We should probably put some checks on the input data here
  # if it doesn't meet requirements, call stop function

  rbind(chemdata_lrm, chemdata_csi) %>%
    arrange(StationID) %>%
    group_by(StationID) %>%
    summarize(
      # we call ceiling because we are always wanting to round up
      # err on the side of determining that a site is more impacted, rather than not
      Score = ceiling(mean(`Category Score`)),
      `Category Score` = ceiling(mean(`Category Score`)) # Category Score is the same as the Score in this case
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
chemdata_prep <- function(chem){
  # Here chemdata consists of data in the same format as our database, with the columns
  # stationid, analytename, result, rl, mdl

  # lowercase column names
  names(chem) <- names(chem) %>% tolower()

  if ('labreplicate' %in% names(chem)) {
    chem <- chem %>% rename(labrep = labreplicate)
  }

  # SampleTypeCode should be "Result" and LabReplicate should be 1
  # This is not explicitly from the SQO Manual, but rather just based on "Common sense" and guidance from Darrin Greenstein
  # We take labreplicate 1 to avoid bias.
  # Darrin said on 7/20/2022 "I'd just pick the first rep. That will randomly even out any bias."
  if (all( c('sampletypecode','labrep') %in% names(chem)  )) {
    chem <- chem %>%
      filter(
        sampletypecode == 'Result'
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




  # result, rl, mdl should be numeric fields
  chem <- chem %>%
     mutate(
      result = as.numeric(result),
      rl = as.numeric(rl),
      mdl = as.numeric(mdl),
    ) %>%
    filter(!is.na(stationid))

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
    'PCB', c(008, 018, 028, 044, 052, 066, 101, 105, 110, 118, 128, 138, 153, 180, 187, 195), sep = "-"
  )

  # The other group is the "DD's" DDT, DDE, DDD
  # Can be done with grepl

  # its hard to think of useful, good variable names. Cut me some slack.
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
        analytename == "PCBs_total", 1.72 * result, result
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


  chemdata <- chemdata %>%
    bind_rows(ddts_total) %>%
    arrange(stationid, compound) %>%
    rename(
      StationID = stationid,
      AnalyteName = compound,
      Result = result
    )

  return(chemdata)

}



