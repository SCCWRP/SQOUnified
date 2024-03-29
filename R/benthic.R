# ---- benthic sqo function ----
#' Get Benthic Related Indices, and the Integrated Benthic SQO Score
#'
#' This function will calculate, RBI, IBI, BRI, RIVPACS, MAMBI,
#' as well as the Integrated Benthic SQO score / category.
#' The ultimate guide for these indices (EXCEPT MAMBI) is the CASQO Technical Manual
#' pages 64 to 94
#'
#' @usage
#' benthic.sqo(benthic_data)
#'
#' @examples
#' data(benthic_sampledata) # load the sample data
#' benthic.sqo(benthic_sampledata) # get scores and see output
#'
#' @param benthic_data a data frame. This data frame must contain the following
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
#'    \code{Stratum} - ;
#'
#'    \code{Exclude} - ;
#'
#' @importFrom dplyr bind_rows
#'
#' @export

benthic.sqo <- function(benthic_data){


  mambi.score <- MAMBI(benthic_data) %>%
    # rename(
    #   Stratum = Stratum
    # ) %>%
    mutate(
      Score = MAMBI_Score,
      Category = New_MAMBI_Condition
    ) %>%
    mutate(
      `Category Score` = case_when(
        Category == "Reference" ~ 1,
        Category == "Low Disturbance" ~ 2,
        Category == "Moderate Disturbance" ~ 3,
        Category == "High Disturbance" ~ 4,
        TRUE ~ NA_real_
      )
    )
  rbi.scores <- RBI(benthic_data)
  ibi.scores <- IBI(benthic_data)
  bri.scores <- BRI(benthic_data)
  rivpacs.score <- RIVPACS(benthic_data) #only SoCal (no SFBay)

  # Integrated Scores
  # CASQO Technical Manual page 87 -
  #     Simply take the ceiling of the median of BRI, RBI, IBI and RIVPACS
  integrated.score <- bind_rows(
      rbi.scores, ibi.scores, bri.scores, rivpacs.score
    ) %>%
    # David says take only where replicate = 1, although other scientists have different opinions
    filter(Replicate == 1) %>%
    select(
      StationID, Replicate, SampleDate, Stratum, Index, `Category Score`
    ) %>%
    group_by(
      StationID, Replicate, SampleDate, Stratum
    ) %>%
    summarize(
      # I asked David Gillett on August 1, 2022 if we should say benthic is unknown if one index is missing, this is his response:
      #   "We don't really say in the guidance document.
      #   There really isn't a reason that one of the 4 indices couldn't be calculated if there is a sample.
      #   My thought is that I wouldn't want to return an unknown if there are fewer than 4 indices, for whatever magic reason that may have occurred."
      # Based on this answer, we will include the keyword argument, "na.rm = T"
      # -Robert Butler, August 1, 2022

      `Category Score` = ceiling(median(`Category Score`, na.rm = T))
    ) %>%
    ungroup() %>%
    mutate(
      Index = 'Integrated SQO',
      Category = case_when(
        `Category Score` == 1 ~ "Reference",
        `Category Score` == 2 ~ "Low Disturbance",
        `Category Score` == 3 ~ "Moderate Disturbance",
        `Category Score` == 4 ~ "High Disturbance",
        TRUE ~ NA_character_
      ),
      Score = `Category Score`
    )

  # will add other scores to this data frame as they are computed
  final.scores <- bind_rows(
      mambi.score,
      rbi.scores,
      ibi.scores,
      bri.scores,
      rivpacs.score,
      integrated.score
    ) %>%
    select(
      StationID, Replicate, SampleDate, Stratum, Index, Score, Category, `Category Score`, Use_MAMBI
    ) %>%
    arrange(StationID, SampleDate, Replicate)

  return(final.scores)

}

