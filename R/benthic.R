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
#'    \code{Stratum} - The stratum under which the station falls (Bays, Estuaries, etc);
#'
#'    \code{Exclude} - Yes or No;
#'
#' @importFrom dplyr bind_rows
#'
#' @export

benthic.sqo <- function(benthic_data, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'benthiclog.Rmd' ), verbose = F, knitlog = F){

  # Initialize Logging
  logfile.type <- ifelse(tolower(tools::file_ext(logfile)) == 'rmd', 'RMarkdown', 'text')
  init.log(logfile, base.func.name = sys.call(), type = logfile.type, current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)


  writelog("\n# Benthic SQO Main Function\n", logfile = logfile, verbose = verbose)

  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'benthic.sqo.input.RData'
  if (verbose) {
    save(benthic_data, file = file.path( dirname(logfile), rawinput.filename ))
  }


  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n### Raw input to benthic.sqo:\n",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') ### This will load a dataframe called 'benthic_data' into your environment"),
    verbose = verbose
  )
  create_download_link(data = benthic_data, logfile = logfile, filename = 'BenthicSQO-RawInput.csv', linktext = 'Download Raw Input to Benthic SQO Function', verbose = verbose)



  writelog('\n  \n## Calling M-AMBI within Benthic SQO function', logfile = logfile, verbose = verbose)
  mambi.score <- MAMBI(benthic_data, logfile = logfile , verbose = verbose) %>%
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
  # Create code block and download link to RBI output
  writelog(
    "\n## MAMBI output - Here is its output along with a code block (for R Studio users)\n",
    logfile = logfile,
    code = '
      mambi.score <- MAMBI(benthic_data, verbose = FALSE) %>%
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
    ',
    data = mambi.score,
    verbose = verbose
  )
  create_download_link(data = mambi.score, logfile = logfile, filename = 'MAMBI-output.csv', linktext = 'MAMBI output', verbose = verbose)

  writelog('\n#### ***As of now, since M-AMBI is not adopted into the official Benthic Integrated SQO score, details on M-AMBI will be excluded from the audit logs**', logfile = logfile, verbose = verbose)


  writelog('\n## Calling IBI within Benthic SQO function', logfile = logfile, verbose = verbose)
  ibi.scores <- IBI(benthic_data, logfile = logfile, verbose = verbose)
  # Create code block and download link to IBI output
  writelog(
    "\n## IBI function is finished executing - Here is its output along with a code block (for R Studio users):",
    logfile = logfile,
    code = "ibi.scores <- IBI(benthic_data, verbose = FALSE)",
    data = ibi.scores,
    verbose = verbose
  )
  create_download_link(data = ibi.scores, logfile = logfile, filename = 'IBI-output.csv', linktext = 'Download IBI function output', verbose = verbose)




  writelog('\n## Calling RBI within Benthic SQO function', logfile = logfile, verbose = verbose)
  rbi.scores <- RBI(benthic_data, logfile = logfile, verbose = verbose)
  # Create code block and download link to RBI output
  writelog(
    "\n## RBI function is finished executing - Here is its output along with a code block (for R Studio users):",
    logfile = logfile,
    code = "rbi.scores <- RBI(benthic_data, verbose = FALSE)",
    data = rbi.scores,
    verbose = verbose
  )
  create_download_link(data = rbi.scores, logfile = logfile, filename = 'RBI-output.csv', linktext = 'Download RBI output', verbose = verbose)



  writelog('\n## Calling BRI within Benthic SQO function', logfile = logfile, verbose = verbose)
  bri.scores <- BRI(benthic_data, logfile = logfile, verbose = verbose)
  # Create code block and download link to BRI output
  writelog(
    "\n## BRI function is finished executing - Here is its output along with a code block (for R Studio users):",
    logfile = logfile,
    code = "bri.scores <- BRI(benthic_data, verbose = FALSE)",
    data = bri.scores,
    verbose = verbose
  )
  create_download_link(data = bri.scores, logfile = logfile, filename = 'BRI-output.csv', linktext = 'Download BRI function output', verbose = verbose)

  writelog('\n## Calling RIVPACS within Benthic SQO function\n', logfile = logfile, verbose = verbose)
  rivpacs.score <- RIVPACS(benthic_data, logfile = logfile, verbose = verbose) #only SoCal (no SFBay)
  # Create code block and download link to RIVPACS output
  writelog(
    "\n## RIVPACS function is finished executing - Here is its output along with a code block (for R Studio users):",
    logfile = logfile,
    code = "rivpacs.score <- RIVPACS(benthic_data, verbose = FALSE)",
    data = rivpacs.score,
    verbose = verbose
  )
  create_download_link(data = rivpacs.score, logfile = logfile, filename = 'RIVPACS-output.csv', linktext = 'Download RIVPACS function output', verbose = verbose)

  # Integrated Scores
  # CASQO Technical Manual page 73 -
  #     Simply take the ceiling of the median of BRI, RBI, IBI and RIVPACS
  writelog('\n## Calculate Benthic Integrated scores\n  These will be visible in the final scores dataframe (benthic.sqo-final.csv)', logfile = logfile, verbose = verbose)
  writelog('Line that calculates integrated benthic score: `Category Score` = ceiling(median(`Category Score`, na.rm = T))', logfile = logfile, verbose = verbose)
  writelog('In other words - group the data by StationID, Replicate, SampleDate, Stratum and then take the ceiling of the mean - excluding missing values', logfile = logfile, verbose = verbose)

  # Stack the different scores dataframes on top of each other
  integrated.score.step1 <- bind_rows(
      rbi.scores, ibi.scores, bri.scores, rivpacs.score
    )
  # Create code block and download link
  writelog(
    "\n### Step 1 combining data for the Integrated Benthic Score - stack the dataframes",
    logfile = logfile,
    code = "integrated.score.step1 <- bind_rows(
      rbi.scores, ibi.scores, bri.scores, rivpacs.score
    )",
    data = integrated.score.step1,
    verbose = verbose
  )
  create_download_link(data = integrated.score.step1, logfile = logfile, filename = 'IntegratedBenthicStep1.csv', linktext = 'Download Integrated Benthic Step 1', verbose = verbose)


  # Filter to rep 1 and select only needed columns
  integrated.score.step2 <- integrated.score.step1 %>%
    # David says take only where replicate = 1, although other scientists have different opinions
    filter(Replicate == 1) %>%
    select(
      StationID, Replicate, SampleDate, Stratum, Index, `Category Score`
    )
  # Create code block and download link
  writelog(
    "\n### Step 2 combining data for the Integrated Benthic Score - filter to replicate 1 and select necessary columns in order",
    logfile = logfile,
    code = "
      # Filter to rep 1 and select only needed columns
      integrated.score.step2 <- integrated.score.step1 %>%
        filter(Replicate == 1) %>%
        select(
          StationID, Replicate, SampleDate, Stratum, Index, `Category Score`
        )
    ",
    data = integrated.score.step2,
    verbose = verbose
  )
  create_download_link(data = integrated.score.step2, logfile = logfile, filename = 'IntegratedBenthicStep2.csv', linktext = 'Download Integrated Benthic Step 2', verbose = verbose)


  # Group and get the ceiling of the median of all category scores for each station/rep/sampledate/stratum combo
  integrated.score.step3 <- integrated.score.step2 %>%
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
    ungroup()
  # Create code block and download link
  writelog(
    "\n### Step 3 Group and get the ceiling of the median of all category scores for each station/rep/sampledate/stratum combo",
    logfile = logfile,
    code = "
      # Group and get the ceiling of the median of all category scores for each station/rep/sampledate/stratum combo
      integrated.score.step3 <- integrated.score.step2 %>%
        group_by(
          StationID, Replicate, SampleDate, Stratum
        ) %>%
        summarize(
          `Category Score` = ceiling(median(`Category Score`, na.rm = T))
        ) %>%
        ungroup()
    ",
    data = integrated.score.step3,
    verbose = verbose
  )
  create_download_link(data = integrated.score.step3, logfile = logfile, filename = 'IntegratedBenthicStep3.csv', linktext = 'Download Integrated Benthic Step 2', verbose = verbose)


  # Last integrated score step
  integrated.score <- integrated.score.step3 %>%
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
  # Create code block and download link
  writelog(
    "\n### Final benthic integrated scores df",
    logfile = logfile,
    code = '
      # Last integrated score step
      integrated.score <- integrated.score.step3 %>%
        mutate(
          Index = "Integrated SQO",
          Category = case_when(
            `Category Score` == 1 ~ "Reference",
            `Category Score` == 2 ~ "Low Disturbance",
            `Category Score` == 3 ~ "Moderate Disturbance",
            `Category Score` == 4 ~ "High Disturbance",
            TRUE ~ NA_character_
          ),
          Score = `Category Score`
        )
    ',
    data = integrated.score,
    verbose = verbose
  )
  create_download_link(data = integrated.score, logfile = logfile, filename = 'IntegratedBenthic-Final.csv', linktext = 'Download Integrated Benthic Final Step', verbose = verbose)

  # Final dataframe
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

  # Create code block and download link
  writelog(
    "\n### Final benthic scores df (all scores)",
    logfile = logfile,
    code = '
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
    ',
    data = integrated.score,
    verbose = verbose
  )
  create_download_link(data = integrated.score, logfile = logfile, filename = 'AllBenthicScores-Final.csv', linktext = 'Download Final Benthic Scores df', verbose = verbose)


  writelog('\n# END: Benthic SQO function.\n', logfile = logfile, verbose = verbose)

  if (verbose && knitlog) {
    if ( tolower(tools::file_ext(logfile)) =='rmd' ) {

      html_file <- sub("\\.Rmd$", ".html", logfile, ignore.case = TRUE)

      print(paste0("Rendering ", logfile, " to ", html_file))
      rmarkdown::render(
        input = logfile,
        output_file = html_file,
        output_format = "html_document",
        quiet = TRUE
      )
      print("Done")

    } else {
      fn_name <- as.character(sys.call()[[1]])
      warning(paste0("In '", fn_name, "': knitlog = TRUE but the logfile is not an R Markdown (.Rmd) file. Skipping knitting."))
    }
  }

  return(final.scores)

}

