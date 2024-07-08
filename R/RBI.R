#' Compute the relative benthic index (RBI) score.
#'
#' @description
#'   The RBI is the weighted sum of: (a) four community metrics related to biodiversity (total number of taxa, number of crustacean taxa, abundance
#'   of crustacean individuals, and number of mollusc taxa), (b) abundances of three positive indicator taxa, and (c) the presence of two negative
#'   indicator species.
#'
#' @details
#'   The RBI is the weighted sum of: (a) four community metrics related to biodiversity (total number of taxa, number of crustacean taxa, abundance
#'   of crustacean individuals, and number of mollusc taxa), (b) abundances of three positive indicator taxa, and (c) the presence of two negative
#'   indicator species.
#'
#'   The data needed to calculate the RBI are:
#'   (1) Total number of taxa,
#'   (2) Number of mollusc taxa,
#'   (3) Number of crustacean individuals,
#'   (4) Number of individuals of \emph{Monocorophium insidiosum},
#'   (5) Number of individuals of \emph{Asthenothaerus diegensis},
#'   (6) Number of individuals of \emph{Goniada littorea},
#'   (7) Whether the data has the presence of \emph{Capitella capitata} complex, and
#'   (8) Whether the data has the presence of Oligochaeta.
#'
#'   To compute the RBI, the first step is to normalize the values for the benthic community metrics relative to maxima for the data used to develop
#'   the RBI for the Southern California Marine Bays habitat, to produce values relative to the maxima that are referred to as scaled values. The
#'   scaled value calculations use the following formulae:
#'
#'   Total Number of Taxa / 99
#'   Number of Mollusc Taxa / 28
#'   Number of Crustacean Taxa / 29
#'
#'
#'   The next step is to calculate the Taxa Richness Weighted Value (TWV) from the scaled values by the equation:
#'
#'   TWV = Scaled Total Number of Taxa + Scaled Number of Mollusc Taxa + Scaled Number of Crustacean Taxa + (0.25 * Scaled Abundance of Crustacea)
#'
#'   Next, the value for the two negative indicator taxa (NIT) is calculated. The two negative indicator taxa are \emph{Capitella capitata} complex
#'   and Oligochaeta. For each of these taxa that are present, in any abundance, the NIT is decreased by 0.1. Therefore, if neither were found the
#'   NIT = 0, if both are found the NIT = -0.2.
#'
#'   The next step is to calculate the value for the three positive indicator taxa (PIT). The positive indicator taxa are \emph{Monocorophium insidiosum,
#'   Asthenothaerus diegensis}, and \emph{Goniada littorea}. First, the PIT value is calculated for each species using the following equations:
#'
#'   \deqn{\frac{\sqrt[4]{Monocorophium~ insidiosum \textrm{abundance}}}{\sqrt[4]{473}}}
#'   \deqn{\frac{\sqrt[4]{Asthenothaerus~ diegensis \textrm{abundance}}}{\sqrt[4]{27}}}
#'   \deqn{\frac{\sqrt[4]{Goniada littorea~ \textrm{abundance}}}{\sqrt[4]{15}}}
#'
#'   The three species PIT values are then summed to calculate the PIT value for the sample. If none of the three species is present, then the sample
#'   PIT = 0.
#'
#'   The next step is to calculate the Raw RBI:
#'
#'   \deqn{\textrm{Raw RBI} = \textrm{TWV + NIT + } (2 \times \textrm{PIT})}
#'
#'   The final calculation is for the RBI score, normalizing the Raw RBI by the minimum and maximum Raw RBI values in the index development data:
#'
#'   \deqn{\textrm{RBI Score} = (\textrm{Raw RBI} - 0.03)/4.69}
#'
#'   The last step in the RBI process is to compare the RBI Score to a set of thresholds to determine the RBI category (Table 4).
#'
#'   <Insert Table 4>
#'
#'   For the function to run, the following packages NEED to be installed:  tidyverse, reshape2, vegan, and readxl.
#'   Additionally the EQR.R function must also be installed and is included with this code.
#'
#'   The output of the function will be a dataframe with StationID, Replicate, SampleDate, Latitude, Longitude,
#'   SalZone (The Salinity Zone assigned by M-AMBI), AMBI_Score, S (Species Richness), H (Species Diversity),
#'   Oligo_pct (Relative Abundance of Oligochaetes), MAMBI_Score, Orig_MAMBI_Condition, New_MAMBI_Condition,
#'   Use_MAMBI (Can M-AMBI be applied?), Use_AMBI (Can AMBI be applied?), and YesEG (% of Abundance with a EG value)
#'
#' @usage BRI(benthic_data)
#'
#' @param BenthiCData a data frame stored in the R environment. Note that this data frame MUST contain the following
#'                    information with these headings:
#'
#'                         \code{StationID} - an alpha-numeric identifier of the location;
#'
#'                         \code{Replicate} - a numeric identifying the replicate number of samples taken at the location;
#'
#'                         \code{SampleDate} - the date of sample collection;
#'
#'                         \code{Latitude} - latitude in decimal degrees;
#'
#'                         \code{Longitude} - longitude in decimal degrees. Make sure there is a negative sign for the Western coordinates;
#'
#'                         \code{Taxon} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#'                         \code{Abundance} - the number of each Species observed in a sample;
#'
#'                         \code{Salinity} - the salinity observed at the location in PSU, ideally at time of sampling;
#'
#'                         \code{Stratum} - ;
#'
#'                         \code{Exclude} - ;
#'
#' @return The output of the function will be a dataframe with
#'
#'    StationID,
#'
#'    Replicate,
#'
#'    SampleDate,
#'
#'    Latitude,
#'
#'    Longitude,
#'
#'   SalZone (The Salinity Zone assigned by M-AMBI),
#'
#'   AMBI_Score,
#'
#'   S (Species Richness),
#'
#'   H (Species Diversity),
#'
#'   Oligo_pct (Relative Abundance of Oligochaetes),
#'
#'   MAMBI_Score,
#'
#'   Orig_MAMBI_Condition,
#'
#'   New_MAMBI_Condition,
#'
#'   Use_MAMBI (Can M-AMBI be applied?),
#'
#'   Use_AMBI (Can AMBI be applied?),
#'
#'   YesEG (% of Abundance with a EG value)
#'
#' @examples
#'   RBI(benthic_data)
#'   RBI(BenthicData)
#'
#'
#' @import dplyr
#' @importFrom tidyr replace_na
#'
#' @export

# RBI ----
RBI <- function(BenthicData, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = F)
{

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  writelog('\n## BEGIN: RBI function.\n', logfile = logfile, verbose = verbose)


  # ---- Save the raw input to an RData file (for the sake of those who want the auditing logs) ----
  rawinput.filename <- 'benthic.rbi.input.RData'
  if (verbose) {
    save(BenthicData, file = file.path( dirname(logfile), rawinput.filename ))
  }


  # Display raw input data, create a download link for the knitted final RMarkdown output
  writelog(
    "\n### Raw input to RBI:",
    logfile = logfile,
    code = paste0("load('", rawinput.filename, "') ### This will load a dataframe called 'BenthicData' into your environment"),
    verbose = verbose
  )
  create_download_link(data = BenthicData, logfile = logfile, filename = 'RBI-RawInput.csv', linktext = 'Download Raw Input to RBI Function', verbose = verbose)




  # Create code block and download link to RBI input
  writelog(
    "\nInput to RBI function - RBI-step0.csv",
    logfile = logfile,
    data = BenthicData %>% head(25),
    verbose = verbose
  )
  create_download_link(data = BenthicData, logfile = logfile, filename = 'RBI-step0.csv', linktext = 'Download RBI initial input', verbose = verbose)

  # RBI Step 1
  # Prepare the given data frame so that we can compute the RBI score and categories
  rbi_data <- BenthicData %>%
    dplyr::filter(Exclude!="Yes") %>%
    dplyr::left_join(sqo.list.new, by = c("Taxon"="TaxonName")) %>%
    dplyr::mutate_if(is.numeric, list(~na_if(., -88))) %>%
    dplyr::select("StationID","SampleDate", "Replicate","Taxon","Abundance","Stratum", "Phylum", "Mollusc", "Crustacean") %>%
    #dplyr::rename(B13_Stratum = Stratum) %>%
    dplyr::mutate(n=if_else(Taxon=="NoOrganismsPresent", 0,1))

  writelog(
    "\nRBI Step 1 - Prepare the given data frame so that we can compute the RBI score and categories",
    logfile = logfile,
    code = '
      # Prepare the given data frame so that we can compute the RBI score and categories
      rbi_data <- BenthicData %>%
        dplyr::filter(Exclude!="Yes") %>%
        dplyr::left_join(sqo.list.new, by = c("Taxon"="TaxonName")) %>%
        dplyr::mutate_if(is.numeric, list(~na_if(., -88))) %>%
        dplyr::select("StationID","SampleDate", "Replicate","Taxon","Abundance","Stratum", "Phylum", "Mollusc", "Crustacean") %>%
        #dplyr::rename(B13_Stratum = Stratum) %>%
        dplyr::mutate(n=if_else(Taxon=="NoOrganismsPresent", 0,1))
    ',
    data = rbi_data %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_data, logfile = logfile, filename = 'RBI-step1.csv', linktext = 'Download RBI step 1', verbose = verbose)


  # to be used later
  rbi_data_tmp <- rbi_data %>%
    dplyr::group_by(Stratum, SampleDate, StationID, Replicate) %>%
    dplyr::summarise(NumOfTaxa = sum(n))

  writelog(
    "\nA temp dataframe that is used later in the function",
    logfile = logfile,
    code = '
      # to be used later
      rbi_data_tmp <- rbi_data %>%
        dplyr::group_by(Stratum, SampleDate, StationID, Replicate) %>%
        dplyr::summarise(NumOfTaxa = sum(n))
    ',
    data = rbi_data_tmp %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_data_tmp, logfile = logfile, filename = 'RBI-tmp-dataframe.csv', linktext = 'Download RBI tmp dataframe', verbose = verbose)

  # RBI Step 2
  # Filter to molluscs and group by Stratum, StationID, Sampledate Replicate - group and get the number of mollusc taxa
  rbi2 <- rbi_data %>%
    dplyr::filter(Mollusc=="Mollusc") %>%
    dplyr::group_by(Stratum, StationID, SampleDate, Replicate) %>%
    dplyr::summarise(NumOfMolluscTaxa = sum(n))

  writelog(
    "\nRBI Step 2 - Filter to molluscs and group by Stratum, StationID, Sampledate Replicate - group and get the number of mollusc taxa",
    logfile = logfile,
    code = '
      rbi2 <- rbi_data %>%
        dplyr::filter(Mollusc=="Mollusc") %>%
        dplyr::group_by(Stratum, StationID, SampleDate, Replicate) %>%
        dplyr::summarise(NumOfMolluscTaxa = sum(n))
    ',
    data = rbi2 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi2, logfile = logfile, filename = 'RBI-step2.csv', linktext = 'Download RBI Step 2', verbose = verbose)




  ### SQO RBI -3
  rbi3 <- rbi_data %>%
    dplyr::filter(Crustacean=="Crustacean") %>%
    dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
    dplyr::summarise(NumOfCrustaceanTaxa = sum(n))

  writelog(
    "\nFilter to Crustaceans and group by Stratum, StationID, Sampledate, Replicate - group and get the number of crustacean taxa",
    logfile = logfile,
    code = '
      rbi3 <- rbi_data %>%
        dplyr::filter(Crustacean=="Crustacean") %>%
        dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
        dplyr::summarise(NumOfCrustaceanTaxa = sum(n))
    ',
    data = rbi3 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi3, logfile = logfile, filename = 'RBI-step3.csv', linktext = 'Download RBI Step 3', verbose = verbose)



  ### SQO RBI -4
  rbi4 <- rbi_data %>%
    dplyr::filter(Crustacean=="Crustacean") %>%
    dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
    dplyr::summarise(CrustaceanAbun = sum(Abundance))

  writelog(
    "\nFilter to Crustaceans and group by Stratum, StationID, Sampledate, Replicate - group and get the abundance of crustacean taxa",
    logfile = logfile,
    code = '
      rbi4 <- rbi_data %>%
        dplyr::filter(Crustacean=="Crustacean") %>%
        dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
        dplyr::summarise(CrustaceanAbun = sum(Abundance))
    ',
    data = rbi4 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi4, logfile = logfile, filename = 'RBI-step4.csv', linktext = 'Download RBI Step 4', verbose = verbose)



  ### SQO RBI -5
  rbi5 <- rbi_data %>%
    dplyr::filter(Taxon == "Monocorophium insidiosum") %>%
    dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
    dplyr::summarise(M_insidiosumAbun = sum(Abundance))

  writelog(
    "\nFilter to Monocorophium insidiosum and group by Stratum, StationID, Sampledate, Replicate - group and get the abundance of Monocorophium insidiosum",
    logfile = logfile,
    code = '
      rbi5 <- rbi_data %>%
        dplyr::filter(Taxon == "Monocorophium insidiosum") %>%
        dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
        dplyr::summarise(M_insidiosumAbun = sum(Abundance))
    ',
    data = rbi5 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi5, logfile = logfile, filename = 'RBI-step5.csv', linktext = 'Download RBI Step 5', verbose = verbose)


  ### SQO RBI -6
  rbi6 <- rbi_data %>%
    dplyr::filter(Taxon == "Asthenothaerus diegensis") %>%
    dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
    dplyr::summarise(A_diegensisAbun = sum(Abundance))

  writelog(
    "\nFilter to Asthenothaerus diegensis and group by Stratum, StationID, Sampledate, Replicate - group and get the abundance of Asthenothaerus diegensis",
    logfile = logfile,
    code = '
      rbi6 <- rbi_data %>%
        dplyr::filter(Taxon == "Asthenothaerus diegensis") %>%
        dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
        dplyr::summarise(A_diegensisAbun = sum(Abundance))
    ',
    data = rbi6 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi6, logfile = logfile, filename = 'RBI-step6.csv', linktext = 'Download RBI Step 6', verbose = verbose)



  ### SQO RBI -7
  rbi7 <- rbi_data %>%
    dplyr::filter(Taxon == "Goniada littorea") %>%
    dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
    dplyr::summarise(G_littoreaAbun = sum(Abundance))

  writelog(
    "\nFilter to Goniada littorea and group by Stratum, StationID, Sampledate, Replicate - group and get the abundance of Goniada littorea",
    logfile = logfile,
    code = '
      rbi7 <- rbi_data %>%
        dplyr::filter(Taxon == "Goniada littorea") %>%
        dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
        dplyr::summarise(G_littoreaAbun = sum(Abundance))
    ',
    data = rbi7 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi7, logfile = logfile, filename = 'RBI-step7.csv', linktext = 'Download RBI Step 7', verbose = verbose)



  ### SQO RBI -8
  rbi8 <- rbi_data %>%
    dplyr::filter(Taxon %in% c( "Capitella capitata Cmplx","Oligochaeta")) %>%
    dplyr::mutate(badness=-0.1) %>%
    dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
    dplyr::summarise(NIT = sum(badness))

  writelog(
    "\nFilter to Capitella capitata Cmplx and Oligochaeta and group by Stratum, StationID, Sampledate, Replicate - group and get the sum of 'badness'",
    logfile = logfile,
    code = '
      rbi8 <- rbi_data %>%
        dplyr::filter(Taxon %in% c( "Capitella capitata Cmplx","Oligochaeta")) %>%
        dplyr::mutate(badness=-0.1) %>%
        dplyr::group_by(Stratum, StationID, Replicate, SampleDate) %>%
        dplyr::summarise(NIT = sum(badness))
    ',
    data = rbi8 %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi8, logfile = logfile, filename = 'RBI-step8.csv', linktext = 'Download RBI Step 8', verbose = verbose)



  writelog('\n', logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)

  ### RBI Metrics
  # We are using a full join because if there are missing values, we might just get an empty data frame.
  rbi_metrics <- rbi_data_tmp %>%
    dplyr::full_join(rbi2, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::full_join(rbi3, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::full_join(rbi4, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::full_join(rbi5, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::full_join(rbi6, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::full_join(rbi7, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::full_join(rbi8, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
    dplyr::select(Stratum, StationID, SampleDate, Replicate, NumOfTaxa, NumOfMolluscTaxa, NumOfCrustaceanTaxa, CrustaceanAbun, M_insidiosumAbun, A_diegensisAbun, G_littoreaAbun, NIT)

  writelog(
    "\nRBI Metrics - join the dataframes",
    logfile = logfile,
    code = '
      ### RBI Metrics
      # We are using a full join because if there are missing values, we might just get an empty data frame.
      rbi_metrics <- rbi_data_tmp %>%
        dplyr::full_join(rbi2, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::full_join(rbi3, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::full_join(rbi4, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::full_join(rbi5, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::full_join(rbi6, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::full_join(rbi7, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::full_join(rbi8, by = c("Stratum", "StationID", "Replicate", "SampleDate")) %>%
        dplyr::select(Stratum, StationID, SampleDate, Replicate, NumOfTaxa, NumOfMolluscTaxa, NumOfCrustaceanTaxa, CrustaceanAbun, M_insidiosumAbun, A_diegensisAbun, G_littoreaAbun, NIT)
    ',
    data = rbi_metrics %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_metrics, logfile = logfile, filename = 'RBI-Metrics.csv', linktext = 'Download RBI Metrics dataframe', verbose = verbose)


  ### RBI Category Thresholds for Southern California Marine Bays
  RBI_category_thresholds <- data.frame(
    ref_low = c(0.27, 0.16, 0.08, 0.08),
    ref_high = c(0.27, 0.27, 0.16, 0.08),
    category = as.factor(c("Reference","Low Disturbance","Moderate Disturbance","High Disturbance")),
    category_score = c(1, 2, 3, 4)
  )

  writelog(
    "\nRBI Category Threshold Table",
    logfile = logfile,
    code = '
      ### RBI Category Thresholds for Southern California Marine Bays
      RBI_category_thresholds <- data.frame(
        ref_low = c(0.27, 0.16, 0.08, 0.08),
        ref_high = c(0.27, 0.27, 0.16, 0.08),
        category = as.factor(c("Reference","Low Disturbance","Moderate Disturbance","High Disturbance")),
        category_score = c(1, 2, 3, 4)
      )
    ',
    data = RBI_category_thresholds,
    verbose = verbose
  )

  # Log the RBI calculation code
  rbi_code <- "
    rbi_scores <- rbi_metrics %>%
        replace(., is.na(.), 0) %>%
        mutate(scaled_NumTaxa = NumOfTaxa / 99) %>%
        mutate(scaled_NumMolluscTaxa = NumOfMolluscTaxa / 28) %>%
        mutate(scaled_NumCrustaceanTaxa = NumOfCrustaceanTaxa / 29) %>%
        mutate(scaled_CrustaceanAbun = CrustaceanAbun / 1693) %>%

        # TWV = Taxa Richness Weighted Value
        mutate(TWV = scaled_NumTaxa + scaled_NumMolluscTaxa + scaled_NumCrustaceanTaxa + (0.25 * scaled_CrustaceanAbun)) %>%

        # PIT = Positive Indicator Taxa
        mutate(PIT = ((M_insidiosumAbun)^(1/4) / (473)^(1/4)) + ((A_diegensisAbun)^(1/4) / (27)^(1/4)) + ((G_littoreaAbun)^(1/4) / (15)^(1/4))) %>%

        mutate(Raw_RBI = TWV + NIT + (2 * PIT)) %>%
        dplyr::mutate(Score = (Raw_RBI - 0.03) / 4.69) %>%

        # RBI Categories based on RBI scores
        dplyr::mutate(Category = case_when(
            (Score > 0.27) ~ 'Reference',
            (Score > 0.16 & Score <= 0.27) ~ 'Low Disturbance',
            (Score > 0.08 & Score <= 0.16) ~ 'Moderate Disturbance',
            (Score <= 0.08) ~ 'High Disturbance'
        )) %>%

        # RBI Category Scores based on RBI scores
        dplyr::mutate(`Category Score` = case_when(
            (Category == 'Reference') ~ 1,
            (Category == 'Low Disturbance') ~ 2,
            (Category == 'Moderate Disturbance') ~ 3,
            (Category == 'High Disturbance') ~ 4
        )) %>%
        dplyr::mutate(Index = 'RBI')
    "
  writelog("RBI is very math intensive and difficult to explain with writing - instead I will log the code that calculates the final dataframe", logfile = logfile, verbose = verbose, prefix = hyphen.log.prefix)
  writelog(rbi_code, logfile = logfile, verbose = verbose)

  # Compute the RBI scores.
  # This was not included in the queries that D. Gillet listed. We went through the Technical Manual (p. 69-71)
  # to find the appropriate calculations.

  # NOTE: The above comment ^
  # "This was not included in the queries that D. Gillet listed. We went through the Technical Manual (p. 69-71) to find the appropriate calculations."
  # is from an unknown author. To me it seems that this code may be based on an outdated edition of the manual
  # RBI instructions are on page 62 of the June 2021 edition (https://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf)
  # Scaled values do not seem to match
  # June 2024 - Robert Butler
  rbi_scores <- rbi_metrics %>%
    replace(., is.na(.), 0) %>%
    mutate(scaled_NumTaxa = NumOfTaxa / 99) %>%
    mutate(scaled_NumMolluscTaxa = NumOfMolluscTaxa / 28) %>%
    mutate(scaled_NumCrustaceanTaxa = NumOfCrustaceanTaxa / 29) %>%
    mutate(scaled_CrustaceanAbun = CrustaceanAbun / 1693) %>%

    # TWV = Taxa Richness Weighted Value
    mutate(TWV = scaled_NumTaxa + scaled_NumMolluscTaxa + scaled_NumCrustaceanTaxa + (0.25 * scaled_CrustaceanAbun)) %>%

    # PIT = Positive Indicator Taxa
    mutate(PIT = ((M_insidiosumAbun)^(1/4) / (473)^(1/4)) + ((A_diegensisAbun)^(1/4) / (27)^(1/4)) + ((G_littoreaAbun)^(1/4) / (15)^(1/4))) %>%

    mutate(Raw_RBI = TWV + NIT + (2 * PIT)) %>%
    dplyr::mutate(Score = (Raw_RBI - 0.03) / 4.69) %>%

    # RBI Categories based on RBI scores
    dplyr::mutate(Category = case_when(
      (Score > 0.27) ~ 'Reference',
      (Score > 0.16 & Score <= 0.27) ~ 'Low Disturbance',
      (Score > 0.08 & Score <= 0.16) ~ 'Moderate Disturbance',
      (Score <= 0.08) ~ 'High Disturbance'
    )) %>%

    # RBI Category Scores based on RBI scores
    dplyr::mutate(`Category Score` = case_when(
      (Category == 'Reference') ~ 1,
      (Category == 'Low Disturbance') ~ 2,
      (Category == 'Moderate Disturbance') ~ 3,
      (Category == 'High Disturbance') ~ 4
    )) %>%
    dplyr::mutate(Index = 'RBI')

  writelog(
    "\nRBI Scores Final",
    logfile = logfile,
    code = "
      # Compute RBI Scores
      rbi_scores <- rbi_metrics %>%
        replace(., is.na(.), 0) %>%
        mutate(scaled_NumTaxa = NumOfTaxa / 99) %>%
        mutate(scaled_NumMolluscTaxa = NumOfMolluscTaxa / 28) %>%
        mutate(scaled_NumCrustaceanTaxa = NumOfCrustaceanTaxa / 29) %>%
        mutate(scaled_CrustaceanAbun = CrustaceanAbun / 1693) %>%

        # TWV = Taxa Richness Weighted Value
        mutate(TWV = scaled_NumTaxa + scaled_NumMolluscTaxa + scaled_NumCrustaceanTaxa + (0.25 * scaled_CrustaceanAbun)) %>%

        # PIT = Positive Indicator Taxa
        mutate(PIT = ((M_insidiosumAbun)^(1/4) / (473)^(1/4)) + ((A_diegensisAbun)^(1/4) / (27)^(1/4)) + ((G_littoreaAbun)^(1/4) / (15)^(1/4))) %>%

        mutate(Raw_RBI = TWV + NIT + (2 * PIT)) %>%
        dplyr::mutate(Score = (Raw_RBI - 0.03) / 4.69) %>%

        # RBI Categories based on RBI scores
        dplyr::mutate(Category = case_when(
          (Score > 0.27) ~ 'Reference',
          (Score > 0.16 & Score <= 0.27) ~ 'Low Disturbance',
          (Score > 0.08 & Score <= 0.16) ~ 'Moderate Disturbance',
          (Score <= 0.08) ~ 'High Disturbance'
        )) %>%

        # RBI Category Scores based on RBI scores
        dplyr::mutate(`Category Score` = case_when(
          (Category == 'Reference') ~ 1,
          (Category == 'Low Disturbance') ~ 2,
          (Category == 'Moderate Disturbance') ~ 3,
          (Category == 'High Disturbance') ~ 4
        )) %>%
        dplyr::mutate(Index = 'RBI')
    ",
    data = rbi_scores %>% head(25),
    verbose = verbose
  )
  create_download_link(data = rbi_scores, logfile = logfile, filename = 'RBI-Scores-Final.csv', linktext = 'Download RBI Scores Final dataframe', verbose = verbose)


  writelog('\nEND: RBI function.\n', logfile = logfile, verbose = verbose)

  return(rbi_scores)

}

