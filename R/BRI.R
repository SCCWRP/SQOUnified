#' Compute the benthic response index (BRI) score and BRI condition category.
#'
#' @description
#'   The BRI is the abundance weighted pollution tolerance score of the organisms present in a benthic sample. The higher
#'   the BRI score, the more degraded the benthic community represented by the sample.
#'
#' @details
#'   The BRI is the abundance weighted pollution tolerance score of the organisms present in a benthic sample. The higher
#'   the BRI score, the more degraded the benthic community represented by the sample.
#'
#'   Two types of data are needed to calculate the BRI:
#'
#'   (1) the abundance of each species and
#'   (2) its pollution tolerance score, P.
#'
#'   The P values are available for most species present in the assemblage. Only species for which P values are avialable
#'   are used in the BRI calculations. P values showld be obtained for the appropriate habitat and from the most
#'   up-to-date list available.
#'
#'   The first step in the BRI calculation is to compute the 4th root of the abundance of each taxon in the sample for
#'   which P values are available. The next step is to multiply the 4th root abundance value by the P value, for each
#'   taxon.
#'
#'   Next, separately sum all of the 4th roots of the abundances and all of the products of the 4th roots of abundance
#'   and P values. Taxa that lack P values are not included in either sum. The next step is to calculate the BRI score
#'   as:
#'
#'   \deqn{ \frac{\sum \left(\sqrt[p]{\textrm{Abundance}} \right) \times P}{\sum \sqrt[p]{\textrm{Abundance}}} }
#'
#'   The last step is to compare the BRI score to the BRI threshold values in Table 5 to determine the BRI category and
#'   category score.
#'
#'   <Table 5. To be included in R markdown file>
#'
#'
#'
#' @param BenthicData a data frame with the following headings
#'
#'    \strong{\code{StationID}} - an alpha-numeric identifier of the location;
#'
#'    \strong{\code{Replicate}} - a numeric identifying the replicate number of samples taken at the location;
#'
#'    \strong{\code{SampleDate}} - the date of sample collection;
#'
#'    \strong{\code{Latitude}} - latitude in decimal degrees;
#'
#'    \strong{\code{Longitude}} - longitude in decimal degrees.
#'    Make sure there is a negative sign for the Western coordinates;
#'
#'    \strong{\code{Species}} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#' @usage
#' BRI(benthic_data)
#'
#' @examples
#' data(benthic_sampledata) # load sample data
#' BRI(benthic_sampledata) # see the output
#'
#' @import vegan
#' @import reshape2
#' @importFrom dplyr left_join filter rename select mutate group_by summarize summarise case_when
#'
#' @export


BRI <- function(BenthicData, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), 'log.txt' ), verbose = T)
{

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)
  hyphen.log.prefix <- rep('-', (2 * (length(sys.calls))) - 1)

  out <- BenthicData %>%
  left_join(sqo.list.new, by = c('Taxon' = 'TaxonName')) %>%
  #dplyr::right_join(assignment, by = 'stationid') %>%
  # I assume that the next line is something they had in there as a method of removing duplicates
  # for this reason, this next line will likely be eliminated.
  # They grouped by all the columns that were selected (In query BRI - 1)
  # Instead, if need be we can use something from dplyr that deals with duplicates
  # I actually found that it didn't appear to make a difference
  #dplyr::group_by(stratum, stationid, replicate, taxon, abundance, `B-CodeScore`) %>%
  #dplyr::filter(B13_Stratum %in% c("Estuaries", "Marinas", "Bays", "Ports")) %>%
  filter(!is.na(ToleranceScore)) %>%
  #rename(Stratum) %>%
  select(Stratum, StationID, SampleDate, Replicate, Taxon, Abundance, ToleranceScore)  %>%
  # End of BRI - 1 query. Begin BRI - 2 query
  mutate(
    fourthroot_abun = Abundance ** 0.25,
    tolerance_score = fourthroot_abun * ToleranceScore
  ) %>%
  # End of BRI - 2. Begin BRI - 3
  group_by(
    Stratum, StationID, SampleDate, Replicate
  ) %>%
  summarize(
    Score = sum(tolerance_score, na.rm = T) / sum(fourthroot_abun, na.rm = T)
  ) %>%
    # Output the BRI category given the BRI score and the thresholds for Southern California Marine Bays
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

  return(out)
}


