#' Compute the benthic response index (BRI) score and BRI condition category.
#'
#' @description
#'   The BRI is the abundance-weighted pollution tolerance score of the organisms present in a benthic sample.
#'   The higher the BRI score, the more degraded the benthic community represented by the sample.
#'
#' @details
#'   Two types of data are needed to calculate the BRI:
#'   (1) the abundance of each species and
#'   (2) its pollution tolerance score, P.
#'
#'   The P values are available for most species present in the assemblage. Only species for which P values are available
#'   are used in the BRI calculations. P values should be obtained for the appropriate habitat and from the most
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
#'   The last step is to compare the BRI score to the BRI threshold values to determine the BRI category and
#'   category score.
#'
#' @param benthic_data a data frame with the following headings
#'
#'    \strong{\code{station_id}} - an alpha-numeric identifier of the location;
#'    \strong{\code{replicate}} - a numeric identifying the replicate number of samples taken at the location;
#'    \strong{\code{sample_date}} - the date of sample collection;
#'    \strong{\code{latitude}} - latitude in decimal degrees;
#'    \strong{\code{longitude}} - longitude in decimal degrees.
#'    Make sure there is a negative sign for the Western coordinates;
#'    \strong{\code{species}} - name of the fauna, ideally in SCAMIT ed12 format. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#' @usage
#' BRI(benthic_data)
#'
#' @examples
#' data(benthic_sampledata) # load sample data
#' BRI(benthic_sampledata) # see the output
#'
#' @importFrom dplyr left_join filter select mutate group_by summarize case_when
#' @importFrom tidyr pivot_longer
#'
#' @export

BRI <- function(benthic_data) {

  out <- benthic_data %>%
    left_join(sqo.list.new, by = c("taxon" = "taxon_name")) %>%
    filter(!is.na(tolerance_score)) %>%
    select(stratum, station_id, sample_date, replicate, taxon, abundance, tolerance_score) %>%
    mutate(
      fourthroot_abun = abundance ^ 0.25,
      weighted_score = fourthroot_abun * tolerance_score
    ) %>%
    group_by(stratum, station_id, sample_date, replicate) %>%
    summarize(
      score = sum(weighted_score, na.rm = TRUE) / sum(fourthroot_abun, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      category = case_when(
        score < 39.96 ~ "Reference",
        score >= 39.96 & score < 49.15 ~ "Low Disturbance",
        score >= 49.15 & score < 73.27 ~ "Moderate Disturbance",
        score >= 73.27 ~ "High Disturbance"
      ),
      category_score = case_when(
        category == "Reference" ~ 1,
        category == "Low Disturbance" ~ 2,
        category == "Moderate Disturbance" ~ 3,
        category == "High Disturbance" ~ 4
      ),
      index = "BRI"
    )

  return(out)
}
