# Documentation for exported package datasets.
# Source files are in data/. Run data-raw/build_data.R to regenerate .rda files.

#' Edition 14 complex taxon name changes (one-to-many type)
#'
#' Maps original SQO taxon names to their SCAMIT Edition 14 replacements for
#' cases where a single old name expands to multiple new names (complex splits).
#' Used by \code{\link{benthicdata_prep}} during taxonomy retrofitting.
#'
#' @format A data frame with 42 rows and 5 columns:
#' \describe{
#'   \item{Original.SQO.Taxon}{Original SQO-compatible taxon name}
#'   \item{change}{Type of taxonomic change}
#'   \item{type}{Change category (e.g. "complex")}
#'   \item{Ed.14.Taxon}{SCAMIT Edition 14 replacement taxon name}
#'   \item{Priority}{Priority order when multiple Edition 14 names map to one SQO name}
#' }
#' @seealso \code{\link{ed14.rollups}}, \code{\link{SoCal.SQO.xls.ed14.link}}, \code{\link{benthicdata_prep}}
#' @source SCCWRP internal reference data. Stored in \code{data/ed14_to_SQO_v2.RData}.
"ed.14.complex"


#' Edition 14 taxon rollup lookup for SQO retrofitting
#'
#' Maps Edition 14 daughter taxa back to their parent SQO-compatible names.
#' Used by \code{\link{benthicdata_prep}} to aggregate modern taxonomy into
#' the coarser SQO name space before index calculations.
#'
#' @format A data frame with 988 rows and 3 columns:
#' \describe{
#'   \item{sqo.name}{SQO-compatible parent taxon name}
#'   \item{ed.14_daughters}{SCAMIT Edition 14 taxon name that rolls up to the SQO parent}
#'   \item{join_level}{Taxonomic level at which the rollup join is applied}
#' }
#' @seealso \code{\link{ed.14.complex}}, \code{\link{SoCal.SQO.xls.ed14.link}}, \code{\link{benthicdata_prep}}
#' @source SCCWRP internal reference data. Stored in \code{data/ed14_to_SQO_v2.RData}.
"ed14.rollups"


#' SQO Excel tool Edition 14 taxon link table
#'
#' Master lookup table linking original SQO taxon names to SCAMIT Edition 14
#' equivalents, with index-specific designations (IBI sensitivity, mollusc/crustacean
#' flags, BRI tolerance score, RIVPACS inclusion). This is the primary reference
#' used by \code{\link{benthicdata_prep}} to retrofit modern taxonomy.
#'
#' @format A data frame with 1218 rows and 9 columns:
#' \describe{
#'   \item{original_sqo_taxon}{Original SQO-compatible taxon name}
#'   \item{change}{Type of taxonomic change applied}
#'   \item{type}{Change category (e.g. "one-to-one", "complex", "rollup")}
#'   \item{ed_14_taxon}{Corresponding SCAMIT Edition 14 taxon name}
#'   \item{ibi_sensitive_taxa}{IBI sensitivity designation ("S" = sensitive)}
#'   \item{mollusc_designation}{Mollusc flag ("Mollusc" or blank)}
#'   \item{crustacean_designation}{Crustacean flag}
#'   \item{bri_tolerance_score}{BRI tolerance score for the taxon}
#'   \item{used_in_rivpacs}{Whether the taxon is included in the RIVPACS model}
#' }
#' @seealso \code{\link{ed.14.complex}}, \code{\link{ed14.rollups}}, \code{\link{benthicdata_prep}}
#' @source SCCWRP internal reference data. Stored in \code{data/ed14_to_SQO_v2.RData}.
"SoCal.SQO.xls.ed14.link"


#' SQO Excel tool taxon lookup list
#'
#' Lookup table used by IBI, BRI, and RBI to assign phylum, class, order, family,
#' species-level flag, tolerance score, and index-specific designations to each taxon.
#' Derived from the Southern California SQO Excel tool maintained by SCCWRP.
#'
#' @format A data frame with 1193 rows and 10 columns:
#' \describe{
#'   \item{TaxonName}{Taxon name (SQO-compatible)}
#'   \item{Phylum}{Phylum classification}
#'   \item{Class}{Class classification}
#'   \item{Order}{Order classification}
#'   \item{Family}{Family classification}
#'   \item{SpeciesLevel}{Whether identification is to species level}
#'   \item{ToleranceScore}{BRI tolerance score}
#'   \item{IBISensitive}{IBI sensitivity designation ("S" = sensitive)}
#'   \item{Mollusc}{Mollusc designation ("Mollusc" or blank)}
#'   \item{Crustacean}{Crustacean designation}
#' }
#' @seealso \code{\link{SoCal.SQO.xls.ed14.link}}, \code{\link{BRI}}, \code{\link{IBI}}, \code{\link{RBI}}
#' @source SCCWRP internal reference data. Stored in \code{data/SoCal SQO xls LU.RData}.
"xl_tool.SoCalLUList"


#' M-AMBI reference standards for saline stations
#'
#' Reference condition values for the M-AMBI index at saline (marine bay) stations
#' in Southern California. Used by \code{\link{MAMBI}} to establish the good/high
#' boundary and bad/bad boundary for index normalisation.
#'
#' @format A data frame with 16 rows and 7 columns:
#' \describe{
#'   \item{stationid}{Station identifier}
#'   \item{replicate}{Replicate number}
#'   \item{sampledate}{Sample collection date}
#'   \item{ambi_score}{AMBI score at the reference station}
#'   \item{S}{Species richness at the reference station}
#'   \item{H}{Shannon diversity (H') at the reference station}
#'   \item{SalZone}{Salinity zone designation}
#' }
#' @seealso \code{\link{TidalFresh_Standards}}, \code{\link{MAMBI}}
#' @source SCCWRP internal reference data. Stored in \code{data/Saline_Standards.RData}.
"Saline_Standards"


#' M-AMBI reference standards for tidal-fresh stations
#'
#' Reference condition values for the M-AMBI index at tidal-freshwater stations
#' in Southern California. Used by \code{\link{MAMBI}} alongside
#' \code{\link{Saline_Standards}} for index normalisation.
#'
#' @format A data frame with 2 rows and 7 columns:
#' \describe{
#'   \item{stationid}{Station identifier}
#'   \item{replicate}{Replicate number}
#'   \item{sampledate}{Sample collection date}
#'   \item{ambi_score}{AMBI score at the reference station}
#'   \item{H}{Shannon diversity (H') at the reference station}
#'   \item{oligo_pct}{Percentage oligochaetes at the reference station}
#'   \item{SalZone}{Salinity zone designation}
#' }
#' @seealso \code{\link{Saline_Standards}}, \code{\link{MAMBI}}
#' @source SCCWRP internal reference data. Stored in \code{data/TidalFresh_Standards.RData}.
"TidalFresh_Standards"


#' Edition 14 taxon-to-p-code lookup for offshore BRI
#'
#' Maps SCAMIT Edition 14 taxon names to p-codes used in the offshore Benthic
#' Response Index (BRI) calculation. Each p-code links to depth-zone tolerance
#' scores in \code{\link{pcodes}}.
#'
#' @format A data frame with 1196 rows and 2 columns:
#' \describe{
#'   \item{p_code}{Integer p-code identifying the tolerance score group}
#'   \item{taxon}{SCAMIT Edition 14 taxon name}
#' }
#' @seealso \code{\link{pcodes}}, \code{\link{BRI.Offshore}}
#' @source SCCWRP internal reference data. Stored in \code{data/pcode_14.RData}.
"ed.14.ptaxa"


#' Depth-zone pollution tolerance scores (p-codes) for offshore BRI
#'
#' Tolerance scores by depth zone for each p-code, used in the offshore Benthic
#' Response Index (BRI) calculation. Taxon-to-p-code assignments are in
#' \code{\link{ed.14.ptaxa}}.
#'
#' @format A data frame with 520 rows and 6 columns:
#' \describe{
#'   \item{p_code}{Integer p-code}
#'   \item{shallow}{Tolerance score for the shallow depth zone (<25 m)}
#'   \item{mid}{Tolerance score for the mid depth zone (35-110 m)}
#'   \item{deep}{Tolerance score for the deep depth zone (130-324 m)}
#'   \item{shallow_mid}{Tolerance score for the shallow-mid overlap zone (25-35 m)}
#'   \item{mid_deep}{Tolerance score for the mid-deep overlap zone (110-130 m)}
#' }
#' @seealso \code{\link{ed.14.ptaxa}}, \code{\link{BRI.Offshore}}
#' @source SCCWRP internal reference data. Stored in \code{data/pcode_14.RData}.
"pcodes"


#' US M-AMBI ecological group reference values (April 2023 update)
#'
#' Ecological group (EG) assignments and reference values for taxa used in the
#' US M-AMBI calculation. Covers taxa from US marine, estuarine, and tidal-fresh
#' habitats. Updated April 2023.
#'
#' @format A data frame with 7185 rows and 10 columns:
#' \describe{
#'   \item{Taxon}{Taxon name}
#'   \item{Exclude}{Whether the taxon is excluded from index calculations}
#'   \item{Hybrid}{Whether the taxon is a hybrid designation}
#'   \item{US}{Ecological group assignment for US waters (overall)}
#'   \item{Standard}{Standard EG value}
#'   \item{US_East}{EG assignment for US East Coast}
#'   \item{US_Gulf}{EG assignment for US Gulf Coast}
#'   \item{US_West}{EG assignment for US West Coast}
#'   \item{Oligochaeta}{Whether the taxon is an oligochaete (relevant for tidal-fresh)}
#'   \item{note}{Additional notes on the EG assignment}
#' }
#' @seealso \code{\link{Saline_Standards}}, \code{\link{TidalFresh_Standards}}, \code{\link{MAMBI}}
#' @source SCCWRP internal reference data. Stored in \code{data/us mambi egvalues 04-23-24.RData}.
"us.mambi.eg.values.04_23_24"


#' RIVPACS reference group assignments for Southern California stations
#'
#' Integer vector of length 72 assigning each of the 72 Southern California
#' reference stations to one of 12 reference groups. Used by \code{\link{RIVPACS}}
#' alongside \code{\link{socal.reference.group.means}} and
#' \code{\link{socal.reference.covariance}} to predict expected fauna.
#'
#' @format An integer vector of length 72 with values 1–12.
#' @seealso \code{\link{socal.reference.group.means}}, \code{\link{socal.reference.covariance}},
#'   \code{\link{socal.reference.taxa}}, \code{\link{RIVPACS}}
#' @source SCCWRP internal reference data. Stored in \code{data/SoCalReference.RData}.
"socal.reference.groups"


#' RIVPACS reference group centroid means for Southern California
#'
#' A 12-by-3 matrix of environmental predictor means for each of the 12
#' Southern California RIVPACS reference groups. Columns are
#' \code{Latitude}, \code{Longitude}, and \code{SampleDepth}. Used with
#' \code{\link{socal.reference.covariance}} to compute Mahalanobis distances
#' when assigning a new station to reference groups.
#'
#' @format A numeric matrix with 12 rows (one per reference group) and 3 columns:
#' \describe{
#'   \item{Latitude}{Mean latitude of the reference group}
#'   \item{Longitude}{Mean longitude of the reference group}
#'   \item{SampleDepth}{Mean sample depth (m) of the reference group}
#' }
#' @seealso \code{\link{socal.reference.groups}}, \code{\link{socal.reference.covariance}},
#'   \code{\link{socal.reference.taxa}}, \code{\link{RIVPACS}}
#' @source SCCWRP internal reference data. Stored in \code{data/SoCalReference.RData}.
"socal.reference.group.means"


#' RIVPACS pooled covariance matrix for Southern California environmental predictors
#'
#' A 3-by-3 pooled covariance matrix for the environmental predictors
#' (\code{Latitude}, \code{Longitude}, \code{SampleDepth}) used to compute
#' Mahalanobis distances in the Southern California RIVPACS model.
#'
#' @format A numeric matrix with 3 rows and 3 columns
#'   (\code{Latitude}, \code{Longitude}, \code{SampleDepth}).
#' @seealso \code{\link{socal.reference.groups}}, \code{\link{socal.reference.group.means}},
#'   \code{\link{socal.reference.taxa}}, \code{\link{RIVPACS}}
#' @source SCCWRP internal reference data. Stored in \code{data/SoCalReference.RData}.
"socal.reference.covariance"


#' RIVPACS reference station taxa abundances for Southern California
#'
#' A data frame of taxa abundances across the 72 Southern California RIVPACS
#' reference stations. Each row is a reference station; each column is a taxon.
#' Used by \code{\link{RIVPACS}} to compute expected taxon probabilities
#' (O/E ratio) for new stations.
#'
#' @format A data frame with 72 rows (reference stations) and 457 columns (taxa).
#'   Each cell contains the abundance of the column taxon at the row station.
#' @seealso \code{\link{socal.reference.groups}}, \code{\link{socal.reference.group.means}},
#'   \code{\link{socal.reference.covariance}}, \code{\link{RIVPACS}}
#' @source SCCWRP internal reference data. Stored in \code{data/SoCalReference.RData}.
"socal.reference.taxa"
