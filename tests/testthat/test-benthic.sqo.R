# Tests for the inshore (bay/estuary) Benthic LOE — benthic.sqo()
#
# Fixtures (generated from the unified database by ignore/benthic_unit_test_data_creation.R,
# then written with tests/testthat/fixtures/save_fixtures.R):
#   inshorebenthic       - raw input, from vw_unified_benthic_sqo_rawdata
#   inshore_benthic_sqo  - gold-standard benthic.sqo(inshorebenthic) output
#
# benthic.sqo() returns a named list of data frames (sqo_bloe, sqo_bloe_w_mambi,
# sqo_bri, ibi, rbi, rivpacs, mambi, sqo_potential_mismatches, and the long-format
# all_benthic_sqo_scores_long), so the comparison is element-wise and order-robust
# (see expect_equal_result_list in helper-fixtures.R).



####  expected output is produced by the R package as of version 0.10.18 ###

# test input data is from this query vw_unified_benthic_sqo_rawdata and vw_unified_benthic_bri_offshore_data
# Only difference between these views are the strata. One's where clause will say NOT and the other wont have the NOT
# SELECT
#     b.stationid, b.replicate, b.sampledate, b.latitude, b.longitude, b.depth, b.taxon, b.abundance, b.salinity, st.stratum AS stratum, b.exclude, b.surveyyear
# FROM
# 	-- June 25 2026 - staging_benthicunifiedpublish shall be used, as this is the latest version sent from David
# 	--   supposed to have been confirmed as the accurate benthic taxonomy data by the benthic committee
# 	-- Not moved to "tbl_" yet, as we must wait for the report to be published.
# 	staging_benthicinfaunaunifiedpublish

#     LEFT JOIN {{ source( 'sde', 'lu_bight_stations' ) }} st ON b.stationid = st.stationid
	
# -- mexico stations all had depth of 15 meters or more, indicating that they were likely offshore sites
# -- canyons have really high depth, and are not typically close to the shore. 
# --Technically they are not the same kind of habitat as either embayments or offshore sites, 
# -- but they are closer to offshore sites and therefore lumped into that category

# -- For inshore benthic, remove the "NOT"
# WHERE NOT (
# 	LOWER ( st.stratum ) ~ 'slope' 
# 	OR LOWER( st.stratum ) ~ 'shelf' 
# 	OR LOWER( st.stratum ) ~ 'islands' 
# 	OR LOWER( st.stratum ) ~ 'canyon' 
# 	OR LOWER( st.stratum ) ~ 'mexico' 
# )
# AND UPPER(b.stationid) NOT LIKE '%BOEM%'
# AND UPPER(b.stationid) NOT LIKE '%SONGS%'
# AND UPPER(b.stationid) NOT LIKE '%SAMPLE%'



test_that("benthic.sqo reproduces the expected gold-standard output", {
  inshorebenthic      <- load_fixture("inshorebenthic.rds")
  inshore_benthic_sqo <- load_fixture("inshore_benthic_sqo.rds")

  actual <- benthic.sqo(inshorebenthic, verbose = FALSE, knitlog = FALSE)

  expect_equal_result_list(actual, inshore_benthic_sqo, tolerance = 1e-8)
})

test_that("benthic.sqo output has the documented list structure and BLOE schema", {
  inshorebenthic <- load_fixture("inshorebenthic.rds")

  actual <- benthic.sqo(inshorebenthic, verbose = FALSE, knitlog = FALSE)

  expect_type(actual, "list")
  expect_true(all(c(
    "sqo_bloe", "sqo_bloe_w_mambi", "sqo_bri", "ibi", "rbi", "rivpacs",
    "mambi", "sqo_potential_mismatches", "all_benthic_sqo_scores_long"
  ) %in% names(actual)))

  bloe <- actual$sqo_bloe
  expect_true(all(c(
    "stationid", "sampledate", "replicate",
    "BLOE_score", "BLOE_category", "notes"
  ) %in% names(bloe)))

  # The integrated BLOE score is the 1-4 condition bin.
  expect_true(all(bloe$BLOE_score %in% 1:4 | is.na(bloe$BLOE_score)))

  # Condition labels come from the documented vocabulary.
  expect_true(all(stats::na.omit(bloe$BLOE_category) %in% c(
    "Reference", "Low Disturbance", "Moderate Disturbance", "High Disturbance"
  )))

  # The long-format table tags the integrated row with its full index name.
  long <- actual$all_benthic_sqo_scores_long
  expect_true(all(c("stationid", "sampledate", "replicate", "index") %in% names(long)))
  expect_true("Integrated Benthic LOE SQO Assessment Score" %in% long$index)
})
