# Tests for the offshore Benthic Response Index (BRI.Offshore)
# offshore_data is a joined benthic+station fixture defined in helper-fixtures.R

test_that("BRI.Offshore returns the documented list structure", {
  res <- BRI.Offshore(offshore_data)

  expect_type(res, "list")
  expect_named(
    res,
    c("bri_scores", "taxa_with_pcode", "taxa_without_pcode", "all_taxa_by_sample"),
    ignore.order = TRUE
  )
  expect_s3_class(res$bri_scores, "data.frame")
  expect_s3_class(res$taxa_with_pcode, "data.frame")
  expect_s3_class(res$taxa_without_pcode, "data.frame")
  expect_s3_class(res$all_taxa_by_sample, "data.frame")
})

test_that("BRI.Offshore (wide) has the documented schema and never drops a sample", {
  res <- BRI.Offshore(offshore_data)
  scores <- res$bri_scores

  expect_true(all(c(
    "stationid", "sampledate", "replicate", "depth", "latitude", "longitude",
    "depth_zone", "bri_score", "bri_cond", "bri_class", "note"
  ) %in% names(scores)))

  # The skeleton-join must retain every submitted sample (no silent drops).
  n_in <- length(unique(paste(
    offshore_data$stationid,
    offshore_data$sampledate,
    offshore_data$replicate
  )))
  expect_equal(nrow(scores), n_in)

  # depth is carried on every row (regression guard for the old NA-depth bug)
  expect_false(any(is.na(scores$depth)))

  # condition class is in the documented 1-5 range (NA allowed for "BRI not Applicable")
  expect_true(all(scores$bri_class %in% c(1, 2, 3, 4, 5, NA)))

  # condition labels are from the documented vocabulary
  expect_true(all(stats::na.omit(scores$bri_cond) %in% c(
    "Reference", "Marginal Deviation", "Biodiversity Loss",
    "Function Loss", "Defaunation", "BRI not Applicable"
  )))
})

test_that("BRI.Offshore condition categories match the score thresholds", {
  scores <- BRI.Offshore(offshore_data)$bri_scores
  scored <- scores[!is.na(scores$bri_score) & !is.nan(scores$bri_score), ]

  expected_class <- cut(
    scored$bri_score,
    breaks = c(-Inf, 25, 34, 44, 72, Inf),
    labels = c(1, 2, 3, 4, 5),
    right = FALSE
  )
  expect_equal(as.integer(as.character(expected_class)), as.numeric(scored$bri_class))
})

test_that("BRI.Offshore long format exposes the columns SQOUnified consumes", {
  long <- BRI.Offshore(offshore_data, output_format = "long")$bri_scores

  expect_true(all(c("stationid", "index", "score", "category", "category_score") %in% names(long)))
  expect_true(all(long$index == "Offshore BRI"))
})

test_that("BRI.Offshore taxa_without_pcode has a stable schema", {
  tw <- BRI.Offshore(offshore_data)$taxa_without_pcode

  expect_named(tw, c("submitted_taxon", "did_you_mean_taxon", "p_code"), ignore.order = TRUE)
  # the NoOrganismsPresent sentinel is not reported as an unmatched taxon
  expect_false("NoOrganismsPresent" %in% tw$submitted_taxon)
})

test_that("BRI.Offshore rejects an invalid output_format", {
  expect_error(
    BRI.Offshore(offshore_data, output_format = "tall"),
    "Invalid output format"
  )
})
