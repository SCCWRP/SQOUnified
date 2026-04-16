test_that("chemdata_prep reproduces the expected preprocessed output", {
  rawchem     <- load_fixture("rawchem.rds")
  preppedchem <- load_fixture("preppedchem.rds")

  actual <- chemdata_prep(rawchem)

  # Output is arranged by (stationid, compound) before being renamed to
  # (StationID, AnalyteName) — sort both sides the same way to guarantee
  # stable row order regardless of fixture generation environment.
  key_cols <- c("StationID", "AnalyteName")
  actual   <- actual[do.call(order, actual[key_cols]), , drop = FALSE]
  expected <- preppedchem[do.call(order, preppedchem[key_cols]), , drop = FALSE]
  rownames(actual)   <- NULL
  rownames(expected) <- NULL

  expect_s3_class(actual, "data.frame")
  expect_setequal(names(actual), names(expected))

  # Column order doesn't matter — line them up before comparing.
  actual <- actual[, names(expected), drop = FALSE]

  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("chemdata_prep returns the documented key columns", {
  rawchem <- load_fixture("rawchem.rds")

  actual <- chemdata_prep(rawchem)

  expect_true(all(c("StationID", "AnalyteName", "Result") %in% names(actual)))
  expect_gt(nrow(actual), 0)
})
