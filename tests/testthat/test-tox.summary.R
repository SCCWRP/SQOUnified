test_that("tox.summary reproduces the expected summary output", {
  tox        <- load_fixture("tox.rds")
  toxsummary <- load_fixture("toxsummary.rds")

  # Matches the call used in scripts/create_tox_sampledata.R.
  actual <- tox.summary(tox, results.sampletypes = c("Grab"), verbose = FALSE, knitlog = FALSE)

  # Sort both sides by a stable composite key.
  key_cols <- intersect(
    c("lab", "stationid", "species", "toxbatch", "sampletypecode"),
    names(actual)
  )
  actual   <- actual[do.call(order, actual[key_cols]), , drop = FALSE]
  expected <- toxsummary[do.call(order, toxsummary[key_cols]), , drop = FALSE]
  rownames(actual)   <- NULL
  rownames(expected) <- NULL

  expect_s3_class(actual, "data.frame")
  expect_setequal(names(actual), names(expected))

  actual <- actual[, names(expected), drop = FALSE]

  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("tox.summary output has the documented schema", {
  tox <- load_fixture("tox.rds")

  actual <- tox.summary(tox, results.sampletypes = c("Grab"), verbose = FALSE, knitlog = FALSE)

  expect_true(all(
    c(
      "stationid", "species", "toxbatch", "sampletypecode",
      "P Value", "Mean", "Control Adjusted Mean",
      "Score", "Category"
    ) %in% names(actual)
  ))

  expect_true(all(actual$Score %in% 1:4 | is.na(actual$Score)))
  expect_true(all(
    actual$Category %in% c(
      "Nontoxic", "Low Toxicity", "Moderate Toxicity", "High Toxicity", NA
    )
  ))
})
