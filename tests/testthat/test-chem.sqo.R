test_that("chem.sqo reproduces the expected final SQO output", {
  rawchem      <- load_fixture("rawchem.rds")
  chem_sqo_out <- load_fixture("chem_sqo_out.rds")

  actual <- chem.sqo(rawchem, verbose = FALSE, knitlog = FALSE)

  # chem.sqo returns combined.final arranged by (StationID, Index); enforce
  # the same order on both sides to be robust to locale/collation.
  key_cols <- c("StationID", "Index")
  actual   <- actual[do.call(order, actual[key_cols]), , drop = FALSE]
  expected <- chem_sqo_out[do.call(order, chem_sqo_out[key_cols]), , drop = FALSE]
  rownames(actual)   <- NULL
  rownames(expected) <- NULL

  expect_s3_class(actual, "data.frame")
  expect_setequal(names(actual), names(expected))

  actual <- actual[, names(expected), drop = FALSE]

  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("chem.sqo output has the documented score schema", {
  rawchem <- load_fixture("rawchem.rds")

  actual <- chem.sqo(rawchem, verbose = FALSE, knitlog = FALSE)

  expect_true(all(
    c("StationID", "Index", "Score", "Category", "Category Score") %in% names(actual)
  ))

  # Every station gets three Index rows: LRM, CSI, and the integrated score.
  expect_setequal(
    unique(actual$Index),
    c("LRM", "CSI", "Integrated Chemistry LOE SQO Assessment Score")
  )

  # `Category Score` is the 1–4 bin for every row. `Score` is the raw
  # continuous index value for LRM/CSI rows, and equals `Category Score`
  # only for the Integrated row.
  expect_true(all(actual$`Category Score` %in% 1:4 | is.na(actual$`Category Score`)))
  expect_true(all(
    actual$Category %in% c(
      "Minimal Exposure", "Low Exposure", "Moderate Exposure", "High Exposure", NA
    )
  ))

  integrated <- actual[actual$Index == "Integrated Chemistry LOE SQO Assessment Score", ]
  expect_true(all(integrated$Score %in% 1:4 | is.na(integrated$Score)))
  expect_equal(integrated$Score, integrated$`Category Score`)
})
