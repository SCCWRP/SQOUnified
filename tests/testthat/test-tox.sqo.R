test_that("tox.sqo reproduces the expected SQO output", {
  tox           <- load_fixture("tox.rds")
  toxsqo_scores <- load_fixture("toxsqo_scores.rds")

  actual <- tox.sqo(tox, verbose = FALSE, knitlog = FALSE)

  # tox.sqo arranges by (StationID, Index, Category) before return; enforce
  # the same order on both sides so a locale/collation quirk won't fail the test.
  key_cols <- c("StationID", "Index", "Category")
  actual   <- actual[do.call(order, actual[key_cols]), , drop = FALSE]
  expected <- toxsqo_scores[do.call(order, toxsqo_scores[key_cols]), , drop = FALSE]
  rownames(actual)   <- NULL
  rownames(expected) <- NULL

  expect_s3_class(actual, "data.frame")
  expect_setequal(names(actual), names(expected))

  actual <- actual[, names(expected), drop = FALSE]

  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("tox.sqo output has the documented score schema", {
  tox <- load_fixture("tox.rds")

  actual <- tox.sqo(tox, verbose = FALSE, knitlog = FALSE)

  expect_true(all(
    c("StationID", "Index", "Score", "Category", "Category Score") %in% names(actual)
  ))

  # Category Score is the 1–4 bin for every row.
  expect_true(all(actual$`Category Score` %in% 1:4 | is.na(actual$`Category Score`)))

  # Tox categories per the manual.
  expect_true(all(
    actual$Category %in% c(
      "Nontoxic", "Low Toxicity", "Moderate Toxicity", "High Toxicity", NA
    )
  ))

  # On the integrated row Score equals Category Score.
  integrated <- actual[actual$Index == "Integrated Toxicity LOE SQO Assessment Score", ]
  expect_gt(nrow(integrated), 0)
  expect_equal(integrated$Score, integrated$`Category Score`)
})
