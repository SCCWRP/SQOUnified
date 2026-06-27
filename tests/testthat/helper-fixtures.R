# Combined offshore sample data used by test-BRI.Offshore.R
.offshore_env <- new.env()
data(offshore_bri_sampledata, package = "SQOUnified", envir = .offshore_env)
offshore_data <- .offshore_env$offshore_bri_sampledata
rm(.offshore_env)

load_fixture <- function(name) {
  path <- testthat::test_path("fixtures", name)
  if (!file.exists(path)) {
    testthat::skip(paste0(
      "Fixture '", name, "' not found. ",
      "Source tests/testthat/fixtures/save_fixtures.R after running ",
      "scripts/create_test_data.R to regenerate fixtures."
    ))
  }
  readRDS(path)
}

# --- Gold-standard list comparison -------------------------------------------
# benthic.sqo() and BRI.Offshore() each return a *named list of data frames*
# rather than a single table. These helpers compare an actual result list to a
# saved expected one in a way that is robust to row order: the machine that
# generated the fixture and the machine running the test can disagree on
# locale/collation and on dplyr group ordering, so we sort both sides the same
# way before comparing.

# Sort a data frame by whatever identifier columns it carries so two copies with
# the same rows in a different order line up for an element-wise comparison.
sort_for_compare <- function(df) {
  key_candidates <- c(
    "sample_id", "stationid", "StationID", "sampledate", "SampleDate",
    "replicate", "Replicate", "index", "Index", "taxon", "Taxon",
    "submitted_taxon", "p_code"
  )
  keys <- intersect(key_candidates, names(df))
  if (length(keys) == 0) keys <- names(df)
  # method = "radix" sorts in the C locale, so the order is identical regardless
  # of the locale the fixture was generated in.
  ord <- do.call(order, c(unname(as.list(df[keys])), list(method = "radix")))
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL
  # Some columns (e.g. M-AMBI's Shannon "H") arrive as *named* numeric vectors;
  # whether those element names survive is environment-dependent and is not
  # tabular data, so strip them so a stray names attribute can't fail the compare.
  df[] <- lapply(df, function(col) {
    names(col) <- NULL
    col
  })
  df
}

# Compare two named lists of data frames produced by the benthic SQO functions.
# For each element: column *names* must match (order need not), then rows are
# sorted the same way on both sides and compared with the given tolerance.
expect_equal_result_list <- function(actual, expected, tolerance = 1e-8) {
  testthat::expect_type(actual, "list")
  testthat::expect_setequal(names(actual), names(expected))

  for (nm in names(expected)) {
    act <- actual[[nm]]
    exp <- expected[[nm]]

    if (is.data.frame(exp)) {
      testthat::expect_s3_class(act, "data.frame")
      testthat::expect_setequal(names(act), names(exp))
      # align column order to the expected frame, then sort rows identically
      act <- sort_for_compare(act[, names(exp), drop = FALSE])
      exp <- sort_for_compare(exp)
    }

    testthat::expect_equal(act, exp, tolerance = tolerance, info = paste0("list element: ", nm))
  }
}
