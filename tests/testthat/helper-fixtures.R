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
