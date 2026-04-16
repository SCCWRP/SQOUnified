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
