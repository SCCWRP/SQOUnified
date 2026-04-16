## Run this ONCE from the R session that has rawchem, preppedchem,
## and chem.sqo.out already loaded (e.g. right after scripts/create_test_data.R).
## It writes the three fixtures the testthat suite loads.

fixtures_dir <- file.path("tests", "testthat", "fixtures")
if (!dir.exists(fixtures_dir)) {
  stop(
    "Could not find '", fixtures_dir, "'. ",
    "Run this script with the package root as your working directory."
  )
}

required <- c("rawchem", "preppedchem", "chem.sqo.out")
missing  <- setdiff(required, ls(envir = .GlobalEnv))
if (length(missing) > 0) {
  stop(
    "Missing objects in the global environment: ",
    paste(missing, collapse = ", "),
    ". Run scripts/create_test_data.R first."
  )
}

saveRDS(get("rawchem",      envir = .GlobalEnv), file.path(fixtures_dir, "rawchem.rds"))
saveRDS(get("preppedchem",  envir = .GlobalEnv), file.path(fixtures_dir, "preppedchem.rds"))
saveRDS(get("chem.sqo.out", envir = .GlobalEnv), file.path(fixtures_dir, "chem_sqo_out.rds"))

message("Wrote fixtures to ", normalizePath(fixtures_dir))
