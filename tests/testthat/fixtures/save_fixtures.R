## Run this from the R session that has the fixture objects loaded.
## Chem set:     rawchem, preppedchem, chem.sqo.out
##               (sourced from scripts/create_test_data.R)
## Tox set:      tox, toxsqo_scores, toxsummary
##               (sourced from scripts/create_tox_sampledata.R)
## Benthic set:  inshorebenthic, inshore_benthic_sqo, offshorebenthic, offshore_bri_output
##               (sourced from ignore/benthic_unit_test_data_creation.R)
##
## Each set is saved independently — you can run this after any script
## and only the objects that are loaded will be written.

fixtures_dir <- file.path("tests", "testthat", "fixtures")
if (!dir.exists(fixtures_dir)) {
  stop(
    "Could not find '", fixtures_dir, "'. ",
    "Run this script with the package root as your working directory."
  )
}

fixture_sets <- list(
  chem = list(
    rawchem        = "rawchem.rds",
    preppedchem    = "preppedchem.rds",
    `chem.sqo.out` = "chem_sqo_out.rds"
  ),
  tox = list(
    tox           = "tox.rds",
    toxsqo_scores = "toxsqo_scores.rds",
    toxsummary    = "toxsummary.rds"
  ),
  benthic = list(
    inshorebenthic      = "inshorebenthic.rds",
    inshore_benthic_sqo = "inshore_benthic_sqo.rds",
    offshorebenthic     = "offshorebenthic.rds",
    offshore_bri_output = "offshore_bri_output.rds"
  )
)

loaded <- ls(envir = .GlobalEnv)
written <- character(0)
skipped_sets <- character(0)

for (set_name in names(fixture_sets)) {
  mapping <- fixture_sets[[set_name]]
  missing <- setdiff(names(mapping), loaded)

  if (length(missing) == length(mapping)) {
    skipped_sets <- c(skipped_sets, set_name)
    next
  }
  if (length(missing) > 0) {
    stop(
      "Partial '", set_name, "' fixture set in the global environment. Missing: ",
      paste(missing, collapse = ", "),
      ". Rerun the matching create_*.R script first."
    )
  }

  for (obj_name in names(mapping)) {
    out_path <- file.path(fixtures_dir, mapping[[obj_name]])
    saveRDS(get(obj_name, envir = .GlobalEnv), out_path)
    written <- c(written, out_path)
  }
}

if (length(written) == 0) {
  stop(
    "No fixture objects found in the global environment. ",
    "Source scripts/create_test_data.R and/or scripts/create_tox_sampledata.R first."
  )
}

message("Wrote ", length(written), " fixture file(s) to ", normalizePath(fixtures_dir))
for (p in written) message("  - ", basename(p))
if (length(skipped_sets) > 0) {
  message("Skipped (objects not in global env): ", paste(skipped_sets, collapse = ", "))
}
