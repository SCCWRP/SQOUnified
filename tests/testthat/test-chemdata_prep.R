# Test data for chem from vw_unified_chem_sqo_rawdata
#  SELECT chem.stationid,
#   st.latitude,
#   st.longitude,
#   st.stratum,
#   chem.analytename,
#   chem.result,
#   chem.mdl,
#   chem.rl,
#   chem.units,
#   chem.sampletype AS sampletypecode,
#   chem.labreplicate AS labrep,
#       CASE
#           WHEN chem.fieldduplicate = 0 THEN 1
#           ELSE chem.fieldduplicate
#       END AS fieldrep,
#   chem.surveyyear
#  FROM tbl_chemistryunifiedpublish chem
#    LEFT JOIN lu_bight_stations st ON chem.stationid::text = st.stationid::text
# WHERE chem.sampletype::text = 'Result'::text AND chem.labreplicate = 1 AND (chem.fieldduplicate = ANY (ARRAY[0, 1])) AND chem.analytename::text <> 'Fines'::text AND lower(chem.analytename::text) !~~ 'phi%'::text AND chem.analytename::text !~ '^SEM'::text AND (chem.stationid::text <> ALL (ARRAY['B23-B1'::character varying::text, 'B23-SONGS Unit 2'::character varying::text, 'B23-12148'::character varying::text, 'B23-12150'::character varying::text, 'B23-13260'::character varying::text, 'B23-12363'::character varying::text, 'B23-12364'::character varying::text, 'B23-12366'::character varying::text, 'B23-12368'::character varying::text, 'B23-12369'::character varying::text]))
# ORDER BY chem.surveyyear, chem.stationid, chem.analytename, chem.result, chem.units;

# expected output is produced by the R package as of version 0.9.5

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
