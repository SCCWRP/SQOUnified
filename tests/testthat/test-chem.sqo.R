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
