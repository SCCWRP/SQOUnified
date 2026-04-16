library(SQOUnified) 
library(dplyr)
library(openxlsx)

unifiedhost <- Sys.getenv("UNIFIED_DB_HOST")
unifieduser <- Sys.getenv("UNIFIED_DB_USER")
unifiedpass <- Sys.getenv("UNIFIED_DB_PASSWORD")
unifieddbname <- Sys.getenv("UNIFIED_DB_NAME")

unifiedcon <- DBI::dbConnect(
  RPostgres::Postgres(),
  host = unifiedhost,
  user = unifieduser,
  password = unifiedpass,
  dbname = unifieddbname,
  port = 5432,
  sslmode = "require"
)



# replace the stationid with the one that you are evaluating
qry <- paste0("SELECT * FROM vw_unified_chem_sqo_rawdata")

rawchem <- DBI::dbGetQuery(unifiedcon, qry)

preppedchem <- chemdata_prep(rawchem) # prepped chem will be expected out

chem.sqo.out <- chem.sqo(rawchem, verbose = F, knitlog = F)

# --- Write preppedchem values into the SQO Calc Tool Excel file ---
# NOTE: openxlsx only works with .xlsx — re-save the .xls file as .xlsx first
xlpath <- "testing/CalSQOCalcToolVer5.10 (3).xlsx"

wb <- loadWorkbook(xlpath)

# Mapping: preppedchem AnalyteName -> row in the DataInput sheet
analyte_row_map <- c(
  "Cadmium"         = 11,
  "Copper"          = 12,
  "Lead"            = 13,
  "Mercury"         = 14,
  "Zinc"            = 15,
  "HPAH"            = 16,
  "LPAH"            = 17,
  "alpha-Chlordane" = 18,
  "gamma-Chlordane" = 19,
  "Dieldrin"        = 20,
  "trans-Nonachlor" = 21,
  "DDDs_total"      = 22,
  "DDEs_total"      = 23,
  "DDTs_total"      = 24,
  "4,4'-DDT"        = 25,
  "PCBs_total"      = 26
)

# Each station gets its own column: station 1 -> col B (2), station 2 -> col C (3), etc.
for (i in seq_along(stations)) {
  col <- i + 1  # column B=2, C=3, ..., K=11
  stn <- stations[i]
  stn_chem <- preppedchem %>% filter(StationID == stn)

  # Write Station ID in row 6
  writeData(wb, sheet = "DataInput", x = stn, startCol = col, startRow = 6, colNames = FALSE)

  # Write each analyte value
  for (analyte in names(analyte_row_map)) {
    val <- stn_chem$Result[stn_chem$AnalyteName == analyte]
    
    if (length(val) == 1) {
      writeData(wb, sheet = "DataInput", x = val, startCol = col, startRow = analyte_row_map[analyte], colNames = FALSE)
    } else {
      writeData(wb, sheet = "DataInput", x = "x", startCol = col, startRow = analyte_row_map[analyte], colNames = FALSE)
    }
  }

}

saveWorkbook(wb, xlpath, overwrite = TRUE)
