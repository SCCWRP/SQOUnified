library(SQOUnified)
library(dplyr)
library(janitor)

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

qry <- "SELECT * FROM tbl_toxresults_unified WHERE species IN ('Eohaustorius estuarius','Mytilus galloprovincialis')"

tox <- DBI::dbGetQuery(unifiedcon, qry) |> 
  filter( !grepl('sample|boem|songs', tolower( as.character(stationid) ) ) ) |>
  filter(!is.null(stationid))


toxsqo_scores <- tox.sqo(tox)

toxsummary <- tox.summary(tox, results.sampletypes = c('Grab'))

