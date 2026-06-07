library(tidyverse)
library(naniar)
library(openxlsx)
library(vegan)
library(janitor)
library(messydates)
library(fuzzyjoin)

source('R/Generic SQO BLOE calculator v2-0.R')



conn <- DBI::dbConnect(
  RPostgres::Postgres(),
  host = "geodbinstance.cottkh4djef2.us-west-2.rds.amazonaws.com",
  port = 5432L,
  dbname = "unified",
  user = "benthic_notebook_reader",
  password = 'SQO_audit_readonly_2026!',
  sslmode = "require"
)
on.exit(DBI::dbDisconnect(conn), add = TRUE)


qry <- '
  SELECT 
    stationid,
    sampledate,
    replicate,
    taxon,
    abundance,
    exclude,
    depth,
    latitude,
    longitude,
    salinity
  FROM 
    sde.staging_benthicinfaunaunifiedpublish
'

rawbenthic <- DBI::dbGetQuery(conn, qry)

infolder <- 'real-bight-benthic-data'
infile <- 'raw-bight-benthic.csv'
outfolder <- 'real-bight-benthic-scores'
file_id = 'bight-benthic-scores.csv'

dir.create(infolder, showWarnings = FALSE, recursive = TRUE)
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

rawbenthic.path <- paste(c(infolder, infile), collapse = '/')
output.scores.path <- paste(c(outfolder, file_id), collapse = '/')

# write it out to a path since the function expects a path
write.csv(rawbenthic, rawbenthic.path)

SQO_BLOE.generic_v2(file_id = file_id, infauna_path = rawbenthic.path, output_path = outfolder)

