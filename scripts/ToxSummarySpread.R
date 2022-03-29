library(DBI) # needed to connect to database
library(dbplyr) # needed to connect to database
library(RPostgreSQL) # needed to connect to our database
library(rstudioapi) # just so we can type the password as we run the script, so it is not written in the clear
library(tidyverse)

# con is short for connection
# Create connection to the database
con <- DBI::dbConnect(
  PostgreSQL(),
  host = "192.168.1.16",
  dbname = 'bight2018',
  user = 'b18read',
  password = '1969$Harbor' # if we post to github, we might want to do rstudioapi::askForPassword()
)

tox_sampledata <- dbGetQuery(con, '
  SELECT 
    tbl_toxresults.stationid,
    tbl_toxresults.toxbatch,
    tbl_toxresults.species,
    tbl_toxresults.sampletypecode,
    tbl_toxresults.matrix,
    tbl_toxresults.fieldreplicate,
    tbl_toxresults.labrep,
    tbl_toxresults.result
  FROM
    tbl_toxresults 
  WHERE
    tbl_toxresults.sampletypecode NOT IN  (\'QA\')
  AND
    tbl_toxresults.matrix != \'Reference Toxicant\'
    ;
  ') %>% as_tibble()

demo_tox_controls <- dbGetQuery(con, '
  SELECT 
    tbl_toxresults.stationid,
    tbl_toxresults.toxbatch,
    tbl_toxresults.species,
    tbl_toxresults.sampletypecode,
    tbl_toxresults.matrix,
    tbl_toxresults.fieldreplicate,
    tbl_toxresults.labrep,
    tbl_toxresults.result
  FROM
    tbl_toxresults 
  WHERE
    tbl_toxresults.stationid = \'0000\'
  AND
    tbl_toxresults.sampletypecode = \'CNEG\'
    ;
  ') 


summary <- dbGetQuery(con, '
  SELECT 
    tbltoxicitysummaryresults.stationid, 
    tbltoxicitysummaryresults.latitude, 
    tbltoxicitysummaryresults.longitude, 
    tbltoxicitysummaryresults.stratum, 
    field_assignment_table.gisdepth as depth, 
    field_assignment_table.region,
    tbltoxicitysummaryresults.lab,
    tbltoxicitysummaryresults.sampletypecode,
    tbltoxicitysummaryresults.sqocategory,
    tbltoxicitysummaryresults.species,
    tbltoxicitysummaryresults.pctcontrol
  FROM
    tbltoxicitysummaryresults 
  LEFT JOIN 
    field_assignment_table
  ON
    tbltoxicitysummaryresults.stationid = field_assignment_table.stationid
  WHERE
    tbltoxicitysummaryresults.stratum != \'QC\'
  AND
    tbltoxicitysummaryresults.sampletypecode != \'QA\'
    ;
  ') %>%
  dplyr::arrange(stationid, stratum, lab, sampletypecode)
save(tox_sampledata, file = "data/tox_sampledata.RData")

summary <- summary %>%
  dplyr::select(-lab) %>%
  dplyr::mutate(
    IntegratedSQO = dplyr::case_when(
      sqocategory == 'Nontoxic' ~ 1,
      sqocategory == 'Low Toxicity' ~ 2,
      sqocategory == 'Moderate Toxicity' ~ 3,
      sqocategory == 'High Toxicity' ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::group_by(stationid) %>%  
  dplyr::mutate(
    `Integrated Toxicity Category` = mean(IntegratedSQO, na.rm = FALSE)
  ) %>%
  dplyr::select(-c(sqocategory, IntegratedSQO)) %>%
  tidyr::spread(key = species, value = pctcontrol) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    # We use ceiling because rounding in R would round 2.5 down to 2, since it goes to the nearest even number
    `Integrated Toxicity Category` = ceiling(`Integrated Toxicity Category`)
  ) %>%
  dplyr::select(
    stationid, latitude, longitude, stratum, depth, region, 
    `Eohaustorius estuarius`, `Mytilus galloprovincialis`, 
    `Integrated Toxicity Category`
  ) %>% 
  dplyr::rename(
    Station = stationid,
    `Latitude (north)` = latitude,
    `Longitude (west)` = longitude,
    Stratum = stratum,
    `Depth (m)` = depth,
    Region = region,
    `Amphipod Survival (% control)` = `Eohaustorius estuarius`,
    `Mussel Embryo Percent Normal Alive (% control)` = `Mytilus galloprovincialis`
  ) 

write.csv(summary, "ToxResultsByStation.csv", row.names=F)


