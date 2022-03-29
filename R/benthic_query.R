#' Pull the necessary benthic data from SCCWRP database to compute SQO Benthic indices.
#'
#' @description
#'     This function is to be used internally at SCCWRP.
#'     The naming of the tables is based on how we named tables in our own database.
#'     Notice that any IP addresses, passwords, usernames for a DB Connection must be entered manually.
#'
#'     When this function is run, the latest and greatest benthic data will be pulled in, and a dataframe will be returned.
#'     It will also query the grabevent, field assignment, and station occupation table for other necessary information,
#'     such as LatLongs, or Salinity.
#'
#'     The function returns a dataframe that can be plugged into any of the benthic functions:
#'     RBI, BRI, IBI, RIVPACS, MAMBI, benthic.sqo, and SQOUnified.
#'
#' @param
#'     \strong{db_connection} - a DBI database connection object, created with the DBI::dbConnect function
#'          \link{https://db.rstudio.com/getting-started/connect-to-database/}
#'
#'          If the function is run with no arguments, user will be prompted for database connection information.
#'          hostname or IP, username, password, and database name.
#'
#'
#'
#' @return
#'   \strong{benthic_data} - a dataframe with the following columns
#'
#'     \strong{\code{StationID}} - an alpha-numeric identifier of the location;
#'
#'     \strong{\code{Replicate}} - a numeric identifying the replicate number of samples taken at the location;
#'
#'     \strong{\code{SampleDate}} - the date of sample collection;
#'
#'     \strong{\code{Latitude}} - latitude in decimal degrees;
#'
#'     \strong{\code{Longitude}} - longitude in decimal degrees. Make sure there is a negative sign for the Western coordinates;
#'
#'     \strong{\code{Species}} - name of the fauna, ideally in SCAMIT ed12 format, do not use sp. or spp.,
#'        use sp only or just the Genus. If no animals were present in the sample use
#'        NoOrganismsPresent with 0 abundance;
#'
#'     \strong{\code{Abundance}} - the number of each Species observed in a sample;
#'
#'     \strong{\code{Salinity}} - the salinity observed at the location in PSU, ideally at time of sampling.
#'
#'     \strong{\code{Stratum}} - the stratum of the station (e.g. Port, Bay, Estuary, etc)
#'
#'     \strong{\code{Exclude}} - To be completely honest, I don't know what this means. You'd have to ask Joana or David Gillett.
#'
#' @importFrom DBI dbConnect dbGetQuery
#' @importFrom RPostgreSQL PostgreSQL
#' @importFrom rstudioapi showPrompt askForPassword
#' @import dplyr
#'
#' @export

benthic_query <- function(db_connection = NULL, shelves = T) {
  # con is short for connection
  # Create connection to the database
  if (is.null(db_connection)) {
    con <- DBI::dbConnect(
      PostgreSQL(),
      host = showPrompt("username", "Please enter the hostname of IP address of the database"),
      dbname = showPrompt("dbname", "Please enter the name of the database"),
      user = showPrompt("username", "Please enter the username for the database"),
      askForPassword()
    )
  } else {
    con <- db_connection
  }
  # Bring in our tables from the database ----
  infauna <- dbGetQuery(con, "SELECT * FROM tbl_infaunalabundance_final") %>% as_tibble

  grabqry <- 'SELECT
        stationid, latitude, longitude, depth AS sampledepth, salinity, finalstratum AS stratum
      FROM
        stations_grab_final
      '
  if (!shelves) {
    grabqry <- paste(grabqry, 'WHERE finalstratum IN (\'Bays\', \'Marinas\', \'Ports\', \'Estuaries\', \'Brackish Estuaries\', \'Inner\')')
  }

  grab <- dbGetQuery(con, grabqry) %>% as_tibble()


  # Create the dataset needed to compute all the SQO benthic indices ----
  benthic_data <- infauna  %>%
    inner_join(grab, by = c('stationid')) %>%
    select(
      'stationid','replicate','sampledate','latitude','longitude', 'sampledepth', 'taxon','abundance','salinity', 'stratum', 'exclude'
    ) %>%
    mutate_if(is.numeric, list(~na_if(., -88))) %>%
    rename(
      StationID = stationid,
      Replicate = replicate,
      SampleDate = sampledate,
      Latitude = latitude,
      Longitude = longitude,
      SampleDepth = sampledepth,
      Taxon = taxon,
      Abundance = abundance,
      Salinity = salinity,
      Stratum = stratum,
      Exclude = exclude
    )

  #mutate(EG_Test=ifelse(is.na(EG),"NoEG", "YesEG"))

  return(benthic_data)

}

