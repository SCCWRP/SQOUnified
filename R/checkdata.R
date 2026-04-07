#'
#'
#'
#'
#'
#'
#' @export
checkdata <- function(benthic = NULL, chem = NULL, tox = NULL, stations = NULL, logfile = file.path(getwd(), 'logs', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), 'log.Rmd' ), verbose = F){

  # Initialize Logging
  init.log(logfile, base.func.name = sys.call(), current.time = Sys.time(), is.base.func = length(sys.calls()) == 1, verbose = verbose)


  # check - the LOE must be valid. Benthic, Chemistry or Toxicity
  # certain abbvreviations are allowed


  if (is.null(benthic)){
    warning("Benthic data was not provided.")
  } else {
    if (!is.null(stations)) {
      missing_stations <- setdiff(unique(benthic$stationid), unique(stations$stationid))
      if (length(missing_stations) > 0) {
        stop(
          "The following stationid(s) are in the benthic data but not in the stations data: ",
          paste(missing_stations, collapse = ", ")
        )
      }
    }
  }

  if (is.null(chem)){
    warning("Chemistry data was not provided.")
  } else {
    # perform chemistry checks
  }

  if (is.null(tox)){
    warning("Toxicity data was not provided.")
  } else {
    # perform toxicity checks
  }

  # other checks
  # Make sure the column names of each respective dataframe are correct - i.e. the functions wont blow up
  # For Tox make sure that each sample has the associated control sample within the same dataframe
  # Of course there are others, but we just dont know them at this time



}
