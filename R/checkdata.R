#'
#'
#'
#'
#'
#'
#' @export
checkdata <- function(benthic = NULL, chem = NULL, tox = NULL, logfile = file.path(getwd(), 'logs', paste0(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), '-log.txt') ), verbose = T){

  # check - the LOE must be valid. Benthic, Chemistry or Toxicity
  # certain abbvreviations are allowed


  if (is.null(benthic)){
    warning("Benthic data was not provided.")
  } else {
    # perform benthic checks
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
