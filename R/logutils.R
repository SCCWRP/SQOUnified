#' Function to write to a log file
#'
#' @param content Content to write to the logfile - text or a dataframe/tibble-like object
#' @param logfile Path to the file where the content will be written
#' @param filetype 'txt' or 'csv' - the type of file
#' @param append Append to existing file? or overwrite existing contents
#' @param verbose This is what it will be called in parent functions - verbose tells whether to actually write to the logfile or not
#'                This will allow us to write the code in the other functions and files without a bunch of if statements .
#'
#' @examples
#' writelog(df, 'logs/test.csv', filetype = 'csv', append = T)
#' writelog("# --- Test Log Statement --- #", 'logs/log.txt')

#' @export
writelog <- function(content, logfile, filetype = 'txt', append = T, verbose = T, prefix = NULL, include.row.names = F) {


  # Only text and csv (for now - ...... don't see any reason why we would support any other type at this point)
  if ( !(filetype %in% c('txt','csv')) ) {
    stop("In logging utility function - filetype must be txt or csv")
  }

  if (filetype == 'csv') {
    write.csv(content, logfile, row.names = include.row.names)
  } else {
    if (!file.exists(logfile)) {
      stop(paste0("Error in ",  sys.call(), " - File '", logfile, "' not found!"))
    }
    if (!is.null(prefix)) {
      write(paste(prefix, content), logfile, append = append)
    } else {
      write(content, logfile, append = append)
    }
  }

}



#' Function to write to a log file
#'
#' @param logfile Path to the file where the content will be written
#' @param base.func.name the function name where the call originated - to be written at the top of the log file
#' @param current.time time the function got called
#' @param is.base.func log file shouldnt get initialized if the function gets called within a function that is not the originating one
#'                    Just tells us whether to actually create the log file or not
#' @param verbose This is what it will be called in parent functions - verbose tells whether to actually write to the logfile or not
#'                This will allow us to write the code in the other functions and files without a bunch of if statements .
#'
#' @examples
#' init.log(logfile, sys.call())
#' @export
init.log <- function(logfile, base.func.name, current.time = Sys.time(), is.base.func = T, verbose = T) {

  logfile = path.expand(logfile)
  logdir = dirname(logfile)

  if (!verbose) return(NULL);

  # Create the parent directories if they don't exist
  if (!dir.exists(logdir)) {
    dir.create(logdir, recursive = TRUE)
  }
  # Create the file if it doesn't exist
  if (!file.exists(logfile)) {
    file.create(logfile)
  }

  write(
    paste(
      "----------------------------------------------------------- BEGIN LOG -",
      current.time,
      "----------------------------------------------------------- "
    ),
    logfile
  )

  return(NULL)

}
