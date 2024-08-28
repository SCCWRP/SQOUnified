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
writelog <- function(content, logfile, filetype = 'text', append = T, verbose = F, prefix = NULL, include.row.names = F, code = NULL, data = NULL, echo.code = TRUE, pageLength = 10) {

  if (!verbose) return(NULL)

  # Only text and csv (for now - ...... don't see any reason why we would support any other type at this point)
  if ( !(filetype %in% c('text','csv')) ) {
    stop("In logging utility function - filetype must be text or csv")
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

    # Mainly for R Markdown - write code snippets to be able to view - also data table sections
    if (!is.null(code)) {
      if (echo.code) {
        write("```{r}\n", file = logfile, append = append)
      } else {
        write("```{r echo=FALSE}\n", file = logfile, append = append)
      }
      write(code, file = logfile, append = append)
      write("```\n", file = logfile, append = append)
    }
    if (!is.null(data)) {
      write("```{r echo=FALSE}\n", file = logfile, append = append)

      data_name <- deparse(substitute(data))
      write(paste0("datatable(", data_name, ", options = list(pageLength = ", pageLength, ", autoWidth = TRUE))\n"), file = logfile, append = append)

      write("```\n", file = logfile, append = append)
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
init.log <- function(logfile, base.func.name, type = 'text', current.time = Sys.time(), is.base.func = T, verbose = F, title = 'Log', libraries = c('rmarkdown')) {

  logfile = path.expand(logfile)
  logdir = dirname(logfile)

  if ( !(type %in% c('RMarkdown','text') )) {
    stop("Logfile type must be RMarkdown or text")
  }

  if (!verbose) return(NULL);

  # Create the parent directories if they don't exist
  if (!dir.exists(logdir)) {
    dir.create(logdir, recursive = TRUE)
  }
  # Create the file if it doesn't exist
  if (!file.exists(logfile)) {
    file.create(logfile)
    if (type == 'text') {
      write(
        paste(
          "----------------------------------------------------------- BEGIN LOG -",
          current.time,
          "----------------------------------------------------------- "
        ),
        logfile
      )
    } else {
      write( paste0("---\ntitle: \"", title ,"\"\noutput: html_document\n---", collapse = ''), file = logfile)
      write('```{r setup, include=FALSE, message=FALSE, warning=FALSE}\n', file = logfile, append = TRUE)
      for (lib in libraries) {
        write(paste0('library(', lib, ')'), file = logfile, append = TRUE)
      }
      write('```', file = logfile, append = TRUE)
    }
  }

  return(NULL)

}



#' Create download link within an RMarkdown log file
#' @param data data to be downloaded
#' @param logfile the main log file filepath
#' @param filename name of the csv file to be downloaded
#' @param linktext text on the download link
#' @param include.rownames should rownames of the dataframe be included in the CSV file?
#'
#' @export
create_download_link <- function(data, logfile, filename, linktext = 'Download the data', include.row.names = FALSE, verbose = F){
  if(verbose){
    csvfile <- file.path(dirname(logfile), paste0(filename))
    write.csv(data, csvfile, row.names = include.row.names)
    write(paste0("\n[", linktext , "](./", basename(csvfile), ")\n    \n<br>"), file = logfile, append = TRUE)
  }
}


#' Round numbers using the more conventional rounding method - halfway goes up
#'
#' @description
#'   This is a function that rounds numbers using the more conventional rounding method
#'
#' @details
#'   The R programming language rounds numbers using "Bankers Rounding" where if the number is halfway, it goes to the nearest even number.
#'   For example, 4.5 rounds to 4 rather than 5. This will lead to inconsistencies between the SQOUnified R package and the calculations of others.
#'   As of today, August 28, 2024 the only SQO Calculation method adopted by the State Waterboard is the Excel calculator created by Darrin Greenstein at SCCWRP.
#'   Excel rounds using the traditional arithmetic method, as do most people in the world.
#'   Therefore this function is added as a utility to this R package to ensure consistency of rounding methods.
#'   It is given the name "arithmetic round" because that is the name for this type of rounding, which is the most traditional, familiar method to most people.
#'   R's default rounding scheme is called "bankers rounding"
#'
#' @param x a number
#' @param digits the number of digits we will round to
#'
#' @export
arithmetic.round <- function(x, digits = 0) {
  adjust <- ifelse(x >= 0, 0.5, -0.5) * 10^(-digits)
  round(x + adjust, digits = digits)
}

