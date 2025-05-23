% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{writelog}
\alias{writelog}
\title{Function to write to a log file}
\usage{
writelog(
  content,
  logfile,
  filetype = "text",
  append = T,
  verbose = F,
  prefix = NULL,
  include.row.names = F,
  code = NULL,
  data = NULL,
  echo.code = TRUE,
  pageLength = 10
)
}
\arguments{
\item{content}{Content to write to the logfile - text or a dataframe/tibble-like object}

\item{logfile}{Path to the file where the content will be written}

\item{filetype}{'txt' or 'csv' - the type of file}

\item{append}{Append to existing file? or overwrite existing contents}

\item{verbose}{This is what it will be called in parent functions - verbose tells whether to actually write to the logfile or not
This will allow us to write the code in the other functions and files without a bunch of if statements .}
}
\description{
Function to write to a log file
}
\examples{
writelog(df, 'logs/test.csv', filetype = 'csv', append = T)
writelog("# --- Test Log Statement --- #", 'logs/log.txt')
}
