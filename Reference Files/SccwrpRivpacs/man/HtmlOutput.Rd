\name{HtmlOutput}
\alias{HtmlOutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write RIVPACS Output to HTML Files
%%  ~~function to do ... ~~
}
\description{
Takes RIVPACS results and writes 3 HTML files intended for the webserver.
}
\usage{
HtmlOutput(rivpacs, timestamp, user.filename, path = "/var/www/sqo/files/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rivpacs}{
The output from either \code{SoCalRivpacs} or \code{SFBayRivpacs}.
}
\item{timestamp}{
The timestamp from the Rscript command. Used to name the output files.
}
\item{user.filename}{
The user's filename from the Rscript command. Printed on the output files.
}
\item{path}{
The directory path for the output files. The default is that used on the server.
}
}
\value{
Three HTML files containing the observed/expected table, probability comparison tables, and probability matrix.
}
