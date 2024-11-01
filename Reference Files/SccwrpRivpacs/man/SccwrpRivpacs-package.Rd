\name{SccwrpRivpacs-package}
\alias{SccwrpRivpacs-package}
\docType{package}

\title{RIVPACS}

\description{Calculation of River Invertebrate Prediction and Classification System (RIVPACS) index scores. This package is a sediment quality assessment tool intended to be used in the State of California's Sediment Quality Objectives (SQO). The RIVPACS calculations in this package are a modification of those obtained from http://www.epa.gov/wed/pages/models/rivpacs/rivpacs.htm

This package is used on a webserver and therefore has features intended for that use. For example, timestamps for file handling. An R script on the webserver takes the user's \code{.csv} files, splits the data in to Southern California Bight (SoCal) and San Francisco Bay (SF Bay) portions, and sends those to the package functions SoCalRivpacs() and SFBayRivpacs(). Those functions perform the RIVPACS calculations. The output is then sent to the function HtmlOutput() which produces 3 \code{.html} files for each of the 2 geographic regions. Those \code{.html} files are displayed by the webserver.

Data from one or both of the geographic regions may be supplied. Data from other regions is not valid input because reference data is avaiable for those two regions only.

}

\keyword{package}
