% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toxicity.R
\name{tox.sqo}
\alias{tox.sqo}
\title{Get Tox SQO Scores and Categories}
\usage{
tox.sqo(toxresults)
}
\arguments{
\item{toxresults}{a dataframe with the following columns: stationid, toxbatch, species, sampletypecode
   matrix, labrep, result. This data must also include the control samples
   (stationcode 0000, sampletypecode CNEG etc)

   The input dataframe is structured as follows:

   \strong{\code{stationid}} - an alpha-numeric identifier of the location;

   \strong{\code{lab}} -  (optional) The laboratory which performed the test

   \strong{\code{toxbatch}} - the toxbatch id - used to join with the control sample

   \strong{\code{species}} - The Genus and species of the animale that was tested

   \strong{\code{sampletypecode}} - The sampletype used Grab, CNEG etc. Control samples must be included

   \strong{\code{matrix}} - (optional) Whole Sediment, Sediment Water Interface, etc. Be sure to not include Reference Toxicants

   \strong{\code{labrep}} - There should be 5 per station, species pair

   \strong{\code{result}} - the percentage that survived the test, or had normal development}
}
\description{
This function will calculate the tox SQO scores for stations given a dataframe structured as described
in the details or the Arguments section. The funtion will get SQO scores for each individual test
conducted for a station, as well as the integrated Toxicity LOE SQO score and category
}
\examples{
data(tox_sampledata)
tox.sqo(tox_sampledata)

}
