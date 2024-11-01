\name{SoCalRivpacs}
\alias{SoCalRivpacs}
\alias{SFBayRivpacs}
\title{
RIVPACS for the Southern California Bight and San Francisco Bay
}
\description{
Perform the RIVPACS calculations using reference data for either the Southern California Bight (SoCal) or San Francisco Bay (SFBay).
}
\usage{
SoCalRivpacs(Pcutoff = 0.5,
reference.groups = socal.reference.groups,
observed.predictors = socal.example.habitat,
reference.taxa = socal.reference.taxa,
group.means = socal.reference.group.means,
reference.cov = socal.reference.covariance,
observed.taxa = socal.example.taxa)

SFBayRivpacs(Pcutoff = 0.5,
reference.groups = sfbay.reference.groups,
observed.predictors = sfbay.example.habitat,
reference.taxa = sfbay.reference.taxa,
group.means = sfbay.reference.group.means,
reference.cov = sfbay.reference.covariance,
observed.taxa = sfbay.example.taxa)
}
\arguments{
  \item{Pcutoff}{
Probability cutoff.
}
 \item{reference.groups}{
Reference groups.
}
 \item{observed.predictors}{
Example field habitat data supplied by the user.
}
 \item{reference.taxa}{
Reference taxa.
}
 \item{group.means}{
Group means.
}
 \item{reference.cov}{
Reference covariance matrix.
}
\item{observed.taxa}{
Example field taxa abundance supplied by the user. 
}
}
\value{
A list containing the observed/expected data frame, observed data, expected data, probability cutoff, and region. 
}
