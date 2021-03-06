% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigSAME.R
\name{sigSAME}
\alias{sigSAME}
\title{Compute the likelihood for the observed somatic mutation and the read-depth data for gSAME.
It can significantly reduce computational time for assication analysis
between one gene-level mutation and lots of outcomes.}
\usage{
sigSAME(O, D, A, d0 = 20, min_altcount = 1, spesen = NULL,
  mix_4bb = NULL, mix_2bb = NULL)
}
\arguments{
\item{O}{A matrix for the observed somatic mutation.}

\item{D}{A matrix for the total read-depth.}

\item{A}{A matrix for the number of alternative number matrix.}

\item{d0}{The minimum of the total read-depth for obtaining the observed somatic mutation value. The default value is 20.}

\item{min_altcount}{The mimimum of the number of alternative reads that the somatic mutation could acutally orrur. The default value is 1.}

\item{spesen}{A dataframe specifying the specificity and sensitivity for all the somatic mutations.}

\item{mix_4bb}{A dataframe indicating the parameters of four beta-binomial distributions depending on the values of the observed somatic mutaton and the true somatic mutation when the read-depth is high.}

\item{mix_2bb}{A dataframe indicating the parameters of two beta-binomial distributions depending on the true value of the somatic mutation when the read-depth is low.}
}
\value{
A matrix with the likelihood of O, D, A conditioning on the true value of the somatic mutation (S=0/1) for all the samples.
}
\description{
Compute the likelihood for the observed somatic mutation and the read-depth data for gSAME.
It can significantly reduce computational time for assication analysis
between one gene-level mutation and lots of outcomes.
}
