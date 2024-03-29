% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesNormalPrediction.R
\name{bayesNormalPrediction}
\alias{bayesNormalPrediction}
\title{Obtain a Bayesian prediction interval}
\usage{
bayesNormalPrediction(Y, alpha = 0.15, mu = 0, tau2 = 1)
}
\arguments{
\item{Y}{Observed data vector}

\item{alpha}{Prediction error rate}

\item{mu}{Prior expected mean of the population mean}

\item{tau2}{Prior expected variance of the population mean}
}
\value{
pred object
}
\description{
This function computes a Bayesian prediction interval based on a normal model.
}
