% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.pred}
\alias{plot.pred}
\title{Plot a `pred` object constructed for a categorical response}
\usage{
\method{plot}{pred}(x, ...)
}
\arguments{
\item{x}{pred object- a list classified as pred containing objects data and bound}

\item{...}{additional parameters passed to the default plot method}
}
\value{
capability to plot pred object. More details: the command `plot(obj)` 
plots the empirical densities of each category. Mass denoted in red indicates
inclusion in the prediction set
}
\description{
Plot a `pred` object constructed for a categorical response
}
\concept{pred plots}
