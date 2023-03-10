% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{annograph}
\alias{annograph}
\title{annograph: matrisome annotations with nested pie charts}
\usage{
annograph(result.obj = NULL, select = c("sender", "receiver"))
}
\arguments{
\item{result.obj}{character. The results from a previous call to matricom.}

\item{select}{character. Annotate matrisome families and categories for the sender or the receiver population. Default is "sender".}
}
\value{
a PieDonut (a nested pie chart) made with ggplot2
}
\description{
annograph: matrisome annotations with nested pie charts
}
\examples{
annograph(results)
}
