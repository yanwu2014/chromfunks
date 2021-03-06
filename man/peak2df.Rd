% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{peak2df}
\alias{peak2df}
\title{Function for converting peak names into a dataframe}
\usage{
peak2df(peak.names, keep.colnames = F, metadata.df = NULL, delim = c(":", "-"))
}
\arguments{
\item{peak.names}{Peak names}

\item{keep.colnames}{Keep column names in resulting dataframe}

\item{metadata.df}{Extra metadata to annotate the peaks with}

\item{delim}{Peak name delimiter}
}
\value{
Peak info in dataframe format
}
\description{
Function for converting peak names into a dataframe
}
