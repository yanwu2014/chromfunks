% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coaccess.R
\name{FilterConns}
\alias{FilterConns}
\title{Reformat coaccessible peaks in list format}
\usage{
FilterConns(conns, min.coaccess)
}
\arguments{
\item{conns}{Dataframe of peak to peak connections}

\item{peaks}{to find coaccessible peaks for}
}
\value{
List of coaccessible peaks (in granges format) for each peak
}
\description{
Reformat coaccessible peaks in list format
}
