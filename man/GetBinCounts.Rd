% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chrom.R
\name{GetBinCounts}
\alias{GetBinCounts}
\title{Get binned counts matrix
Adapted from the Cicero package (https://cole-trapnell-lab.github.io/cicero-release/)}
\usage{
GetBinCounts(counts, reduced_coordinates, k = 50, silent = T)
}
\arguments{
\item{counts}{Counts matrix}

\item{reduced_coordinates}{Dimensionality reduction used to define bins}

\item{k}{Bin size}

\item{silent}{If true will print bin stats}
}
\value{
Binned counts matrix
}
\description{
Get binned counts matrix
Adapted from the Cicero package (https://cole-trapnell-lab.github.io/cicero-release/)
}
