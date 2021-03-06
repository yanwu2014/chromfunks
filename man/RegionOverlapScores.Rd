% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region_enrich.R
\name{RegionOverlapScores}
\alias{RegionOverlapScores}
\title{Get overlaps between list of granges and granges in matrix format}
\usage{
RegionOverlapScores(ranges, regionsList)
}
\arguments{
\item{ranges}{GenomicRanges of dataset peaks (will be rows of output matrix)}

\item{regionsList}{List of regions in GenomicRanges format (will be columns of output matrix)}
}
\value{
Binary Matrix of ranges by regions overlaps (1 for overlap, 0 otherwise)
}
\description{
Get overlaps between list of granges and granges in matrix format
}
