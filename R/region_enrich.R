## Functions for computing enrichment in arbitrary sets of genomic regions


#' Get overlaps between list of granges and granges in matrix format
#'@import Matrix
#'
#'@param ranges GenomicRanges of dataset peaks (will be rows of output matrix)
#'@param regionsList List of regions in GenomicRanges format (will be columns of output matrix)
#'
#'@return Binary Matrix of ranges by regions overlaps (1 for overlap, 0 otherwise)
#'@export
#'
RegionOverlapScores <- function(ranges, regionsList) {
  la <- sapply(regionsList, function(hitG){
    # Get overlaps / summarized statistic summarized by function
    hitG$score <- 1
    ov <- findOverlaps(ranges, hitG)
    score <- rep(0, length(ranges))
    ss <- tapply(hitG$score[subjectHits(ov)], queryHits(ov), FUN = sum)
    score[as.integer(names(ss))] <- unname(ss)
    score
  })

  # Make a new SummarizedExperiment and export
  colnames(la) <- names(regionsList)
  as(la, "dgCMatrix")
}


#' Compute enrichment in each column of the chromatin accessibility matrix for all regions in regionsList
#'@import chromVAR
#'@import SummarizedExperiment
#'
#'@param SE Peaks accessibility matrix in SummarizedExperiment format
#'@param regionsList List of regions in GenomicRanges format
#'@param n.runs Number of replicates to run
#'
#'@return Matrix of z-scores representing the enrichment of each set of regions in each column of SE
#'@export
#'
ChromRegionEnrich <- function(SE, regionsList, n.runs = 5) {
  scores.mat <- RegionOverlapScores(rowRanges(SE), regionsList)

  dev_list <- lapply(1:n.runs, function(i) {
    # wDEV <- computeWeightedDeviations(SE, scores.mat)
    wDEV <- computeDeviations(SE, scores.mat)
    assays(wDEV)[["z"]]
  })

  apply(abind::abind(dev_list, along = 3), c(1,2), function(x) {
    x <- x[!(is.na(x) | is.infinite(x) | is.nan(x))]
    mean(x)
  })
}


#' Run chromVAR on a subset of peaks overlapping a list of regions
#'@import chromVAR
#'@import SummarizedExperiment
#'
#'@param SE Peaks accessibility matrix in SummarizedExperiment format
#'@param motif_ix Transcription factor motif data in SummarizedExperiment format (from motifmatchr)
#'@param regionsList List of regions in GenomicRanges format
#'@param n.runs Number of replicates to run
#'
#'@return For each set of regions z-scores representing the enrichment of TF motifs in peaks overlapping those regions for each column in SE
#'@export
#'
RegionChromVAR <- function(SE, motif_ix, regionList, n.runs) {
  region.peak.scores <- RegionOverlapScores(rowRanges(SE), regionList)
  peak.maps <- lapply(colnames(region.peak.scores), function(i) {
    ix <- region.peak.scores[,i]
    names(ix[ix > 0])
  })
  names(peak.maps) <- colnames(region.peak.scores)

  lapply(peak.maps, function(peaks) {
    dev_list <- lapply(1:n.runs, function(i) {
      dev <- computeDeviations(object = SE[peaks,], annotations = motif_ix[peaks,])
      assays(dev)[["z"]]
    })

    apply(abind::abind(dev_list, along = 3), c(1,2), function(x) {
      x <- x[!(is.na(x) | is.infinite(x) | is.nan(x))]
      mean(x)
    })
  })
}
