## Functions for computing the enrichment of genomic regions (i.e. GWAS-linked SNPs)
## with cluster specific diff accessible sites


#' Calculate p-values linking phenotypes to clusters using permutation testing
#'
#' @param pheno.blocks A list of genomic regions for each GWAS phenotype in granges format
#' @param top_diff_peaks A list of top accessible sites for each cluster in granges format
#' @param background_peaks A granges of background accessible peaks to sample from
#' @param n.rand Number of permutations
#' @param n.cores Number of cores for parallel processing
#' @param score.col Whether or not to use z-scores or just the number of overlapping regions
#'
#' @return A matrix of empirical p-values for each phenotype/cluster enrichment test
#' @import snow
#' @import GenomicRanges
#' @import abind
#' @export
#'
cluster_pheno_pvals <- function(pheno.blocks, top_diff_peaks, background_peaks,
                                n.rand = 1e4, n.cores = 8, z.score = T) {
  scores <- cluster_pheno_scores(pheno.blocks, top_diff_peaks, z.score)

  n_diff_peaks <- sapply(top_diff_peaks, length)
  cl <- snow::makeCluster(n.cores, type = "SOCK")
  invisible(snow::parLapply(cl, 1:n.cores, function(i) library(GenomicRanges)))
  snow::clusterExport(cl, c("background_peaks", "cluster_pheno_scores", "n_diff_peaks",
                            "pheno.blocks", "z.score"),
                      envir = environment())
  rand.scores <- snow::parLapply(cl, 1:n.rand, function(i) {
    idx.draw <- 1:length(background_peaks)
    rand_diff_peaks <- lapply(n_diff_peaks, function(n) {
      idx.sample <- sample(idx.draw, size = n)
      background_peaks[idx.sample]
    })
    cluster_pheno_scores(pheno.blocks, rand_diff_peaks, z.score)
  })
  on.exit(snow::stopCluster(cl))
  rand.scores <- abind::abind(rand.scores, along = 3)

  t(calc_emp_pvals(scores, rand.scores))
}



#' Intersect phenotype-linked LD blocks with peaks
#'
#' @param pheno.blocks A list of genomic regions for each GWAS phenotype in granges format
#' @param top_diff_peaks A list of top accessible sites for each cluster in granges format
#' @param z.score Whether or not to use z-scores or just the number of overlapping regions
#'
#' @return A matrix of overlap scores for each phenotype/cluster combination
#' @import GenomicRanges
#' @export
#'
cluster_pheno_scores <- function(pheno.blocks, top_diff_peaks, z.score = T) {
  sapply(pheno.blocks, function(blocks) {
    sapply(top_diff_peaks, function(gr) {
      hits <- suppressWarnings(GenomicRanges::findOverlaps(blocks, gr, ignore.strand = T))
      blocks.ix <- unique(hits@from)
      if (z.score) {
        return(sum(blocks[blocks.ix]$z))
      } else {
        return(length(blocks.ix))
      }
    })
  })
}



#' Parse bed files
#' @param dir.path Path to directory with bed files
#'
#' @return List of genomic ranges, one for each file in the directory
#'
#' @import ChIPseeker
#' @import GenomicRanges
#' @export
#'
parse_bed_directory <- function(dir.path) {
  require(ChIPseeker)
  bed.files <- list.files(path = dir.path, pattern = ".bed")
  peaks_list <- lapply(bed.files, function(fi) {
    fi <- paste0(dir.path, "/", fi)
    gr <- readPeakFile(fi)
    names(gr) <- paste(seqnames(gr), start(gr), end(gr), sep = "_")
    gr
  })
  names(peaks_list) <- bed.files
  peaks_list
}


## Helper functions


## Calculate empirical p-values
calc_emp_pvals <- function(cfs, cfs.rand) {
  p.vals <- matrix(1, nrow(cfs), ncol(cfs))
  colnames(p.vals) <- colnames(cfs)
  rownames(p.vals) <- rownames(cfs)
  for (i in 1:nrow(cfs)) {
    for (j in 1:ncol(cfs)) {
      v <- sort(na.omit(cfs.rand[i,j,]))
      b <- cfs[i,j]
      p.vals[i,j] <- (sum(v >= b) + 1)/(length(v) + 1)
    }
  }
  return(p.vals)
}


