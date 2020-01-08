## Functions for computing differentially accessible sites/peaks. Adapted from Trapnell lab code


#' Find differentially accessible peaks/genes for each cluster
#'
#' @param counts Binary accessibility matrix (peaks x cells)
#' @param clusters Cluster assignments encoded as a factor
#' @param min.p Minimum p-value possible
#'
#' @return List of dataframes with cluster specificity, logfc, and p-value for each site tested
#' @export
#'
find_all_diff <- function(counts, clusters, min.p = 1e-100) {
  stopifnot(length(clusters) == ncol(counts))
  clusters <- factor(clusters[!is.na(clusters)])
  if (!is.null(names(clusters))) {
    counts <- counts[,names(clusters)]
  } else {
    counts <- counts[,!is.na(clusters)]
  }

  ## Calculate p-values with Fisher's exact test
  site.p <- calc_site_phyper(counts, clusters)
  site.p[site.p < min.p] <- min.p
  site.q <- matrix(p.adjust(as.numeric(site.p), method = "BH"), nrow = nrow(site.p), ncol = ncol(site.p))
  rownames(site.q) <- rownames(site.p); colnames(site.q) <- colnames(site.p);

  ## Calculate site specificity for each cluster
  site.prop <- calc_site_prop(counts, clusters)
  site.logfc <- calc_site_logfc(counts, clusters)

  cl.depth <- Matrix::colSums(counts)
  cl.depth <- as.numeric(tapply(cl.depth, clusters, median))
  names(cl.depth) <- levels(clusters)

  site.spec <- calc_site_specificity(site.prop, cl.depth)

  ## Add logfc and specificity to differential accessibility results
  diff_peaks_list <- lapply(levels(clusters), function(cl) {
    p <- site.p[,cl]
    q <- site.q[,cl]
    logfc <- site.logfc[,cl]
    specificity <- site.spec[, cl]
    df <- data.frame(pval = p, qval = q, logfc = logfc, specificity = specificity)
    rownames(df) <- rownames(counts)
    df
  })
  names(diff_peaks_list) <- levels(clusters)

  diff_peaks_list
}



#' Find top differentially accessible peaks for each cluster
#'
#' @param diff_peaks_list List of dA peaks per cluster output by find_all_diff
#' @param qval.cutoff Q-value cutoff to call a site differentially accessible
#' @param specificity.cutoff Specificity score cutoff
#' @param logfc.cutoff Log fold-change cutoff
#' @param max.peaks Max number of peaks to keep for each cluster
#'
#' @return Filtered list of granges with cluster specificity, logfc, and p-value for each site tested
#' @import GenomicRanges
#' @export
#'
top_diff_peaks <- function(diff_peaks_list, qval.cutoff = 0.05, specificity.cutoff = 1e-4,
                           logfc.cutoff = 0.5, max.peaks = 2e4) {
  lapply(diff_peaks_list, function(df) {
    df <- subset(df, qval < qval.cutoff & specificity > specificity.cutoff & logfc > logfc.cutoff)
    df <- df[order(df$specificity, decreasing = T),]
    df <- df[1:min(max.peaks, nrow(df)),]

    peaks.df <- peak2df(rownames(df), metadata.df = df[,c("qval", "logfc")], keep.colnames = T)
    GenomicRanges::makeGRangesFromDataFrame(peaks.df, keep.extra.columns = T)
  })
}



## Helper functions for differential accessibility

## Fisher's exact test
calc_site_phyper <- function(counts, clusters) {
  ## total number of site observations
  pickedN <- Matrix::rowSums(counts)

  ## number of site observations per cluster
  pickedRed <- t(as.matrix(sparse_agg_row(counts, clusters)))
  colnames(pickedRed) <- rownames(counts)
  rownames(pickedRed) <- levels(clusters)

  ## total accessible sites per cluster
  redN <- Matrix::rowSums(pickedRed)

  ## total accessible sites not in cluster
  whiteN <- sum(redN) - redN

  p <- t(do.call(rbind, lapply(1:nrow(pickedRed), function(i) {
    phyper(pickedRed[i,], redN[i], whiteN[i], pickedN + 1, lower.tail = F, log.p = F)
  })))
  rownames(p) <- rownames(counts)
  colnames(p) <- levels(clusters)

  return(p)
}


calc_site_prop <- function(counts, clusters) {
  stopifnot(all(colnames(counts) %in% names(clusters)))
  counts <- counts[,names(clusters)]

  acc.sites <- as.matrix(sparse_agg_row(counts, clusters))
  rownames(acc.sites) <- rownames(counts)
  colnames(acc.sites) <- levels(clusters)

  cluster.sizes <- as.numeric(table(clusters)[colnames(acc.sites)])
  acc.sites/cluster.sizes
}


calc_site_logfc <- function(counts, clusters) {
  stopifnot(all(colnames(counts) %in% names(clusters)))
  counts <- counts[,names(clusters)]

  acc.sites <- as.matrix(sparse_agg_row(counts, clusters))
  rownames(acc.sites) <- rownames(counts)
  colnames(acc.sites) <- levels(clusters)

  cluster.sizes <- table(clusters)[colnames(acc.sites)]
  sapply(levels(clusters), function(cl) {
    cl.prop <- acc.sites[,cl]/cluster.sizes[[cl]]

    if(nlevels(clusters) > 2) {
      ctrl.prop <- rowSums(acc.sites[,colnames(acc.sites) != cl])
      ctrl.prop <- ctrl.prop/sum(cluster.sizes[names(cluster.sizes) != cl])
    } else {
      ctrl.cl <- levels(clusters)[levels(clusters) != cl]
      ctrl.prop <- acc.sites[,ctrl.cl]/cluster.sizes[[ctrl.cl]]
    }
    eps <- 1/cluster.sizes[[cl]]
    log2((cl.prop + eps)/(ctrl.prop + eps))
  })
}



sparse_agg_row <- function(x, fac) {
  # cast into handy Matrix sparse Triplet form
  x.T <- as(x, "dgTMatrix")

  # factor column indexes (compensating for 0 vs 1 indexing)
  x.T@j <- as.integer(as.integer(fac[x.T@j + 1]) - 1)

  # cast back, magically computing factor sums along the way :)
  y <- as(x.T, "matrix.csr")

  # and fix the dimension (doing this on x.T bus errors!)
  y@dimension <- as.integer(c(nrow(y), nlevels(fac)))
  y
}


## Calculate site specificity
calc_site_specificity <- function(propmat, depth) {
  print("Normalizing proportions...")
  depthnorm = mean(depth)/depth
  logdepthnorm = mean(log10(depth))/log10(depth)
  propmatnormbylogdepth = t(t(propmat)*logdepthnorm)

  print("Calculating specificity scores...")
  marker_specificities_out = specificity_scorer(propmatnormbylogdepth)
  markerdup = marker_specificities_out^2
  markering = markerdup * propmatnormbylogdepth
  rownames(markering) = rownames(propmat)
  colnames(markering) = colnames(propmat)
  return(markering)
}


## Functions for calculating marker specificity
## (adapted from specificy_calculator.R from Cusanovich et al, 2018)
makeprobsvec <- function(p){
  phat<-p/sum(p)
  phat[is.na(phat)] = 0
  phat
}


shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
  -sum( log2(p.norm)*p.norm)
}


JSdistVec <- function(p, q){
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) +
                                           shannon.entropy(q)) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  JSdist
}


specificity_scorer <- function(normpropmat) {
  marker_gene_specificities <- lapply(1:ncol(normpropmat), function(cell_type_i){
    perfect_specificity <- rep(0.0, ncol(normpropmat))
    perfect_specificity[cell_type_i] <- 1.0
    apply(normpropmat, 1, function(x) {
      if (sum(x) > 0) 1 - JSdistVec(makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  return(do.call(cbind, marker_gene_specificities))
}

