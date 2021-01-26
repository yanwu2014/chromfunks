## Functions for handling chromatin accessibility matrices

#' Build accessibility matrices for peaks overlapping list of arbitrary regions
#'
#'@param peak.counts Chromatin accessibility counts matrix
#'@param peak.maps List of peaks
#'@param idents Cell identities (defaults to NULL)
#'@param scale Whether or not to scale the single cell matrices
#'
#'@return List of counts matrices
#'@export
#'
GetRegionCounts <- function(peak.counts, peak.maps, idents = NULL, scale = F) {
  lapply(peak.maps, function(peaks) {
    region.counts <- peak.counts[peaks,]
    if (scale) {
      region.counts <- t(apply(region.counts, 1, scale, center = T, scale = T))
    }
    if (!is.null(idents)) {
      idents <- factor(idents)
      region.counts <- do.call(cbind, lapply(levels(idents), function(cl) {
        ix <- names(idents[idents == cl])
        Matrix::rowMeans(region.counts[,ix])
      }))
      colnames(region.counts) <- levels(idents)
    }
    region.counts
  })
}


#' Get binned counts matrix
#' Adapted from the Cicero package (https://cole-trapnell-lab.github.io/cicero-release/)
#'
#' @param counts Counts matrix
#' @param reduced_coordinates Dimensionality reduction used to define bins
#' @param k Bin size
#' @param silent If true will print bin stats
#' @return Binned counts matrix
#' @export
#'
GetBinCounts <- function(counts, reduced_coordinates, k = 50, silent = T) {
  require(FNN)

  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
                   replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]

  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other, ], this_choice))
  }
  while (length(good_choices) > 0 & it < 5000) {
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), size = 1,
                     replace = FALSE)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen, ]
    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample), ]
    shared <- sapply(others, get_shared, this_choice = this_choice)
    if (max(shared) < 0.9 * k) {
      chosen <- new_chosen
    }
  }
  cell_sample <- nn_map[chosen, ]

  if (!silent) {
    combs <- combn(nrow(cell_sample), 2)
    shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
    })
    message(paste0("Overlap QC metrics:\nCells per bin: ",
                   k, "\nMaximum shared cells bin-bin: ", max(shared),
                   "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
                   median(shared)))
    if (mean(shared)/k > 0.1)
      warning("On average, more than 10% of cells are shared between paired bins.")
  }

  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(counts)) %in%
                   cell_sample[x, , drop = FALSE])
  mask <- as(mask, "dgCMatrix")
  bin.counts <- counts %*% mask
  colnames(bin.counts) <- paste0("agg", chosen)
  bin.counts
}



#' Compute the log TF-IDF transform for accessibility matrices
#' Adapted from Signac
#'
#' @param counts Accessibility counts matrix
#' @param scale.factor Factor to scale normalized counts with. Default is 10000.
#' @param verbose Print progress
#'
#' @return Counts matrix normalized with the log(TF-IDF + 1) transform
#' @export
#'
LogTFIDF <- function(counts, scale.factor = 1e4, verbose = T) {
  counts <- as(counts, Class = "dgCMatrix")

  if (verbose) {
    message("Performing TF-IDF normalization")
  }

  ## Term frequency
  npeaks <- colSums(counts)
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }
  tf <- tcrossprod(counts, y = Diagonal(x = 1 / npeaks))

  ## Inverse document frequency
  rsums <- rowSums(counts)
  if (any(rsums == 0)) {
    warning("Some features contain 0 total counts")
  }
  idf <- ncol(counts) / rsums

  ## Compute TF-IDF
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf

  ## Take log transform
  slot(norm.data, name = "x") <- log1p(
    x = slot(norm.data, name = "x") * scale.factor
  )
  colnames(norm.data) <- colnames(counts)
  rownames(norm.data) <- rownames(counts)

  ## Set NA values to 0
  vals <- slot(norm.data, name = "x")
  vals[is.na(x = vals)] <- 0
  slot(norm.data, name = "x") <- vals

  return(norm.data)
}


#' Average counts matrix by cluster/cell type
#'
#' @param counts Counts matrix
#' @param idents Cluster/Cell Type identities
#' @param binarize Binarize the matrix before computing averages
#'
#' @import Matrix
#' @return Counts matrix averaged over clusters/cell types
#' @export
#'
GetClusterAvg <- function(counts, idents, binarize = F) {
  idents <- factor(idents)
  stopifnot(all(colnames(counts) %in% names(idents)))

  counts <- counts[,names(idents)]
  counts <- as(counts, "dgCMatrix")
  if (binarize) {
    counts@x[counts@x > 0] <- 1
  }

  avg.counts <- do.call(cbind, lapply(levels(idents), function(cl) {
    cells <- names(idents[idents == cl])
    Matrix::rowMeans(counts[,cells])
  }))
  colnames(avg.counts) <- levels(idents)

  return(avg.counts)
}
