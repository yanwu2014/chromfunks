## Utility functions


#' Extend granges by buffer length
#'
#' @param gr GRanges object
#' @param buffer bp to extend gr by
#' @param chr.lengths chromosome lengths
#'
#' @return GRanges object with each range extended by the buffer
#' @import GenomicRanges
#' @export
#'
extend_gr <- function(gr, buffer, chr.lengths) {
  gr.names <- names(gr)
  seqlengths(gr) <- chr.lengths[names(seqlengths(gr)),2]
  new_start <- start(gr) - buffer
  new_end <- end(gr) + buffer
  ranges(gr) <- IRanges(new_start, new_end)
  gr <- trim(gr)
  names(gr) <- gr.names
  gr
}



#' Convert gRanges to dataframe
#'
#' @param gr GRanges object
#'
#' @return GRanges info in dataframe format
#' @import GenomicRanges
#' @export
#'
granges2df <- function(gr) {
  df <- data.frame(seqnames = seqnames(gr), starts = start(gr) - 1, ends = end(gr),
                   annot = gr$V4, names = gr$V6, genes = gr$V5)
  df
}




#' Function for converting peak names into a dataframe
#'
#' @param peak.names Peak names
#' @param keep.colnames Keep column names in resulting dataframe
#' @param metadata.df Extra metadata to annotate the peaks with
#'
#' @return Peak info in dataframe format
#' @import GenomicRanges
#' @export
#'
peak2df <- function(peak.names, keep.colnames = F, metadata.df = NULL) {
  chrom <- sapply(peak.names, function(x) strsplit(x, split = ":")[[1]][[1]])
  bp.range <- sapply(peak.names, function(x) strsplit(x, split = ":")[[1]][[2]])
  bp1 <- as.integer(sapply(bp.range, function(x) strsplit(x, split = "-")[[1]][[1]]))
  bp2 <- as.integer(sapply(bp.range, function(x) strsplit(x, split = "-")[[1]][[2]]))

  bed.df <- data.frame("chr" = chrom, "start" = bp1, "end" = bp2)
  if (!is.null(metadata.df)) bed.df <- cbind(bed.df, metadata.df)

  if (!keep.colnames) {
    colnames(bed.df) <- NULL
  }
  bed.df[complete.cases(bed.df),]
}



#' Function for converting peak names into a granges object
#'
#' @param peak.names Peak names
#' @param metadata.df Extra metadata to annotate the peaks with
#'
#' @return Peak info in granges format
#' @import GenomicRanges
#' @export
#'
peak2granges <- function(peak.names, metadata.df = NULL) {
  bed.df <- peak2df(peak.names, keep.colnames = T, metadata.df = metadata.df)
  makeGRangesFromDataFrame(bed.df, keep.extra.columns = T)
}



#' Inverts a list
#'
#' @param my.list List
#'
#' @return Inverted list
#' @export
#'
invertList <- function(my.list) {
  x <- unlist(my.list, F, T)
  a <- names(x); names(a) <- x;
  split(unname(a),names(a))
}



#' Inverts a list of named numeric vectors
#'
#' @param my.list List
#'
#' @return Inverted list
#' @export
#'
invertListVector <- function(my.list) {
  all.vec.names <- unique(unlist(lapply(my.list, names), F, F))
  list.names.to.vec.names <- invertList(lapply(my.list, names))
  inv.vec.list <- lapply(names(list.names.to.vec.names), function(l) {
    vec.names <- list.names.to.vec.names[[l]]
    w <- sapply(vec.names, function(h) my.list[[h]][[l]])
    names(w) <- vec.names
    w
  })
  names(inv.vec.list) <- names(list.names.to.vec.names)
}


#' Capitalize the first letter of every word in a string
#'
#' @param x String
#'
#' @return String where the first letter of every word is capitalized
#' @export
#'
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


#' Reformat rownames
#'
#' @param peak.names Peak names where the chromosomes and base pairs are separated by underscores (chr_start_end)
#'
#' @return Peak names formatted as so: chr:start-end
#' @export
#'
formatCisTopicPeaks <- function(peak.names) {
  chrom <- sapply(peak.names, function(x) strsplit(x, split = "_")[[1]][[1]])
  loc_start <- sapply(peak.names, function(x) strsplit(x, split = "_")[[1]][[2]])
  loc_end <- sapply(peak.names, function(x) strsplit(x, split = "_")[[1]][[3]])
  paste0(chrom, ":", loc_start,"-",loc_end)
}


#' Reformat rownames
#'
#' @param m Dense matrix
#'
#' @return The same matrix in dataframe format
#' @export
#'
FlattenMatrix <- function(m) {
  rows = dim(m)[1]
  cols = dim(m)[2]
  cbind(rowInd = rep(1:rows, times = cols),
        colInd = rep(1:cols, each = rows),
        reshape2::melt(m))
}



#' Helper function for generating a pseudo-bulk matrix
#'
#' @param mat Input single cell matrix
#' @param celltype Cell type anotations
#'
#' @return A pseudo-bulk matrix where the columns have been summed by cell type
#' @export
#'
getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype == ct]
    pseudobulk <- rowSums(mat[,cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}
