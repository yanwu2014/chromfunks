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
#' @param delim Peak name delimiter
#'
#' @return Peak info in dataframe format
#' @import GenomicRanges
#' @export
#'
peak2df <- function(peak.names, keep.colnames = F, metadata.df = NULL,
                    delim = c(":", "-")) {
  if (length(delim) > 1 && delim[[1]] == ":" && delim[[2]] == "-") {
    chrom <- sapply(peak.names, function(x) strsplit(x, split = ":")[[1]][[1]])
    bp.range <- sapply(peak.names, function(x) strsplit(x, split = ":")[[1]][[2]])
    bp1 <- as.integer(sapply(bp.range, function(x) strsplit(x, split = "-")[[1]][[1]]))
    bp2 <- as.integer(sapply(bp.range, function(x) strsplit(x, split = "-")[[1]][[2]]))
  } else if (delim == "-") {
    chrom <- sapply(peak.names, function(x) strsplit(x, split = "-")[[1]][[1]])
    bp1 <- as.integer(sapply(peak.names, function(x) strsplit(x, split = "-")[[1]][[2]]))
    bp2 <- as.integer(sapply(peak.names, function(x) strsplit(x, split = "-")[[1]][[3]]))
  } else {
    stop("Invalid delim")
  }


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
#' @param delim Peak name delimiters
#'
#' @return Peak info in granges format
#' @import GenomicRanges
#' @export
#'
peak2granges <- function(peak.names, metadata.df = NULL, delim = c(":", "-")) {
  bed.df <- peak2df(peak.names, keep.colnames = T, metadata.df = metadata.df,
                    delim = delim)
  makeGRangesFromDataFrame(bed.df, keep.extra.columns = T)
}


#' Function for converting a granges object into a character vector of peak names
#'
#' @param gr GenomicRanges object
#' @param delim Peak name delimiter
#'
#' @return Character vector of peak names
#' @import GenomicRanges
#' @export
#'
granges2peak <- function(gr, delim = c(":", "-")) {
  paste0(seqnames(gr), delim[[1]],
         start(gr), delim[[2]],
         end(gr))
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


#' Flatten Matrix into a vector
#'
#' @param mat Matrix
#' @param delim Delimiter to separate row and column names
#' @return Vector with row and column names concatenated by delim
#' @export
#'
flattenMat <- function(mat, delim = "_") {
  df <- reshape2::melt(mat)
  vec <- df$value
  names(vec) <- paste0(df$Var1, delim, df$Var2)
  vec
}

#' Unflatten a vector into a matrix
#'
#' @param v Vector
#' @param delim Delimiter to separate row and column names
#' @return Matrix
#' @export
#'
unflattenVec <- function(v, delim = "_") {
  x.lab <- sapply(names(v), ExtractField, delim = delim, field = 1)
  y.lab <- sapply(names(v), ExtractField, delim = delim, field = 2)
  df <- data.frame(x.lab, y.lab, v)
  UnflattenDataframe(df, output.name = "v", row.col = "x.lab", col.col = "y.lab")
}



#' Helper function for generating a pseudo-bulk matrix
#'
#' @param mat Input single cell matrix
#' @param celltype Cell type anotations
#'
#' @return A pseudo-bulk matrix where the columns have been summed by cell type
#'
#' @import Matrix
#' @export
#'
getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype == ct]
    pseudobulk <- Matrix::rowSums(mat[,cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}



#' Compute deviations on a ATAC counts matrix
#' @param counts Bulk or pseudo-bulk ATAC counts matrix
#'
#' @return A matrix of deviations from the expected counts
#' @export
#'
getBulkDev <- function(counts) {
  require(SummarizedExperiment)
  require(chromVAR)

  cluster.sums <- colSums(counts)
  peaks <- peak2granges(rownames(counts))
  SE <- SummarizedExperiment(assays = list(counts = counts),
                             rowData = peaks,
                             colData = data.frame(Name = colnames(counts)))
  peak.frac <- computeExpectations(SE)
  expected.counts <- sapply(cluster.sums, function(x) x*peak.frac)

  ## Compute deviation from expected counts
  (counts - expected.counts)/expected.counts
}



#' Get gene body information from a human txdb object
#'
#' @param txdb TxDb object (i.e. TxDb.Hsapiens.UCSC.hg38.knownGene)
#' @return GRanges of gene body annotations
#'
#' @export
#'
geneAnno <- function(txdb) {
  require(org.Hs.eg.db)
  require(annotate)

  genes.anno <- genes(txdb)
  names(genes.anno) <- genes.anno$gene_id
  genes.anno <- genes.anno[!is.na(genes.anno$gene_id)]

  entrez.to.symbol <- getSYMBOL(genes.anno$gene_id, data = "org.Hs.eg.db")
  genes.anno <- genes.anno[names(entrez.to.symbol)]
  genes.anno$gene_symbol <- entrez.to.symbol

  genes.anno
}


#' Extract the desired field from a delimited string
#'
#' @param string Input string
#' @param field Field position
#' @param delim Delimiter (i.e. tabs, spaces, underscores, commas)
#'
#' @return Extracted field
#'
extractField <- function (string, field = 1, delim = "_") {
  fields <- as.numeric(unlist(strsplit(x = as.character(x = field), split = ",")))
  if (length(fields) == 1) {
    return(strsplit(string, split = delim)[[1]][field])
  }
  return(paste(strsplit(string, split = delim)[[1]][fields], collapse = delim))
}
