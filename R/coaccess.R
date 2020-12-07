## Functions for handling coaccessibility (defined using either Cicero or Hi-C)


#' #' Reformat coaccessible peaks in list format
#' #'
#' #' @param peaks to find coaccessible peaks for
#' #' @param conns Dataframe of peak to peak connections
#' #'
#' #' @return List of coaccessible peaks (in granges format) for each peak
#' #' @export
#' #'
#' get_coaccess_peaks <- function(peaks, conns) {
#'   conns.df <- subset(conns, Peak1 %in% peaks)
#'   peaks <- peaks[peaks %in% conns.df$Peak1]
#'   peaks.conns.list <- lapply(peaks, function(peak) {
#'     df <- subset(conns.df, Peak1 == peak)
#'     gr <- peak2granges(as.character(df$Peak2))
#'     gr$coaccess <- df$coaccess
#'     gr
#'   })
#'   names(peaks.conns.list) <- peaks
#'
#'   peaks.conns.list
#' }



#' Reformat coaccessible peaks in list format
#'
#' @param peaks to find coaccessible peaks for
#' @param conns Dataframe of peak to peak connections
#'
#' @return List of coaccessible peaks (in granges format) for each peak
#' @export
#'
FilterConns <- function(conns, min.coaccess) {
  conns$coaccess <- abs(conns$coaccess)
  conns <- subset(conns, coaccess > min.coaccess)
  conns$Peak1 <- as.character(conns$Peak1)
  conns$Peak2 <- as.character(conns$Peak2)
  conns
}



#' Overlap genomic regions
#'
#' @param gr1 First set of regions to overlap (granges format)
#' @param gr2 Second set of regions to overlap (granges format)
#' @param gr1.name Column with the region names in gr1
#'
#' @return List of overlapping regions
#' @import GenomicRanges
#' @export
#'
RegionOverlapList <- function(gr1, gr2, gr1.name = "HAR") {
  if (is.null(gr1.name)) {
    names(gr1) <- granges2peak(gr)
  } else {
    names(gr1) <- gr1@elementMetadata[[gr1.name]]
  }

  overlap.ix <- findOverlaps(gr1, gr2)
  overlap.ix.vector <- as.character(overlap.ix@from)
  names(overlap.ix.vector) <- as.character(overlap.ix@to)

  unique.gr1 <- unique(overlap.ix.vector)
  overlap.ix.list <- vector(mode = "list", length = length(unique.gr1))
  names(overlap.ix.list) <- unique.gr1

  for (g in unique.gr1) {
    gr2.ix <- names(overlap.ix.vector[overlap.ix.vector == g])
    overlap.ix.list[[g]] <- gr2.ix
  }
  names(overlap.ix.list) <- names(gr1)[as.integer(unique.gr1)]

  lapply(overlap.ix.list, function(gr2.ix) {
    gr2[as.integer(gr2.ix)]
  })
}




#' Compute relationship between arbitrary regions and genes using either Hi-C or coaccessibility
#'
#' @param regions Genomic regions in granges format
#' @param conns Dataframe of peak to peak Hi-C/coaccessibility
#' @param link.promoter Include peaks in gene promoters
#' @param promoter.region Specify the window around the TSS that counts as a promoter
#' @param anno.level Specify "gene" or "transcript" for a gene/transcript level annotation
#' @param region.name Name of region identifier
#'
#' @return Matrix of region to gene Hi-C/coaccessibility connections
#'
#' @import GenomicRanges
#' @import ChIPseeker
#' @import Matrix
#' @export
#'
RegionGeneContact <- function(regions, conns, link.promoter = F, promoter.region = c(-3000, 3000),
                              anno.level = "transcript", region.name = "HAR") {
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)

  regions.peaks <- granges2peak(regions)
  if (is.null(region.name)) {
    names(regions) <- regions.peaks
  } else {
    names(regions) <- regions@elementMetadata[[region.name]]
  }

  regions.anno <- annotatePeak(regions, tssRegion = promoter.region, level = anno.level,
                               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                               annoDb = "org.Hs.eg.db")

  peak1.gr <- peak2granges(conns$Peak1)
  regions.conns.ix <- findOverlaps(peak1.gr, regions)

  peak2.gr <- peak2granges(conns$Peak2)
  peak2.gr.anno <- annotatePeak(peak2.gr, tssRegion = promoter.region, level = anno.level,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                annoDb = "org.Hs.eg.db")
  peak2.gr.anno <- peak2.gr.anno@anno
  peak2.gr.anno$coaccess <- conns$coaccess

  regions.conns.ix.vector <- as.character(regions.conns.ix@to)
  names(regions.conns.ix.vector) <- as.character(regions.conns.ix@from)
  regions.conns.ix.list <- lapply(unique(regions.conns.ix.vector), function(x) {
    names(regions.conns.ix.vector[regions.conns.ix.vector == x])
  })
  names(regions.conns.ix.list) <- names(regions)[as.integer(unique(regions.conns.ix.vector))]

  ## Link HAR to a gene if the HAR is in the promoter region
  if (link.promoter) {
    is.promoter <- grepl("Promoter", regions.anno$annotation)
    region.gene.weights.list <- lapply(1:length(regions), function(i) {
      if (is.promoter[[i]]){
        w <- 1
        names(w) <- regions.anno[["SYMBOL"]][[i]]
      } else {
        w <- c()
      }
      w
    })
    names(region.gene.weights.list) <- names(regions)
  } else {
    region.gene.weights.list <- lapply(1:length(regions), function(i) { c() })
    names(region.gene.weights.list) <- names(regions)
  }

  ## Link coaccessibility
  for (h in names(regions.conns.ix.list)) {
    peaks.ix <- as.integer(regions.conns.ix.list[[h]])
    peaks.anno.gr <- peak2.gr.anno[peaks.ix]
    peaks.anno.gr <- peaks.anno.gr[grepl("Promoter", peaks.anno.gr$annotation)]
    if (length(peaks.anno.gr) > 0) {
      w <- peaks.anno.gr$coaccess
      names(w) <- peaks.anno.gr$SYMBOL
      region.gene.weights.list[[h]] <- c(region.gene.weights.list[[h]], w)
    }
  }

  region.gene.weights.list <- lapply(region.gene.weights.list, function(w) {
    if (length(w) > 1) w <- tapply(w, names(w), max)
    w[!(is.na(w) | is.na(names(w)))]
  })
  region.gene.weights.list <- region.gene.weights.list[sapply(region.gene.weights.list, function(x) !is.null(x))]
  region.gene.weights.list <- region.gene.weights.list[sapply(region.gene.weights.list, function(x) length(x) > 0)]

  all.region.genes <- unique(unlist(lapply(region.gene.weights.list, names), F, F))
  region.gene.weights <- matrix(0, length(regions), length(all.region.genes))
  rownames(region.gene.weights) <- names(regions)
  colnames(region.gene.weights) <- all.region.genes
  for(h in names(region.gene.weights.list)) {
    region.gene.weights[h, names(region.gene.weights.list[[h]])] <- region.gene.weights.list[[h]]
  }

  as(region.gene.weights, "dgCMatrix")
}
