## Functions for handling coaccessibility (defined using either Cicero or Hi-C)


#' Reformat coaccessible peaks in list format
#'
#' @param peaks to find coaccessible peaks for
#' @param conns Dataframe of peak to peak connections
#'
#' @return List of coaccessible peaks (in granges format) for each peak
#' @export
#'
get_coaccess_peaks <- function(peaks, conns) {
  conns.df <- subset(conns, Peak1 %in% peaks)
  peaks <- peaks[peaks %in% conns.df$Peak1]
  peaks.conns.list <- lapply(peaks, function(peak) {
    df <- subset(conns.df, Peak1 == peak)
    gr <- peak2granges(as.character(df$Peak2))
    gr$coaccess <- df$coaccess
    gr
  })
  names(peaks.conns.list) <- peaks

  peaks.conns.list
}



#' Reformat coaccessible peaks in list format
#'
#' @param peaks to find coaccessible peaks for
#' @param conns Dataframe of peak to peak connections
#'
#' @return List of coaccessible peaks (in granges format) for each peak
#' @export
#'
filter_conns <- function(conns, min.coaccess) {
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
region_overlap <- function(gr1, gr2, gr1.name = "HAR") {
  overlap.ix <- findOverlaps(gr1, gr2)

  overlap.ix.vector <- as.character(overlap.ix@to)
  names(overlap.ix.vector) <- as.character(overlap.ix@from)

  unique.groups <- unique(overlap.ix.vector)
  overlap.ix.list <- vector(mode = "list", length = length(unique.groups))
  names(overlap.ix.list) <- unique.groups

  for (g in unique.groups) {
    g.cells <- names(overlap.ix.vector[overlap.ix.vector == g])
    overlap.ix.list[[g]] <- g.cells
  }
  names(overlap.ix.list) <- gr1@elementMetadata[[gr1.name]][as.integer(unique.groups)]

  lapply(overlap.ix.list, function(peaks.ix) {
    gr2[as.integer(peaks.ix)]
  })
}



#' Compute relationship between peaks and genes using coaccessibility/promoters
#'
#' @param peaks.gr Peaks in granges format
#' @param conns Dataframe of peak to peak coaccessibility
#' @param link.promoter Include peaks in gene promoters
#' @param promoter.region Specify the window around the TSS that counts as a promoter
#' @param anno.level Specify "gene" or "transcript" for a gene/transcript level annotation
#'
#' @return Matrix of peak to gene coaccessibilities. A score of 1 means the peak is in the gene promoter.
#'
#' @import GenomicRanges
#' @import ChIPseeker
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @export
#'
peak_gene_coaccess <- function(peaks.gr, conns, link.promoter = T,
                               promoter.region = c(-5000, 5000),
                               anno.level = "gene") {
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)

  names(peaks.gr) <- paste0(seqnames(peaks.gr), ":", start(peaks.gr), "-", end(peaks.gr))
  peaks.gr.anno <- annotatePeak(peaks.gr, tssRegion = promoter.region,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                level = anno.level,
                                annoDb = "org.Hs.eg.db")
  peaks.gr.anno <- as.data.frame(peaks.gr.anno)
  rownames(peaks.gr.anno) <- names(peaks.gr)
  peaks.gr.anno <- subset(peaks.gr.anno, grepl("Promoter", annotation))

  peak2gene <- peaks.gr.anno$SYMBOL
  names(peak2gene) <- rownames(peaks.gr.anno)
  peak2gene <- peak2gene[!is.na(peak2gene)]

  gene.peak.promoter.list <- list()
  unique.prom.genes <- unique(peak2gene)
  for (gene in unique.prom.genes) {
    gene.peak.promoter.list[[gene]] <- names(peak2gene[peak2gene == gene])
  }

  conns <- subset(conns, Peak1 %in% names(peaks.gr))
  conns.peak2.gr <- peak2granges(conns$Peak2)
  conns.peak2.gr.anno <- annotatePeak(conns.peak2.gr, tssRegion = promoter.region,
                                      TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                      annoDb = "org.Hs.eg.db")
  conns.peak2.gr.anno <- as.data.frame(conns.peak2.gr.anno)
  conns.peak2.gr.anno$peak1 <- conns$Peak1
  conns.peak2.gr.anno$coaccess <- conns$coaccess
  conns.peak2.gr.anno <- subset(conns.peak2.gr.anno, grepl("Promoter", annotation))

  peak.gene.coaccess.list <- lapply(unique(conns.peak2.gr.anno$peak1), function(p) {
    ix <- which(conns.peak2.gr.anno$peak1 == p)
    w <- conns.peak2.gr.anno$coaccess[ix]; names(w) <- conns.peak2.gr.anno$SYMBOL[ix];
    tapply(w, names(w), max)
  })
  names(peak.gene.coaccess.list) <- unique(unique(conns.peak2.gr.anno$peak1))

  all.genes <- unique(union(peaks.gr.anno$SYMBOL, conns.peak2.gr.anno$SYMBOL))
  peak.gene.weights <- matrix(0, nrow = length(peaks.gr), ncol = length(all.genes))
  rownames(peak.gene.weights) <- names(peaks.gr)
  colnames(peak.gene.weights) <- all.genes

  ## Link coaccessibility
  for (p in names(peak.gene.coaccess.list)) {
    w <- peak.gene.coaccess.list[[p]]
    peak.gene.weights[p, names(w)] <- w
  }

  ## Link HAR to a gene if the HAR is in the promoter region
  if (link.promoter) {
    for(g in names(gene.peak.promoter.list)) {
      n.peaks <- length(gene.peak.promoter.list[[g]])
      peak.gene.weights[gene.peak.promoter.list[[g]], g] <- rep(1, n.peaks)
    }
  }

  peak.gene.weights
}




#' Compute relationship between arbitrary and genes using coaccessibility between peaks/promoters
#'
#' @param regions Genomic regions in granges format
#' @param regions.anno Annnotated genomic regions
#' @param conns Dataframe of peak to peak coaccessibility
#' @param hg38.chr.lengths Chromosome lengths
#' @param link.promoter Include peaks in gene promoters
#' @param promoter.region Specify the window around the TSS that counts as a promoter
#' @param anno.level Specify "gene" or "transcript" for a gene/transcript level annotation
#' @param buffer.size Buffer size around each region
#' @param region.name Name of region identifier
#'
#' @return Matrix of region to gene coaccessibilities. A score of 1 means the region is in the gene promoter.
#'
#' @import GenomicRanges
#' @import ChIPseeker
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @export
#'
region_gene_coaccess <- function(regions, regions.anno, conns, hg38.chr.lengths, link.promoter = F,
                                 promoter.region = c(-5000, 5000), anno.level = "gene",
                                 buffer.size = 1e3, region.name = "HAR") {
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)

  ## Link HAR to more distal genes using coaccessibility
  rownames(hg38.chr.lengths) <- hg38.chr.lengths[[1]]
  regions.buffered <- extend_gr(regions, buffer.size, hg38.chr.lengths)

  conns.from.gr <- peak2granges(conns$Peak1)
  regions.conns.ix <- findOverlaps(conns.from.gr, regions.buffered)

  conns.to.gr <- peak2granges(conns$Peak2)
  conns.to.gr.anno <- annotatePeak(conns.to.gr, tssRegion = promoter.region,
                                   level = anno.level,
                                   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                   annoDb = "org.Hs.eg.db")
  conns.to.gr.anno <- conns.to.gr.anno@anno
  conns.to.gr.anno$coaccess <- conns$coaccess

  regions.conns.ix.vector <- as.character(regions.conns.ix@to)
  names(regions.conns.ix.vector) <- as.character(regions.conns.ix@from)
  regions.conns.ix.list <- UnflattenGroups(regions.conns.ix.vector)
  names(regions.conns.ix.list) <- regions@elementMetadata[[region.name]][as.integer(names(regions.conns.ix.list))]

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
    names(region.gene.weights.list) <- regions@elementMetadata[[region.name]]
  } else {
    region.gene.weights.list <- lapply(1:length(regions), function(i) { c() })
    names(region.gene.weights.list) <- regions@elementMetadata[[region.name]]
  }

  ## Link coaccessibility
  for (h in names(regions.conns.ix.list)) {
    peaks.ix <- as.integer(regions.conns.ix.list[[h]])
    peaks.anno.gr <- conns.to.gr.anno[peaks.ix]
    peaks.anno.gr <- peaks.anno.gr[grepl("Promoter", peaks.anno.gr$annotation)]
    if (length(peaks.anno.gr) > 0) {
      w <- peaks.anno.gr$coaccess
      names(w) <- peaks.anno.gr$SYMBOL
      region.gene.weights.list[[h]] <- c(region.gene.weights.list[[h]], w)
    }
  }

  region.gene.weights.list <- lapply(region.gene.weights.list, function(w) {
    if (length(w) > 1) w <- tapply(w, names(w), max)
    w
  })

  all.region.genes <- unique(unlist(lapply(region.gene.weights.list, names), F, F))
  region.gene.weights <- matrix(0, length(regions), length(all.region.genes))
  rownames(region.gene.weights) <- regions@elementMetadata[[region.name]]
  colnames(region.gene.weights) <- all.region.genes
  for(h in names(region.gene.weights.list)) {
    region.gene.weights[h, names(region.gene.weights.list[[h]])] <- region.gene.weights.list[[h]]
  }

  region.gene.weights
}
