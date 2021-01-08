## Functions for handling coaccessibility (defined using either Cicero or Hi-C)


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
RegionOverlapList <- function(gr1, gr2, gr1.name = NULL) {
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
#' @param region.name Column name of region identifier. Defaults to peak names
#' @param weight.col Column name of weights (i.e. when using coaccessibility). Defaults to 1
#'
#' @return Matrix of region to gene Hi-C/coaccessibility connections
#'
#' @import GenomicRanges
#' @import ChIPseeker
#' @import Matrix
#' @export
#'
RegionGeneLinks <- function(regions,
                            conns,
                            link.promoter = T,
                            promoter.region = c(-3000, 3000),
                            anno.level = "transcript",
                            region.name = NULL,
                            weight.col = NULL)
{
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)

  ## Get region names
  regions <- unique(regions)
  regions.peaks <- granges2peak(regions)
  if (is.null(region.name)) {
    names(regions) <- regions.peaks
  } else {
    names(regions) <- regions@elementMetadata[[region.name]]
  }

  ## Annotate regions
  regions.anno <- annotatePeak(regions, tssRegion = promoter.region, level = anno.level,
                               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                               annoDb = "org.Hs.eg.db")
  regions.anno <- as.data.frame(regions.anno)

  ## Find overlaps between peak1 and regions
  peak1.gr <- peak2granges(conns$Peak1)
  regions.conns.ix <- findOverlaps(peak1.gr, regions)

  ## Annotate peak2
  peak2.gr <- peak2granges(conns$Peak2)
  peak2.gr.anno <- annotatePeak(peak2.gr, tssRegion = promoter.region, level = anno.level,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                annoDb = "org.Hs.eg.db")
  peak2.gr.anno <- peak2.gr.anno@anno

  ## Get weights (otherwise set to 1 if null)
  if (!is.null(weight.col)) {
    stopifnot(weight.col %in% colnames(conns))
    peak2.gr.anno$weights <- conns[[weight.col]]
  } else {
    peak2.gr.anno$weights <- 1
  }

  regions.conns.ix.vector <- as.character(regions.conns.ix@to)
  names(regions.conns.ix.vector) <- as.character(regions.conns.ix@from)
  regions.conns.ix.list <- lapply(unique(regions.conns.ix.vector), function(x) {
    names(regions.conns.ix.vector[regions.conns.ix.vector == x])
  })
  names(regions.conns.ix.list) <- names(regions)[as.integer(unique(regions.conns.ix.vector))]

  ## Link region to a gene if the region is in the promoter region
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
      w <- peaks.anno.gr$weights
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
  region.gene.weights.df <- lapply(names(region.gene.weights.list), function(i) {
    weights <- region.gene.weights.list[[i]]
    data.frame(region = i, gene = names(weights), weight = weights,
               stringsAsFactors = F)
  })
  region.gene.weights.df <- do.call(rbind, region.gene.weights.df)
  rownames(region.gene.weights.df) <- NULL
  return(region.gene.weights.df)
}




#' Link TFs to peaks using motifs
#'
#' @param motif_ix Binary motif matrix of regions by TFs specifying which TFs have motifs in which regions
#' @return TF region motif overlaps in dataframe format
#' @export
#'
TFRegionLinks <- function(motif_ix) {
  tf.region.df <- lapply(colnames(motif_ix), function(i) {
    ix <- motif_ix[, i]
    tf.regions <- names(ix[ix > 0])
    data.frame(TF = i, region = tf.regions,
               stringsAsFactors = F)
  })
  tf.region.df <- do.call(rbind, tf.region.df)
  rownames(tf.region.df) <- NULL
  tf.region.df
}



#' Link TFs to genes
#'
#' @param tf.region.df Dataframe linking TFs to regions using motif data
#' @param region.gene.df Dataframe linking regions to genes using Hi-C or coaccessibility data
#' @return Dataframe linking TFs to genes
#' @export
#'
TFGeneLinks <- function(tf.region.df, region.gene.df) {
  tf.gene.df <- lapply(unique(tf.region.df$TF), function(tf) {
    tf.regions <- unique(subset(tf.region.df, TF == tf)$region)
    tf.genes <- unique(subset(region.gene.df, region %in% tf.regions)$gene)
    if (length(tf.genes) > 0) {
      return(data.frame(TF = tf, gene = tf.genes, stringsAsFactors = F))
    } else {
      return(NULL)
    }
  })
  tf.gene.df <- do.call(rbind, tf.gene.df)
  rownames(tf.gene.df) <- NULL
  tf.gene.df
}




#' Add correlations between region accessiblity and gene expression to region-gene links
#'
#' @param df Dataframe linking two types of features (i.e. TFs & genomic regions, genomic regions & genes, etc)
#' @param col1 Column name for the 1st feature (i.e. "region")
#' @param col2 Column name for the 2nd feature (i.e. "gene")
#' @param mat1 Matrix corresponding to feature 1. Columns must match columns of mat2.
#' @param mat2 Matrix corresponding to feature 2. Columns must match columns of mat1
#' @param n.cores Number of cores for parallel processing
#'
#' @return Dataframe with the "cor" column added
#'
#' @import snow
#' @export
#'
AddLinkCorrelations <- function(df, col1, col2, mat1, mat2, n.cores = 1) {
  stopifnot(all(colnames(mat1) %in% colnames(mat2)))
  mat1 <- mat1[,colnames(mat2)]
  df <- df[df[[col1]] %in% rownames(mat1) & df[[col2]] %in% rownames(mat2),]

  if (n.cores == 1) {
    col1.vec <- df[[col1]]
    col2.vec <- df[[col2]]
    df$cor <- sapply(1:nrow(df), function(i) {
      i1 <- col1.vec[[i]]
      i2 <- col2.vec[[i]]
      cor(mat1[i1,], mat2[i2,])
    })
  } else {
    df.list <- split(df, rep(1:ceiling(nrow(df)/n.cores), each = n.cores, length.out = nrow(df)))
    cl <- snow::makeCluster(n.cores, type = "SOCK")
    snow::clusterExport(cl, "col1", "col2", "mat1", "mat2", envir = environment())
    df.list <- snow::parLapply(cl, df.list, function(df.chunk) {
      col1.vec <- df.chunk[[col1]]
      col2.vec <- df.chunk[[col2]]
      df.chunk$cor <- sapply(1:nrow(df.chunk), function(i) {
        i1 <- col1.vec[[i]]
        i2 <- col2.vec[[i]]
        cor(mat1[i1,], mat2[i2,])
      })
      df.chunk
    })
    df <- do.call(rbind, df.list)
    rownames(df) <- NULL
  }
  df
}



#' Find all genes within specified distance of each genomic region
#'
#' @param gr Genomic regions
#' @param genes.anno Gene annotations generated by the geneAnno function
#' @param distance Distance to search for genes
#' @param chr.lengths Chromosome sizes
#' @param region.name Name of each region if in the GRanges metadata. Otherwise set to null.
#' @return List of genes within specified distance of genomic regions
#'
#' @export
#'
FindNearbyGenes <- function(regions, genes.anno, distance, chr.lengths, region.name = NULL) {
  if(is.null(region.name)) {
    names(regions) <- granges2peak(regions)
  } else {
    names(regions) <- regions@elementMetadata[[region.name]]
  }

  regions <- extend_gr(regions, buffer = distance, chr.lengths = chr.lengths)
  genes.ix <- findOverlaps(regions, genes.anno)
  genes.ix.vec <- as.character(names(regions)[genes.ix@from])
  names(genes.ix.vec) <- as.character(genes.anno$gene_symbol[genes.ix@to])

  regions.genes.list <- UnflattenGroups(genes.ix.vec)
  regions.genes.df <- lapply(names(regions.genes.list), function(i) {
    data.frame(region = i, gene = regions.genes.list[[i]],
               stringsAsFactors = F)
  })
  regions.genes.df <- do.call(rbind, regions.genes.df)
  rownames(regions.genes.df) <- NULL
  regions.genes.df
}
