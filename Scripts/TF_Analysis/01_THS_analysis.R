library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(chromfunks)
library(swne)
library(Seurat)
library(cisTopic)
library(ggplot2)
library(Matrix)
library(ggpubr)
library(cowplot)
library(ChIPseeker)
library(perturbLM)
library(liger)
library(reshape2)
library(TFBSTools)
library(motifmatchr)
library(igraph)
library(Signac)


## Global parameters
output.file <- "../Data/scTHS_analysis.RData"

## Set number of replicates
n.runs <- 20

## Cell type ordering
fetal.cl.order <- c("oRG", "vRG", "CycProg", "ipEx", "ebEx",
                    "ExL2/3", "ExL4", "ExL5/6", 
                    "InCGE", "InMGE", 
                    "OPC", "Mic")
adult.cl.order <- c("Ast", "Ex", "In", "Oli", "Opc", "Mic", "End")


#### Load scTHS-seq data ####

## Import fetal counts and peaks
ths <- readRDS("../Data/scTHS_chromatin_seurat_object.Robj")

fetal.peaks <- peak2granges(rownames(ths), delim = c("-", "-"))
names(fetal.peaks) <- granges2peak(fetal.peaks)

fetal.cell.counts <- GetAssayData(ths, slot = "counts")
fetal.cell.counts@x[fetal.cell.counts@x > 1] <- 1
rownames(fetal.cell.counts) <- names(fetal.peaks)
dim(fetal.cell.counts)

fetal.clusters <- factor(ths$sub.named.ident)
table(fetal.clusters)

fetal.counts <- getPseudobulk(fetal.cell.counts[,names(fetal.clusters)], fetal.clusters)
fetal.counts <- fetal.counts[rowSums(fetal.counts) > 1,]
fetal.peaks <- peak2granges(rownames(fetal.counts))
dim(fetal.counts)

fetal.norm.counts <- RunTFIDF(fetal.cell.counts)
fetal.norm.counts <- fetal.norm.counts[rownames(fetal.counts),]

## Load UMAP data
fetal.umap.emb <- Embeddings(ths, "umap")
dim(fetal.umap.emb)


## Import adult counts and peaks
adult.ths <- readRDS("../Data/adult_scTHS_chromatin_seurat_object.Robj")

adult.peaks <- peak2granges(rownames(adult.ths), delim = c("-", "-"))
names(adult.peaks) <- granges2peak(adult.peaks)

adult.cell.counts <- GetAssayData(adult.ths, slot = "counts")
adult.cell.counts@x[adult.cell.counts@x > 1] <- 1
rownames(adult.cell.counts) <- names(adult.peaks)
dim(adult.cell.counts)

adult.clusters <- adult.ths$ident
table(adult.clusters)

adult.counts <- getPseudobulk(adult.cell.counts, adult.clusters)
adult.counts <- adult.counts[rowSums(adult.counts) > 1,]
adult.peaks <- peak2granges(rownames(adult.counts))
dim(adult.counts)

## Clean up unneeded objects
rm(ths, adult.ths); invisible(gc())


## Build summarized experiment objects and correct for GC bias
fetal.SE <- SummarizedExperiment(assays = list(counts = fetal.counts),
                                 rowData = fetal.peaks, 
                                 colData = DataFrame(names = colnames(fetal.counts)))
fetal.SE <- addGCBias(fetal.SE, genome = BSgenome.Hsapiens.UCSC.hg38)

adult.SE <- SummarizedExperiment(assays = list(counts = adult.counts),
                                 rowData = adult.peaks, 
                                 colData = DataFrame(names = colnames(adult.counts)))
adult.SE <- addGCBias(adult.SE, genome = BSgenome.Hsapiens.UCSC.hg38)




#### Load and annotate HAR/HGEs ####

## Import HARs
# har.files <- list.files(path = "HAR_Regions/Formatted", full.names = T)
har.files <- c("../Data/HAR_HGE_bed_files/HAR_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/HGE_fetal_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/HGE_adult_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/HLE_adult_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/ENCODE_DNAse_regions_hg38.bed")

GetAnnoProp <- function(anno) {
  anno <- as.data.frame(anno)
  anno.type <- as.character(anno$annotation)
  anno.type <- sapply(anno.type, function(x) trimws(ExtractField(x, field = 1, delim = "\\(")))
  anno.df <- data.frame(group = names(table(anno.type)), value = as.integer(table(anno.type)))
  anno.df$percent <- as.numeric(table(anno.type)/sum(table(anno.type)))*100
  anno.df
}

## Custom pie chart function
PlotAnno <- function(anno) {
  anno.df <- GetAnnoProp(anno)
  ggplot(anno.df, aes(x = "", y = value, fill = group)) +
    geom_bar(width = 1, stat = "identity", position = position_stack()) + theme_classic() + 
    theme(axis.title = element_blank(), axis.text = element_blank(),
          legend.position = "none", axis.line = element_blank()) +
    coord_polar("y")
}


## Load HAR/HGE regions
har.list <- lapply(har.files, function(fi) {
  gr <- readPeakFile(fi)
  colnames(gr@elementMetadata) <- c("score", "HAR")
  gr
})
names(har.list) <- gsub("_regions_hg38.bed", "", har.files)
names(har.list) <- sapply(names(har.list), ExtractField, field = 4, delim = "/")

## Annotate regions and create plots
har.anno.list <- lapply(har.list, function(gr) {
  annotatePeak(gr, tssRegion = c(-3e3, 3e3), 
               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
               annoDb = "org.Hs.eg.db")
})

har.anno.list.frac <- lapply(har.anno.list, GetAnnoProp)
har.anno.list.plot <- lapply(har.anno.list, PlotAnno)

# pdf("fCTX-combined_anno_piecharts.pdf", width = 7, height = 3)
plot_grid(plotlist = har.anno.list.plot, nrow = 1)
# dev.off()

## Get proportions
har.anno.list.frac


## Save results
save.image(output.file)



#### Run HAR/HGE Analysis ####

## Run HAR enrichment
fetal.cl.HAR.z <- ChromRegionEnrich(fetal.SE, har.list, n.runs = n.runs)
adult.cl.HAR.z <- ChromRegionEnrich(adult.SE, har.list, n.runs = n.runs)

## Visualize results
HAR.z <- cbind(fetal.cl.HAR.z[,fetal.cl.order], adult.cl.HAR.z[,adult.cl.order])
# rownames(HAR.z) <- gsub("_regions_hg38", "", rownames(HAR.z))
# HAR.z <- HAR.z[c("HAR", "HGE_fetal", "HGE_adult", "HLE_adult", "HCNE.rheMac"),]


max.z <- 20
HAR.z[HAR.z > max.z] <- max.z
HAR.z[HAR.z < -1*max.z] <- -1*max.z
HAR.z <- HAR.z[,rev(colnames(HAR.z))]

min.z <- 6
pdf("fCTX-combined_gChromVar_HAR_HGE_HLE_zscores.pdf", width = 3, height = 6.25)
ggHeat(t(HAR.z), heatscale = c(low = 'deepskyblue', mid = 'white', high = 'darkorchid'),
       x.lab.size = 11, y.lab.size = 11, 
       dot.highlight.cutoff = min.z)
dev.off()

## Save results
saveRDS(HAR.z, file = "../Results/scTHS_HAR_HGE_enrich_zscores.Robj")
save.image(output.file)



#### Load scRNA-seq data ####

## Identify cluster specific expression of ChromVAR TFs

## Load fetal gene expression data
fbcx <- readRDS("../Data/snDrop_RNA_seurat_object.Robj")
fetal.rna.idents <- droplevels(fbcx$annot)
table(fetal.rna.idents)

fetal.rna.expr <- GetAssayData(fbcx, slot = "counts")
fetal.rna.expr[fetal.rna.expr > 0] <- 1
fetal.rna.cl.mean <- sapply(levels(fetal.rna.idents), function(cl) {
  Matrix::rowMeans(fetal.rna.expr[,names(fetal.rna.idents[fetal.rna.idents == cl])])
})

fbcx <- ScaleData(fbcx, features = rownames(fbcx))
fetal.rna.scale <- GetAssayData(fbcx, slot = "scale.data")
fetal.rna.cl.scale <- sapply(levels(fetal.rna.idents), function(cl) {
  Matrix::rowMeans(fetal.rna.scale[,names(fetal.rna.idents[fetal.rna.idents == cl])])
})


## Load adult gene expression data
vctx <- readRDS("../Data/adult_snDrop_seurat_object.Robj")
adult.rna.idents <- vctx$ident.use
table(adult.rna.idents)

adult.rna.expr <- GetAssayData(vctx, slot = "counts")
adult.rna.expr[adult.rna.expr > 0] <- 1
adult.rna.cl.mean <- sapply(levels(adult.rna.idents), function(cl) {
  cl.cells <- names(adult.rna.idents[adult.rna.idents == cl])
  Matrix::rowMeans(adult.rna.expr[,cl.cells])
})

## Clean up objects
rm(fbcx, vctx, fetal.rna.expr, adult.rna.expr); invisible(gc())



#### Run HAR/HGE chromVAR ####

## Run HAR chromVAR

## Get motifs
motifs <- getJasparMotifs()
fetal_motif_ix <- matchMotifs(motifs, fetal.SE, genome = BSgenome.Hsapiens.UCSC.hg38)
adult_motif_ix <- matchMotifs(motifs, adult.SE, genome = BSgenome.Hsapiens.UCSC.hg38)

## Compute TF motif z-scores
n.runs <- 20
register(MulticoreParam(8))
fetal.tf.HAR <- RegionChromVAR(fetal.SE, fetal_motif_ix, har.list, n.runs = n.runs)
adult.tf.HAR <- RegionChromVAR(adult.SE, adult_motif_ix, har.list, n.runs = n.runs)


## Bind fetal and adult HAR chromVAR results
adult.cl.order[grepl("Ex", adult.cl.order)] <- "Ex"
adult.cl.order <- unique(adult.cl.order)

HAR.tf.z <- lapply(names(har.list), function(h) {
  z <- cbind(fetal.tf.HAR[[h]][,fetal.cl.order], adult.tf.HAR[[h]][rownames(fetal.tf.HAR[[h]]), adult.cl.order])
  z[is.infinite(z) | is.nan(z)] <- NA
  rownames(z) <- sapply(rownames(z), ExtractField, delim = "_", field = 2)
  z
})
names(HAR.tf.z) <- names(har.list)


## Get fetal TF expression
fetal.tf.cl.mean <- lapply(rownames(HAR.tf.z[[1]]), function(tf) {
  tfs <- strsplit(tf, split = "::")[[1]]
  tfs <- tfs[!grepl("var", tfs)]
  if (all(tfs %in% rownames(fetal.rna.cl.mean))) {
    if (length(tfs) > 1) return(colMeans(as.matrix(fetal.rna.cl.mean[tfs,])))
    else return(fetal.rna.cl.mean[tfs,])
  }
  else {
    return(NULL)
  }
})
names(fetal.tf.cl.mean) <- rownames(HAR.tf.z[[1]])
fetal.tf.cl.mean <- do.call(rbind, fetal.tf.cl.mean)
fetal.tf.cl.mean <- fetal.tf.cl.mean[,fetal.cl.order]

## Get scaled single cell TF expression
fetal.tf.scale <- lapply(rownames(HAR.tf.z[[1]]), function(tf) {
  tfs <- strsplit(tf, split = "::")[[1]]
  tfs <- tfs[!grepl("var", tfs)]
  if (all(tfs %in% rownames(fetal.rna.scale))) {
    if (length(tfs) > 1) return(colMeans(as.matrix(fetal.rna.scale[tfs,])))
    else return(fetal.rna.scale[tfs,])
  }
  else {
    return(NULL)
  }
})
names(fetal.tf.scale) <- rownames(HAR.tf.z[[1]])
fetal.tf.scale <- do.call(rbind, fetal.tf.scale)


## Get adult TF expression
adult.tf.cl.mean <- lapply(rownames(HAR.tf.z[[1]]), function(tf) {
  tfs <- strsplit(tf, split = "::")[[1]]
  tfs <- tfs[!grepl("var", tfs)]
  if (all(tfs %in% rownames(adult.rna.cl.mean))) {
    if (length(tfs) > 1) return(colMeans(as.matrix(adult.rna.cl.mean[tfs,])))
    else return(adult.rna.cl.mean[tfs,])
  }
  else {
    return(NULL)
  }
})
names(adult.tf.cl.mean) <- rownames(HAR.tf.z[[1]])
adult.tf.cl.mean <- do.call(rbind, adult.tf.cl.mean)
adult.tf.cl.mean <- adult.tf.cl.mean[,adult.cl.order]


## Set non-expressed Z-scores to zero
min.expr <- 0.1
sum(fetal.tf.cl.mean > min.expr)/(length(fetal.tf.cl.mean))
sum(adult.tf.cl.mean > min.expr)/(length(adult.tf.cl.mean))

HAR.tf.z.expr <- lapply(HAR.tf.z, function(z) {
  fetal.z <- z[,fetal.cl.order]
  fetal.z[!rownames(fetal.z) %in% rownames(fetal.tf.cl.mean), ] <- NA
  fetal.z[rownames(fetal.tf.cl.mean),][fetal.tf.cl.mean < min.expr] <- NA
  
  adult.z <- z[,adult.cl.order]
  adult.z[!rownames(adult.z) %in% rownames(adult.tf.cl.mean), ] <- NA
  adult.z[rownames(adult.tf.cl.mean),][adult.tf.cl.mean < min.expr] <- NA
  
  cbind(fetal.z, adult.z)
})


## Save results
save(fetal_motif_ix, fetal.tf.HAR, HAR.tf.z.expr,
     file = "../Results/scTHS_HAR_HGE_chromVAR.RData")

## Save results
save.image(output.file)



#### Analyze HAR/HGE TFs ####

## Look at correlations between HAR types
HAR_cor <- sapply(HAR.tf.z.expr, function(z1) {
  sapply(HAR.tf.z.expr, function(z2) {
    z1 <- flattenMat(z1); z2 <- flattenMat(z2);
    z1 <- z1[!is.na(z1)]; z2 <- z2[!is.na(z2)];
    ix <- intersect(names(z1), names(z2))
    cor(z1[ix], z2[ix])
  })
})
HAR_cor

## Find top TFs for each HAR type
min.tf.z <- 3
n.tfs <- 5

sig.cl <- colnames(HAR.z)[apply(HAR.z, 2, function(x) any(x > min.z))]
top_HAR_TFs <- lapply(names(HAR.tf.z.expr), function(h) {
  z <- HAR.tf.z.expr[[h]]
  z[is.na(z)] <- 0
  
  top.tfs <- lapply(sig.cl, function(cl) {
    z.cl <- z[,cl]
    z.cl <- z.cl[z.cl > min.tf.z]
    head(names(sort(z.cl, decreasing = T)), n = n.tfs)
  })
  
  top.tfs <- unique(unlist(top.tfs, F, F))
  if (length(top.tfs) == 1) {
    top.z <- t(as.matrix(z[top.tfs,]))
    rownames(top.z) <- top.tfs
  } else {
    top.z <- z[top.tfs,]
  }
  top.z
})
names(top_HAR_TFs) <- names(HAR.tf.z.expr)

chrom.heatscale <- c(low = 'deepskyblue', mid = 'white', high = 'darkorchid')
top_HAR_TFs_mat <- do.call(rbind, top_HAR_TFs)

# pdf("fCTX-combined_HAR_top_TF_heatmap.pdf", width = 6, height = 5.75)
ggHeat(top_HAR_TFs_mat[rev(rownames(top_HAR_TFs_mat)),], heatscale = chrom.heatscale, 
       x.lab.size = 11, y.lab.size = 10)
# dev.off()

## Guide for which TFs belong to which HAR/HGE type
lapply(top_HAR_TFs, rownames)

## Save results
save.image(output.file)



#### Link HAR/HGE peaks to genes using Hi-C ####

## Get peaks overlapping HARs/HGEs
fetal.peak.scores <- RegionOverlapScores(rowRanges(fetal.SE), har.list)
fetal.peak.maps <- lapply(colnames(fetal.peak.scores), function(i) {
  ix <- fetal.peak.scores[,i]
  names(ix[ix > 0])
})
names(fetal.peak.maps) <- colnames(fetal.peak.scores)

## Get normalized accessibility of TFs
fetal.region.cl.counts <- GetRegionCounts(fetal.norm.counts[,names(fetal.clusters)], fetal.peak.maps,
                                          idents = fetal.clusters, scale = F)
lapply(fetal.region.cl.counts, dim)


## Load Hi-C peaks
eN.conns <- read.table("../Data/Fetal_HiC_Data/eN.MAPS.peaks.txt", sep = "\t", 
                       header = T, stringsAsFactors = F)
iN.conns <- read.table("../Data/Fetal_HiC_Data/iN.MAPS.peaks.txt", sep = "\t", 
                       header = T, stringsAsFactors = F)
IPC.conns <- read.table("../Data/Fetal_HiC_Data/IPC.MAPS.peaks.txt", sep = "\t", 
                        header = T, stringsAsFactors = F)
RG.conns <- read.table("../Data/Fetal_HiC_Data/RG.MAPS.peaks.txt", sep = "\t", 
                       header = T, stringsAsFactors = F)

## Load chromosome lengths
hg38.chr.lengths <- read.table("../Data/hg38.chr.lengths.txt", header = F, sep = "\t")

## Make list of connections
conns.list <- list(RG = RG.conns,
                   ipEx = IPC.conns,
                   ebEx = eN.conns,
                   In = iN.conns)

conns.list <- lapply(conns.list, function(cl.conns) {
  cl.conns$Peak1 <- paste0(cl.conns$chr1, ":", cl.conns$start1, "-", cl.conns$end1)
  cl.conns$Peak2 <- paste0(cl.conns$chr2, ":", cl.conns$start2, "-", cl.conns$end2)
  cl.conns[,c("Peak1", "Peak2")]
})


## Align scTHS-seq clusters and sorted Hi-C clusters
HiC.cl.map <- c("oRG" = "RG", "vRG" = "RG", "CycProg" = "ipEx", "ipEx"="ipEx",
                "ebEx" = "ebEx", "InMGE"="In", "InCGE"="In")

## Create accessibility matrices with sorted Hi-C clusters
fetal.region.HiC.counts <- lapply(fetal.region.cl.counts, function(cl.mat) {
  cl.mat <- cl.mat[,names(HiC.cl.map)]
  t(apply(cl.mat, 1, function(x) tapply(x, HiC.cl.map, max)))
})

## Select accessibility cutoff
fetal.region.counts.vec <- unlist(lapply(fetal.region.cl.counts, as.vector))
hist(fetal.region.counts.vec)

min.access <- 0.01
sum(fetal.region.counts.vec > min.access)/length(fetal.region.counts.vec)

## Find genes linked to peaks that overlap HARs/HGEs and that are accessible in the relevant cell type
fetal.region.gene.contacts <- lapply(names(fetal.peak.maps), function(h) {
  peaks <- fetal.peak.maps[[h]]
  har.peaks <- fetal.peaks[peaks]
  df.list <- lapply(names(conns.list), function(cl) {
    cl.conns <- conns.list[[cl]]
    
    ## Get links between regions and genes
    df <- RegionGeneLinks(har.peaks, cl.conns, link.promoter = T, promoter.region = c(-3e3, 3e3),
                          region.name = NULL, weight.col = NULL)
    
    ## Add correlations between region accessibility and gene expression
    df <- AddLinkCorrelations(df, "region", "gene", 
                              fetal.region.cl.counts[[h]][,fetal.cl.order], 
                              fetal.rna.cl.scale[,fetal.cl.order])
    
    ## Filter for regions that are accessible in the given cell type
    h.access <- fetal.region.HiC.counts[[h]]
    access.regions <- unique(df$region)
    access.regions <- access.regions[access.regions %in% rownames(h.access)]
    access.regions <- access.regions[h.access[access.regions, cl] > min.access]
    subset(df, region %in% access.regions)
  })
  names(df.list) <- names(conns.list)
  df.list
})
names(fetal.region.gene.contacts) <- names(fetal.peak.maps)

sapply(fetal.region.gene.contacts, function(cl.contacts) {
  sapply(cl.contacts, nrow)
})

## Save results
rm(fetal.norm.counts); invisible(gc());
save.image(output.file)




#### Link TFs to HARs/HGEs using motifs and correlations ####

## Create binned counts for improved correlations
fetal.bin.counts <- GetBinCounts(fetal.cell.counts, fetal.umap.emb)
fetal.bin.norm.counts <- RunTFIDF(fetal.bin.counts)

## Subset to peaks overlapping HAR/HGE regions to save memory
fetal.mapped.peaks <- unique(unlist(fetal.peak.maps, F, F))
fetal.bin.counts <- fetal.bin.counts[fetal.mapped.peaks,]
fetal.bin.norm.counts <- fetal.bin.norm.counts[fetal.mapped.peaks,]

## Run chromVAR at the single cell level for correlating TF activity to HAR/HGE peaks with TF motifs
fetal.bin.SE <- SummarizedExperiment(assays = list(counts = fetal.bin.counts),
                                    rowData = peak2granges(rownames(fetal.bin.counts)), 
                                    colData = DataFrame(names = colnames(fetal.bin.counts)))
fetal.bin.SE <- addGCBias(fetal.bin.SE, genome = BSgenome.Hsapiens.UCSC.hg38)

register(MulticoreParam(6))
fetal.tf.bin.HAR <- RegionChromVAR(fetal.bin.SE, fetal_motif_ix, har.list, n.runs = 3)
fetal.tf.bin.HAR <- lapply(fetal.tf.bin.HAR, function(z) {
  z[is.na(z) | is.nan(z) | is.infinite(z)] <- 0
  rownames(z) <- sapply(rownames(z), ExtractField, 2)
  z
})


## Format TF motif names
fetal_motif_ix_mat <- assay(fetal_motif_ix)
colnames(fetal_motif_ix_mat) <- sapply(colnames(fetal_motif_ix_mat), ExtractField, field = 2)

## Link TFs to regions using motifs
fetal.TF.all.regions <- lapply(names(fetal.peak.maps), function(h) {
  peaks <- fetal.peak.maps[[h]]
  tf.peak.df <- TFRegionLinks(fetal_motif_ix_mat[peaks,])
  AddLinkCorrelations(tf.peak.df, "TF", "region", 
                      fetal.tf.bin.HAR[[h]][,colnames(fetal.bin.counts)],
                      fetal.bin.norm.counts[,colnames(fetal.bin.counts)],
                      n.cores = 4)
})
names(fetal.TF.all.regions) <- names(fetal.peak.maps)
sapply(fetal.TF.all.regions, nrow)

## Distribution of number of regions per TF
hist(unlist(lapply(fetal.TF.all.regions, function(df) as.vector(table(df$TF)))))

## Distribution of correlations between TF activity and region accessibility
hist(unlist(lapply(fetal.TF.all.regions, function(df) df$cor), F, F))
min.tf.peak.cor <- 0.1

fetal.TF.regions <- lapply(fetal.TF.all.regions, function(df) {
  subset(df, cor > min.tf.peak.cor & !is.na(cor))
})
sapply(fetal.TF.regions, nrow)

## Distribution of number of regions per TF
hist(unlist(lapply(fetal.TF.regions, function(df) as.vector(table(df$TF)))))

## Save results
save.image(output.file)



#### Build HAR/HGE TF networks ####

## Link TFs to genes using cell type specific Hi-C
min.peak.gene.cor <- 0
fetal.HiC.tf.expr <- t(apply(fetal.tf.cl.mean, 1, function(x) tapply(x[names(HiC.cl.map)], HiC.cl.map, max)))
fetal.TF.cl.genes <- lapply(names(fetal.region.gene.contacts), function(h) {
  peak.gene.list <- fetal.region.gene.contacts[[h]]
  tf.peak.df <- fetal.TF.regions[[h]]
  
  tf.gene.list <- lapply(names(peak.gene.list), function(cl) {
    peak.gene.df <- subset(peak.gene.list[[cl]], cor > min.peak.gene.cor & !is.na(cor))
    tf.genes <- TFGeneLinks(tf.peak.df, peak.gene.df)
    
    cl.expr <- fetal.HiC.tf.expr[,cl]
    expr.tfs <- names(cl.expr[cl.expr >= min.expr])
    tf.genes <- subset(tf.genes, TF %in% expr.tfs)
    if (nrow(tf.genes) > 0) {
      tf.genes <- AddLinkCorrelations(tf.genes, "TF", "gene", 
                                      fetal.tf.scale[,names(fetal.rna.idents)], 
                                      fetal.rna.scale[,names(fetal.rna.idents)],
                                      n.cores = 8)
      tf.genes <- subset(tf.genes, !is.na(cor))
    } else {
      tf.genes <- NULL
    }
    tf.genes
  })
  names(tf.gene.list) <- names(peak.gene.list)
  tf.gene.list
})
names(fetal.TF.cl.genes) <- names(fetal.region.gene.contacts)
sapply(fetal.TF.cl.genes, function(df.list) sapply(df.list, nrow))

## Export results
saveRDS(fetal.TF.cl.genes, file = "../Results/scTHS_TF_gene_network_dfs.Robj")

## Save full image
save.image(output.file)



#### Analyze HAR/HGE TF networks ####

## Create force directed TF - gene network
fetal.TF.cl.graphs <- lapply(names(fetal.TF.cl.genes), function(h) {
  z <- HAR.tf.z.expr[[h]][,names(HiC.cl.map)]
  z[is.na(z)] <- -Inf
  z <- t(apply(z, 1, function(x) tapply(x, HiC.cl.map, max)))
  
  cl.graphs <- lapply(names(fetal.TF.cl.genes[[h]]), function(cl) {
    z.cl <- z[,cl]
    har.tfs <- names(z.cl[z.cl > min.tf.z])
    
    df <- fetal.TF.cl.genes[[h]][[cl]]
    df$edge.label <- paste0(df$TF, "_", df$gene)
    # df <- subset(df, abs(cor) > min.tf.gene.cor)
    g <- graph_from_data_frame(df, directed = TRUE, vertices = NULL)
    V(g)$color <- sapply(names(V(g)), function(x) {
      if (x %in% har.tfs) return("tomato")
      else return("skyblue")
    })
    V(g)$size <- sapply(names(V(g)), function(x) {
      if (x %in% har.tfs) return(6)
      else return(2)
    })
    V(g)$label <- sapply(names(V(g)), function(x) {
      if (x %in% har.tfs) return(x)
      else return("")
    })
    g
  })
  names(cl.graphs) <- names(fetal.TF.cl.genes[[h]])
  cl.graphs[c("RG", "ipEx", "ebEx", "In")]
})
names(fetal.TF.cl.graphs) <- names(fetal.TF.cl.genes)


fetal.TF.cl.layouts <- lapply(fetal.TF.cl.graphs, function(g.list) {
  lapply(g.list, function(g) {
    layout_with_fr(g)
  })
})


## Plot cell type specific graphs
h.type <- "HGE_fetal"
h.graph.list <- fetal.TF.cl.graphs[[h.type]]
h.layout.list <- fetal.TF.cl.layouts[[h.type]]

pdf(paste0("fCTX-combined_HAR_TF_gene_", h.type, "_network_nolabels.pdf"), width = 6.5, height = 5)
par(mfrow = c(2,2), mar = c(0.25, 0.1, 1, 0.1))
h.plots <- lapply(names(h.graph.list), function(cl) {
  g <- h.graph.list[[cl]]
  l <- h.layout.list[[cl]]
  plot.igraph(g, layout = l, vertex.size = V(g)$size, vertex.label = NA,
              vertex.color = V(g)$color, edge.arrow.size = 0.1, main = cl, 
              vertex.label.cex = 0.6)
})
dev.off()

lapply(h.graph.list, function(g) {
  labels <- V(g)$label
  labels[labels != ""]
})


## Save results
saveRDS(fetal.TF.cl.graphs, file = "../Results/scTHS_TF_gene_igraphs.Robj")
save.image(output.file)



#### Export Supplemental Files ####

## Map peaks to HARs/HGEs
fetal.peaks$Peak <- names(fetal.peaks)
fetal.peak.region.mappings <- lapply(har.list, function(gr) {
  peak.region.list <- RegionOverlapList(fetal.peaks, gr, gr1.name = "Peak")
  peak.region.list <- lapply(peak.region.list, function(gr) {
    paste0(granges2peak(gr), collapse = ",")
  })
  hash.map <- new.env(hash = T)
  for (p in names(peak.region.list)) {
    hash.map[[p]] <- peak.region.list[[p]]
  }
  hash.map
})

## Export fetal HAR/HGE to gene links
for (h in names(fetal.region.gene.contacts)) {
  region.mapping <- fetal.peak.region.mappings[[h]]
  df.list <- fetal.region.gene.contacts[[h]]
  
  df.list <- lapply(names(df.list), function(cl) {
    df <- df.list[[cl]]
    df$cell.type <- cl
    df$accelerated.region <- sapply(df$region, function(p) region.mapping[[p]])
    df
  })
  df <- do.call(rbind, df.list)
  df$ths.region <- df$region
  df <- df[,c("cell.type", "ths.region", "accelerated.region", "gene")]
  print(head(df))
  
  fname <- paste0("../Results/", h, "_", "THS_region_gene_network.tsv")
  write.table(df, file = fname, sep = "\t", row.names = F)
}
saveRDS(fetal.region.gene.contacts, file = "../Results/scTHS_region_gene_dfs.Robj")


## Export fetal TF to HAR/HGE links
for (h in names(fetal.TF.regions)) {
  region.mapping <- fetal.peak.region.mappings[[h]]
  df <- fetal.TF.regions[[h]]
  df$accelerated.region <- sapply(df$region, function(p) region.mapping[[p]])
  df$ths.region <- df$region
  df <- df[,c("TF", "ths.region", "accelerated.region")]
  print(head(df))
  
  fname <- paste0("../Results/", h, "_", "THS_TF_region_network.tsv")
  write.table(df, file = fname, sep = "\t", row.names = F, quote = F)
}
saveRDS(fetal.TF.regions, file = "../Results/scTHS_TF_region_dfs.Robj")


## Export fetal TF to gene links
for (h in names(fetal.TF.cl.genes)) {
  region.mapping <- fetal.peak.region.mappings[[h]]
  df.list <- fetal.TF.cl.genes[[h]]
  df.list <- lapply(names(df.list), function(cl) {
    df <- df.list[[cl]]
    df$cell.type <- cl
    df
  })
  df <- do.call(rbind, df.list)
  df$ths.region <- df$region
  df <- df[,c("cell.type", "TF", "gene")]
  print(head(df))
  
  fname <- paste0("../Results/", h, "_", "THS_TF_gene_network.tsv")
  write.table(df, file = fname, sep = "\t", row.names = F, quote = F)
}


## Get network stats
h.type <- "HGE_fetal"

sapply(fetal.TF.cl.genes[[h.type]], function(df) length(unique(df$TF)))
sapply(fetal.TF.cl.genes[[h.type]], function(df) mean(table(df$TF)))
mean(table(fetal.TF.regions[[h.type]]$TF))

sapply(fetal.region.gene.contacts[[h.type]], function(df) length(unique(df$region)))
sapply(fetal.region.gene.contacts[[h.type]], function(df) mean(table(df$region)))


