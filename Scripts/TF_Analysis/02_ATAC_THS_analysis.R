library(swne)
library(Seurat)
library(cisTopic)
library(chromfunks)
library(cellMapper)
library(ggplot2)
library(ChIPseeker)
library(SummarizedExperiment)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)
library(chromVAR)
library(gchromVAR)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(igraph)
library(reshape2)
library(perturbLM)
library(cowplot)
library(Signac)
library(liger)


## Output file
output.file <- "../Data/ATAC_THS_analysis.RData"

## Replicates for chromVAR runs
n.runs <- 20

## HAR/HGE naming conventions
har.name.map <- c("ENCODE_DNAse" = "ENCODE Ctrl", "HLE_adult" = "Adult HLE",
                  "HGE_adult" = "Adult HGE", "HGE_fetal" = "Fetal HGE")

## Cell type ordering for visualizations
fetal.cl.order <- c("oRG", "vRG", "CycProg", "ipEx", "ebEx",
                    "ExL2/3", "ExL4", "ExL5/6", 
                    "InCGE", "InMGE", 
                    "OPC", "Mic")
adult.cl.order <- c("Ast", "Ex", "In", 
                    "Oli", "Opc", 
                    "Mic", "End")
atac.cl.order <- c("RG", "ipEx", "RG/OPC", "ebEx", "Ex", "In", "End")



#### Load ATAC datasets ####

## Cluster mapping
atac.ths.map <- c("Astrocytes"="RG", 
                  "Astrocytes/Oligodendrocytes"="RG/OPC",
                  "Cerebrum_Unknown.3"="RG/ipEx", 
                  "Excitatory neurons"="Ex",
                  "Inhibitory neurons"="In", 
                  "Limbic system neurons"="LimbicNeurons",
                  "SKOR2_NPSR1 positive cells"="SKOR2_NPSR1_Cells",
                  "Vascular endothelial cells"="End")

## Load atac dataset with peaks lifted to hg38
atac <- readRDS("../Data/sciATAC_seurat_object.Robj")
atac.mapped.ident <- atac$mapped.ident
table(atac.mapped.ident)

## Filter peaks
min.cells <- 50
atac.cell.counts <- GetAssayData(atac, slot = "counts")
atac.cell.counts <- atac.cell.counts[rowSums(atac.cell.counts) > min.cells,]
atac.peaks <- peak2granges(rownames(atac.cell.counts), delim = c("-", "-"))
names(atac.peaks) <- rownames(atac.cell.counts) <- granges2peak(atac.peaks)
dim(atac.cell.counts)

## Get pseudobulk
atac.counts <- getPseudobulk(atac.cell.counts, atac.mapped.ident)

## Get TF-IDF normalized counts
atac.norm.counts <- RunTFIDF(atac.cell.counts)


## Load and format matched RNA data
atac.rna <- readRDS("../Data/sciRNA_seurat_object.Robj")
dim(atac.rna)

atac.rna.ident <- Idents(atac.rna)
atac.rna.ident <- plyr::revalue(atac.rna.ident, replace = c("OPC"="RG/OPC"))
table(atac.rna.ident)

atac.rna.expr <- GetAssayData(atac.rna, slot = "counts")
atac.rna.expr@x[atac.rna.expr@x > 1] <- 1
atac.rna.cl.mean <- sapply(levels(atac.rna.ident), function(cl) {
  Matrix::rowMeans(atac.rna.expr[,names(atac.rna.ident[atac.rna.ident == cl])])
})

## Load THS fetal RNA data
fbcx <- readRDS("../Data/snDrop_RNA_seurat_object.Robj")
fetal.rna.idents <- droplevels(fbcx$annot)
table(fetal.rna.idents)

## Clean up unused objects
rm(atac.rna.expr); invisible(gc());



#### Load HAR regions ####

## HAR files
har.files <- c("../Data/HAR_HGE_bed_files/HAR_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/HGE_fetal_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/HGE_adult_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/HLE_adult_regions_hg38.bed",
               "../Data/HAR_HGE_bed_files/ENCODE_DNAse_regions_hg38.bed")

## Load HAR/HGE regions
har.list <- lapply(har.files, function(fi) {
  gr <- readPeakFile(fi)
  colnames(gr@elementMetadata) <- c("score", "HAR")
  gr
})
names(har.list) <- gsub("_regions_hg38.bed", "", har.files)
names(har.list) <- sapply(names(har.list), ExtractField, field = 4, delim = "/")

## Get peaks overlapping HARs/HGEs
atac.peak.scores <- RegionOverlapScores(atac.peaks, har.list)
atac.peak.maps <- lapply(colnames(atac.peak.scores), function(h) {
  ix <- atac.peak.scores[,h]
  names(ix[ix > 0])
})
names(atac.peak.maps) <- colnames(atac.peak.scores)
sapply(atac.peak.maps, length)


#### Run ATAC HAR enrichment ####

## Create SummarizedExperiment
atac.SE <- SummarizedExperiment(assays = list(counts = atac.counts),
                                rowData = atac.peaks, 
                                colData = DataFrame(names = colnames(atac.counts)))
atac.SE <- addGCBias(atac.SE, genome = BSgenome.Hsapiens.UCSC.hg38)

## Compute multiple samples of weighted deviations
atac.HAR.z <- ChromRegionEnrich(atac.SE, har.list, n.runs = n.runs)

## Visualize results
atac.HAR.z <- atac.HAR.z
rownames(atac.HAR.z) <- gsub("_regions_hg38", "", rownames(atac.HAR.z))
# HAR.z <- HAR.z[c("HAR", "HGE_fetal", "HGE_adult", "HLE_adult", "HCNE.rheMac"),]

ths.HAR.z <- readRDS("../Results/scTHS_HAR_HGE_enrich_zscores.Robj")
HAR.z <- cbind(ths.HAR.z[rownames(atac.HAR.z), rev(colnames(ths.HAR.z))], 
               atac.HAR.z[,atac.cl.order])
# HAR.z <- atac.HAR.z[,rev(colnames(atac.HAR.z))]

max.z <- 20
HAR.z.plot <- HAR.z
HAR.z.plot[HAR.z.plot > max.z] <- max.z
HAR.z.plot[HAR.z.plot < -1*max.z] <- -1*max.z
rownames(HAR.z.plot) <- plyr::revalue(rownames(HAR.z.plot), replace = har.name.map)
HAR.z.plot <- HAR.z.plot[rev(rownames(HAR.z.plot)),]

min.z <- 6
pdf("fCTX-combined_gChromVar_HAR_HGE_HLE_zscores.pdf", width = 7.75, height = 2.25)
ggHeat(HAR.z.plot, heatscale = c(low = 'deepskyblue', mid = 'white', high = 'darkorchid'),
       x.lab.size = 10, y.lab.size = 11, 
       dot.highlight.cutoff = min.z)
dev.off()

## Save results
save.image(output.file)


#### Run ATAC HAR TF analysis ####

## Get motifs
motifs <- getJasparMotifs()
atac_motif_ix <- matchMotifs(motifs, atac.SE, genome = BSgenome.Hsapiens.UCSC.hg38)

## Computing deviations
register(MulticoreParam(4))
atac.HAR.tf.z <- RegionChromVAR(atac.SE, atac_motif_ix, har.list, n.runs = n.runs)
atac.HAR.tf.z <- lapply(atac.HAR.tf.z, function(z) {
  z[is.infinite(z) | is.nan(z)] <- NA
  rownames(z) <- sapply(rownames(z), ExtractField, field = 2)
  z
})

## Get atac TF expression
atac.tf.expr <- lapply(rownames(atac.HAR.tf.z[[1]]), function(tf) {
  tfs <- strsplit(tf, split = "::")[[1]]
  tfs <- tfs[!grepl("var", tfs)]
  if (all(tfs %in% rownames(atac.rna.cl.mean))) {
    if (length(tfs) > 1) return(colMeans(as.matrix(atac.rna.cl.mean[tfs,])))
    else return(atac.rna.cl.mean[tfs,])
  }
  else {
    return(NULL)
  }
})
names(atac.tf.expr) <- rownames(atac.HAR.tf.z[[1]])
atac.tf.expr <- do.call(rbind, atac.tf.expr)
atac.tf.expr <- atac.tf.expr[,atac.cl.order]

## Set non-expressed Z-scores to zero
min.expr <- 0.025
sum(atac.tf.expr > min.expr)/(length(atac.tf.expr))

atac.HAR.tf.z.expr <- lapply(atac.HAR.tf.z, function(z) {
  atac.z <- z[,atac.cl.order]
  atac.z[!rownames(atac.z) %in% rownames(atac.tf.expr), ] <- NA
  atac.z[rownames(atac.tf.expr),][atac.tf.expr < min.expr] <- NA
  atac.z
})



#### Analyze ATAC HAR TF enrichment ####

## Load THS chromVAR results
load("../Results/scTHS_HAR_HGE_chromVAR.RData")

## Merge THS and ATAC results
merged.HAR.tf.z.expr <- lapply(names(HAR.tf.z.expr), function(h) {
  cbind(HAR.tf.z.expr[[h]], atac.HAR.tf.z.expr[[h]][,atac.cl.order])
})
names(merged.HAR.tf.z.expr) <- names(HAR.tf.z.expr)

## TF Z-score cutoff
min.tf.z <- 4
max.tf.z <- 10

## Find number of significant TFs for each HAR type and dataset
sig.tf.summary <- lapply(names(merged.HAR.tf.z.expr), function(h) {
  z <- merged.HAR.tf.z.expr[[h]]
  z[is.na(z)] <- 0
  fetal.ths.n <- sum(apply(z[,fetal.cl.order], 1, function(x) any(x > min.tf.z)))
  adult.ths.n <- sum(apply(z[,adult.cl.order], 1, function(x) any(x > min.tf.z)))
  fetal.atac.n <- sum(apply(z[,atac.cl.order], 1, function(x) any(x > min.tf.z)))
  data.frame(N = c(fetal.ths.n, adult.ths.n, fetal.atac.n),
             Dataset = c("Fetal scTHS", "Adult scTHS", "sciATAC"),
             Region_Type = h)
})
sig.tf.summary <- do.call(rbind, sig.tf.summary)
sig.tf.summary$Region_Type <- plyr::revalue(sig.tf.summary$Region_Type, replace = har.name.map)
# sig.tf.summary$Region_Type <- factor(sig.tf.summary$Region_Type, levels = rownames(HAR.z.plot))

pdf("fCTX-combined_num_top_TFs.pdf", width = 4, height = 3.5)
ggplot(data = sig.tf.summary, aes(x = Region_Type, y = N, fill = Dataset)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 12, angle = 90, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "top")
dev.off()

## Find top TFs for each HAR type and plot
n.tfs <- 5
top_HAR_TFs <- lapply(names(merged.HAR.tf.z.expr), function(h) {
  HAR.z.vec <- HAR.z[h, ]
  sig.cl <- names(HAR.z.vec[HAR.z.vec > min.z])
  # sig.cl <- colnames(HAR.z)
  
  z <- merged.HAR.tf.z.expr[[h]]
  z[is.na(z)] <- 0
  z[z > max.tf.z] <- max.tf.z
  z[z < -1*max.tf.z] <- -1*max.tf.z
  
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
names(top_HAR_TFs) <- names(merged.HAR.tf.z.expr)

chrom.heatscale <- c(low = 'deepskyblue', mid = 'white', high = 'darkorchid')
# top_HAR_TFs_mat <- do.call(rbind, top_HAR_TFs)
top_HAR_TFs_mat <- top_HAR_TFs$HGE_fetal

pdf("fCTX-combined_Fetal_HGE_top_TF_heatmap.pdf", width = 7, height = 5.25)
ggHeat(top_HAR_TFs_mat, heatscale = chrom.heatscale, 
       x.lab.size = 10, y.lab.size = 10)
dev.off()

## Guide for which TFs belong to which HAR/HGE type
lapply(top_HAR_TFs, rownames)

## Save results
save(atac_motif_ix, atac.HAR.tf.z, atac.HAR.tf.z.expr,
     file = "../Results/sciATAC_HAR_HGE_chromVAR.RData")
save.image(output.file)



#### Validate ATAC HAR TF analysis using human/chimp dA peaks in organoids ####
load("../Data/human_chimp_dA_motif_enrich.RData")

human.motif.logFC <- log2(human.motif$fold.enrichment)
names(human.motif.logFC) <- human.motif$motif

min.tf.z <- 0
har.org.cor <- lapply(names(merged.HAR.tf.z.expr), function(h.type) {
  apply(merged.HAR.tf.z.expr[[h.type]], 2, function(z1) {
    z1 <- z1[!is.na(z1) & !is.nan(z1)]
    cor(z1, human.motif.logFC[names(z1)])
    if (any(z1 > min.tf.z)) return(cor(z1, human.motif.logFC[names(z1)]))
    else return(0)
  })
})
har.org.cor <- do.call(rbind, har.org.cor)
rownames(har.org.cor) <- plyr::revalue(names(merged.HAR.tf.z.expr), replace = har.name.map)
har.org.cor[abs(HAR.z[,colnames(har.org.cor)]) < min.z] <- 0

pdf("fCTX-combined_HAR_TF_human_chimp_dA_peak_cor.pdf", width = 7.75, height = 2.25)
ggHeat(har.org.cor[rev(rownames(har.org.cor)),], x.lab.size = 11, y.lab.size = 11,
       dot.highlight.cutoff = 0.3)
dev.off()

h.type <- "HGE_fetal"
cl <- "vRG"
cl.z <- merged.HAR.tf.z.expr[[h.type]][,cl]
cl.z <- cl.z[!is.na(cl.z)]

pts.label <- cl.z > 2.5 | human.motif.logFC[names(cl.z)] > 0.4
pdf(paste0("fCTX-combined_HAR_TF_", h.type, "_", cl, "_cor.pdf"), width = 4, height = 3.75)
PlotCorrelation(cl.z, human.motif.logFC[names(cl.z)], box = F, show.corr = T,
                pts.label = pts.label, use.label = T)
dev.off()

## Save results
save.image(output.file)



#### Compare ATAC and THS TF results ####

min.tf.z <- 4
har.tf.cor.list <- lapply(names(HAR.tf.z.expr), function(h) {
  ths_z_expr <- HAR.tf.z.expr[[h]]
  atac_z_expr <- atac.HAR.tf.z.expr[[h]][,atac.cl.order]
  
  cor.mat <- sapply(1:ncol(ths_z_expr), function(i) {
    sapply(1:ncol(atac_z_expr), function(j) {
      ths.z <- ths_z_expr[,i]
      ths.z <- ths.z[!is.na(ths.z)]
      atac.z <- atac_z_expr[names(ths.z),j]
      atac.z <- atac.z[!is.na(atac.z)]
      tfs.use <- intersect(names(ths.z), names(atac.z))
      
      if(any(ths.z[tfs.use] > min.tf.z) || any(atac.z[tfs.use] > min.tf.z)) {
        return(cor(ths.z[tfs.use], atac.z[tfs.use]))
      } else {
        return(0)
      }
    })
  }) 
  colnames(cor.mat) <- colnames(ths_z_expr)
  rownames(cor.mat) <- colnames(atac_z_expr)
  # cor.mat[rev(rownames(cor.mat)),]
  t(cor.mat[,rev(colnames(cor.mat))])
})
names(har.tf.cor.list) <- names(HAR.tf.z.expr)

pdf("fCTX-combined_ATAC_THS_Fetal_HGE_TF_activity_cor.pdf", width = 3.5, height = 5.75)
ggHeat(har.tf.cor.list$HGE_fetal, dot.highlight.cutoff = 0.3)
dev.off()

h.type <- "HGE_fetal"
atac.cl <- "RG"
ths.cl <- "oRG"

atac.z <- atac.HAR.tf.z.expr[[h.type]][,atac.cl]
atac.z <- atac.z[!is.na(atac.z)]
ths.z <- HAR.tf.z.expr[[h.type]][,ths.cl]
ths.z <- ths.z[!is.na(ths.z)]

pts.ix <- intersect(names(atac.z), names(ths.z))
pts.label <- atac.z[pts.ix] > 2 | ths.z[pts.ix] > 3
pdf(paste0("fCTX-combined_ATAC_THS_Fetal_HGE_RG_oRG_cor.pdf"), width = 4, height = 3.75)
PlotCorrelation(ths.z[pts.ix], atac.z[pts.ix], box = F, show.corr = T,
                pts.label = pts.label, use.label = T)
dev.off()


## Save results
save.image(output.file)



#### Link ATAC HAR/HGE peaks to genes ####

## Get normalized accessibility of TFs
atac.region.cl.counts <- GetRegionCounts(atac.norm.counts[,names(atac.mapped.ident)], atac.peak.maps, 
                                         idents = atac.mapped.ident, scale = F)
lapply(atac.region.cl.counts, dim)

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
atac.HiC.cl.map <- c("RG" = "RG", "ipEx" = "ipEx", "ebEx" = "ebEx", "In"="In")

## Create accessibility matrices with sorted Hi-C clusters
atac.region.HiC.counts <- lapply(atac.region.cl.counts, function(cl.mat) {
  cl.mat <- cl.mat[,names(atac.HiC.cl.map)]
  t(apply(cl.mat, 1, function(x) tapply(x, atac.HiC.cl.map, max)))
})

## Select accessibility cutoff
atac.region.cl.counts.vec <- unlist(lapply(atac.region.cl.counts, as.vector))
hist(atac.region.cl.counts.vec)

min.access <- 0.01
sum(atac.region.cl.counts.vec > min.access)/length(atac.region.cl.counts.vec)

## Find genes linked to peaks that overlap HARs/HGEs and that are accessible in the relevant cell type
atac.region.gene.contacts <- lapply(names(atac.peak.maps), function(h) {
  peaks <- atac.peak.maps[[h]]
  har.peaks <- atac.peaks[peaks]
  df.list <- lapply(names(conns.list), function(cl) {
    cl.conns <- conns.list[[cl]]
    
    ## Get links between regions and genes
    df <- RegionGeneLinks(har.peaks, cl.conns, link.promoter = T, promoter.region = c(-3e3, 3e3),
                          region.name = NULL, weight.col = NULL)
    
    ## Filter for regions that are accessible in the given cell type
    h.access <- atac.region.HiC.counts[[h]]
    access.regions <- unique(df$region)
    access.regions <- access.regions[access.regions %in% rownames(h.access)]
    access.regions <- access.regions[h.access[access.regions, cl] > min.access]
    subset(df, region %in% access.regions)
  })
  names(df.list) <- names(conns.list)
  df.list
})
names(atac.region.gene.contacts) <- names(atac.peak.maps)


## Get scaled single cell RNA-seq dataset
all.contact.genes <- unlist(lapply(atac.region.gene.contacts, function(df.list) {
  unlist(lapply(df.list, function(df) df$gene), F, F)
}), F, F)
all.contact.genes <- unique(all.contact.genes[all.contact.genes %in% rownames(atac.rna)])
length(all.contact.genes)

all.tfs <- unlist(lapply(rownames(atac.HAR.tf.z.expr[[1]]), function(x) strsplit(x, split = "::")[[1]]), F, F)
all.tfs <- unique(all.tfs[all.tfs %in% rownames(atac.rna)])
length(all.tfs)

atac.rna <- ScaleData(atac.rna, features = c(all.contact.genes, all.tfs))
atac.rna.scale <- GetAssayData(atac.rna, slot = "scale.data")
atac.rna.cl.scale <- sapply(levels(atac.rna.ident), function(cl) {
  Matrix::rowMeans(atac.rna.scale[,names(atac.rna.ident[atac.rna.ident == cl])])
})

## Add HAR/HGE gene correlations
atac.region.gene.contacts <- lapply(names(atac.peak.maps), function(h) {
  df.list <- atac.region.gene.contacts[[h]]
  lapply(df.list, function(df) {
    ## Add correlations between region accessibility and gene expression
    df <- AddLinkCorrelations(df, "region", "gene", 
                              atac.region.cl.counts[[h]][,atac.cl.order], 
                              atac.rna.cl.scale[,atac.cl.order])
  })
})
names(atac.region.gene.contacts) <- names(atac.peak.maps)

sapply(atac.region.gene.contacts, function(cl.contacts) {
  sapply(cl.contacts, nrow)
})

## Save results
# rm(atac.rna); invisible(gc());
save.image(output.file)



#### Link ATAC TFs to HARs/HGEs using motifs and correlations ####

## Create binned counts for improved correlations
atac.bin.counts <- GetBinCounts(atac.cell.counts, atac.umap.emb)
atac.bin.norm.counts <- RunTFIDF(atac.bin.counts)

## Run chromVAR at the single cell level for correlating TF activity to HAR/HGE peaks with TF motifs
atac.bin.SE <- SummarizedExperiment(assays = list(counts = atac.bin.counts),
                                    rowData = peak2granges(rownames(atac.bin.counts)), 
                                    colData = DataFrame(names = colnames(atac.bin.counts)))
atac.bin.SE <- addGCBias(atac.bin.SE, genome = BSgenome.Hsapiens.UCSC.hg38)

register(MulticoreParam(6))
atac.tf.bin.HAR <- RegionChromVAR(atac.bin.SE, atac_motif_ix, har.list, n.runs = 3)
atac.tf.bin.HAR <- lapply(atac.tf.bin.HAR, function(z) {
  z[is.na(z) | is.nan(z) | is.infinite(z)] <- 0
  rownames(z) <- sapply(rownames(z), ExtractField, 2)
  z
})


## Subset to peaks overlapping HAR/HGE regions to save memory
atac.mapped.peaks <- unique(unlist(atac.peak.maps, F, F))
atac.bin.counts <- atac.bin.counts[atac.mapped.peaks,]
atac.bin.norm.counts <- atac.bin.norm.counts[atac.mapped.peaks,]

## Format TF motif names
atac_motif_ix_mat <- assay(atac_motif_ix)
colnames(atac_motif_ix_mat) <- sapply(colnames(atac_motif_ix_mat), ExtractField, field = 2)

# atac.tf.npeaks <- sapply(1:ncol(atac_motif_ix_mat), function(i) sum(atac_motif_ix_mat[,i]))
# hist(atac.tf.npeaks)
# median(atac.tf.npeaks)

## Link TFs to regions using motifs
atac.TF.all.regions <- lapply(names(atac.peak.maps), function(h) {
  peaks <- atac.peak.maps[[h]]
  tf.peak.df <- TFRegionLinks(atac_motif_ix_mat[peaks,])
  AddLinkCorrelations(tf.peak.df, "TF", "region",
                      atac.tf.bin.HAR[[h]][,colnames(atac.bin.counts)],
                      atac.bin.norm.counts[peaks, colnames(atac.bin.counts)],
                      n.cores = 2)
})
names(atac.TF.all.regions) <- names(atac.peak.maps)
sapply(atac.TF.all.regions, nrow)

hist(unlist(lapply(atac.TF.all.regions, function(df) df$cor), F, F))
min.tf.peak.cor <- 0.1

atac.TF.regions <- lapply(atac.TF.all.regions, function(df) {
  subset(df, cor > min.tf.peak.cor)
})
sapply(atac.TF.regions, nrow)


## Save results
save.image(output.file)




#### Build ATAC HAR/HGE TF networks ####

## Get scaled single cell TF expression
atac.tf.scale <- lapply(rownames(atac.HAR.tf.z[[1]]), function(tf) {
  tfs <- strsplit(tf, split = "::")[[1]]
  tfs <- tfs[!grepl("var", tfs)]
  if (all(tfs %in% rownames(atac.rna.scale))) {
    if (length(tfs) > 1) return(colMeans(as.matrix(atac.rna.scale[tfs,])))
    else return(atac.rna.scale[tfs,])
  }
  else {
    return(NULL)
  }
})
names(atac.tf.scale) <- rownames(atac.HAR.tf.z[[1]])
atac.tf.scale <- do.call(rbind, atac.tf.scale)

## Link TFs to genes using cell type specific Hi-C
min.peak.gene.cor <- 0
atac.HiC.tf.expr <- t(apply(atac.tf.expr, 1, function(x) tapply(x[names(atac.HiC.cl.map)], atac.HiC.cl.map, max)))
atac.TF.cl.genes <- lapply(names(atac.region.gene.contacts), function(h) {
  peak.gene.list <- atac.region.gene.contacts[[h]]
  tf.peak.df <- atac.TF.regions[[h]]
  
  tf.gene.list <- lapply(names(peak.gene.list), function(cl) {
    peak.gene.df <- subset(peak.gene.list[[cl]], cor > min.peak.gene.cor & !is.na(cor))
    tf.genes <- TFGeneLinks(tf.peak.df, peak.gene.df)
    
    cl.expr <- atac.HiC.tf.expr[,cl]
    expr.tfs <- names(cl.expr[cl.expr >= min.expr])
    tf.genes <- subset(tf.genes, TF %in% expr.tfs)
    if (nrow(tf.genes) > 0) {
      tf.genes <- AddLinkCorrelations(tf.genes, "TF", "gene", 
                                      atac.tf.scale[,names(atac.rna.ident)], 
                                      atac.rna.scale[,names(atac.rna.ident)],
                                      n.cores = 2)
      tf.genes <- subset(tf.genes, !is.na(cor))
    } else {
      tf.genes <- NULL
    }
    tf.genes
  })
  names(tf.gene.list) <- names(peak.gene.list)
  tf.gene.list
})
names(atac.TF.cl.genes) <- names(atac.region.gene.contacts)
sapply(atac.TF.cl.genes, function(df.list) sapply(df.list, nrow))

# ## Export results
# save(atac.TF.cl.genes, file = tf.gene.contacts.file)

## Save full image
save.image(output.file)




#### Analyze ATAC HAR/HGE TF networks ####

## Subset to HAR TFs that are active/expressed in relevant cell types
min.tf.z <- 4
atac.TF.cl.graphs <- lapply(names(atac.TF.cl.genes), function(h) {
  z <- atac.HAR.tf.z.expr[[h]][,names(atac.HiC.cl.map)]
  # z[is.na(z)] <- -Inf
  z <- t(apply(z, 1, function(x) tapply(x, atac.HiC.cl.map, max)))
  
  cl.graphs <- lapply(names(atac.TF.cl.genes[[h]]), function(cl) {
    z.cl <- z[,cl]
    har.tfs <- names(z.cl[z.cl > min.tf.z])
    df <- atac.TF.cl.genes[[h]][[cl]]
    df$edge.label <- paste0(df$TF, "_", df$gene)
    
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
  names(cl.graphs) <- names(atac.TF.cl.genes[[h]])
  cl.graphs
})
names(atac.TF.cl.graphs) <- names(atac.TF.cl.genes)


atac.TF.cl.layouts <- lapply(atac.TF.cl.graphs, function(g.list) {
  lapply(g.list, function(g) {
    # layout_with_fr(g, weights = abs(E(g)$cor))
    layout_with_fr(g)
  })
})


## Plot cell type specific graphs
h.type <- "HGE_fetal"
h.graph.list <- atac.TF.cl.graphs[[h.type]]
h.layout.list <- atac.TF.cl.layouts[[h.type]]

pdf(paste0("fCTX-combined_ATAC_HAR_TF_gene_", h.type, "_network_nolabels.pdf"), width = 6.5, height = 5)
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
save.image(output.file)



#### Compare ATAC and THS TF network nodes ####

## Load THS TF gene graphs
fetal.TF.cl.genes <- readRDS("../Results/scTHS_TF_gene_network_dfs.Robj")
fetal.TF.cl.graphs <- readRDS("../Results/scTHS_TF_gene_igraphs.Robj")

## Compute centrality
ths.cl.centrality <- lapply(fetal.TF.cl.graphs, function(g.list) {
  lapply(g.list, function(g) {
    # gene.centrality <- page_rank(g, directed = F)$vector
    gene.centrality <- eigen_centrality(g, directed = F)$vector
    # gene.centrality <- degree(g, mode = "out")
    sort(gene.centrality, decreasing = T)
  })
})

atac.cl.centrality <- lapply(atac.TF.cl.graphs, function(g.list) {
  lapply(g.list, function(g) {
    # gene.centrality <- page_rank(g, directed = F)$vector
    gene.centrality <- eigen_centrality(g, directed = F)$vector
    # gene.centrality <- degree(g, mode = "out")
    sort(gene.centrality, decreasing = T)
  })
})


## Correlate gene centralities
h.type <- "HGE_fetal"
h.ths.cl.centrality <- ths.cl.centrality[[h.type]]
h.atac.cl.centrality <- atac.cl.centrality[[h.type]]

h.cl.corr.graphs <- lapply(names(h.ths.cl.centrality), function(cl) {
  h.ths.centrality <- h.ths.cl.centrality[[cl]]
  h.atac.centrality <- h.atac.cl.centrality[[cl]]
  
  genes.use <- intersect(names(h.ths.centrality), names(h.atac.centrality))
  genes.use <- genes.use[!(h.ths.centrality[genes.use] == 0 & h.atac.centrality[genes.use] == 0)]
  h.ths.centrality <- h.ths.centrality[genes.use]
  h.atac.centrality <- h.atac.centrality[genes.use]
  
  pts.label <- names(c(head(sort(h.ths.centrality, decreasing = T), n = 4), 
                       head(sort(h.atac.centrality, decreasing = T), n = 4)))
  pts.label <- genes.use %in% pts.label
  # pts.label <- h.ths.centrality[genes.use] > 0.25 | h.atac.centrality[genes.use] > 0.25
  PlotCorrelation(h.ths.centrality, h.atac.centrality, box = F,
                  use.label = T, show.corr = T, pts.label = pts.label, title = cl,
                  label.font.size = 3, pt.size = 1) + 
    theme(title = element_text(size = 11))
})
names(h.cl.corr.graphs) <- names(h.ths.cl.centrality)

pdf("fCTX-combined_ATAC_THS_Fetal_HGE_graph_centrality_corr.pdf", width = 6, height = 3)
plot_grid(plotlist = h.cl.corr.graphs[c("RG", "ipEx")], ncol = 2)
dev.off()


## Plot number of overlapping genes for each cell type specific network
h.cl.nGenes <- lapply(names(h.ths.cl.centrality), function(cl) {
  h.ths.centrality <- h.ths.cl.centrality[[cl]]
  h.atac.centrality <- h.atac.cl.centrality[[cl]]
  
  genes.use <- intersect(names(h.ths.centrality), names(h.atac.centrality))
  nGenes = c(length(h.atac.centrality), length(h.ths.centrality), length(genes.use))
  dataset = c("sciATAC", "scTHS", "Intersect")
  data.frame(dataset, cl = cl, nGenes)
})
h.cl.nGenes <- do.call(rbind, h.cl.nGenes)

nGenes.plot <- ggplot(data = h.cl.nGenes, aes(x = cl, y = nGenes, fill = dataset)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none") + coord_flip()

## Plot number of overlapping edges for each cell type specific network
h.ths.TF.genes <- fetal.TF.cl.genes[[h.type]]
h.atac.TF.genes <- atac.TF.cl.genes[[h.type]]
h.cl.nEdges <- lapply(names(h.ths.TF.genes), function(cl) {
  h.ths.df <- h.ths.TF.genes[[cl]]
  h.atac.df <- h.atac.TF.genes[[cl]]
  ths.edges <- unique(paste0(h.ths.df$TF, "_", h.ths.df$gene))
  atac.edges <- unique(paste0(h.atac.df$TF, "_", h.atac.df$gene))
  
  common.edges <- intersect(ths.edges, atac.edges)
  nEdges = c(length(atac.edges), length(ths.edges), length(common.edges))
  dataset = c("sciATAC", "scTHS", "Intersect")
  data.frame(dataset, cl = cl, nEdges)
})
h.cl.nEdges <- do.call(rbind, h.cl.nEdges)

nEdges.plot <- ggplot(data = h.cl.nEdges, aes(x = cl, y = nEdges, fill = dataset)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
        axis.text.y = element_blank(),
        legend.title = element_blank()) + coord_flip()

pdf("fCTX-combined_ATAC_THS_network_gene_edge_intersect_barplot.pdf", width = 4, height = 3.5)
plot_grid(nGenes.plot, nEdges.plot, align = "h", rel_widths = c(0.4, 0.55))
dev.off()

h.gene.frac <- tapply(h.cl.nGenes$nGenes, h.cl.nGenes$cl, function(x) x[[3]]/x[[2]])
h.edge.frac <- tapply(h.cl.nEdges$nEdges, h.cl.nEdges$cl, function(x) x[[3]]/x[[2]])
h.frac.df <- data.frame(cl = factor(c(names(h.gene.frac), names(h.edge.frac)), levels = names(h.gene.frac)), 
                        frac = c(h.gene.frac,h.edge.frac),
                        frac.type = c(rep("Gene", 4), rep("Edge", 4)))

# pdf("fCTX-combined_ATAC_THS_network_gene_edge_intersect_frac.pdf", width = 2.5, height = 3.5)
ggplot(data = h.frac.df, aes(x = cl, y = frac, fill = frac.type)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_blank()) + coord_flip()
# dev.off()

## Save results
save.image(output.file)



#### Compare ATAC and THS TF network edges ####

## Compute edge betweenness
ths.cl.edges <- lapply(fetal.TF.cl.graphs, function(g.list) {
  lapply(g.list, function(g) {
    w <- edge_betweenness(g, directed = F)
    names(w) <- E(g)$edge.label
    w
  })
})

atac.cl.edges <- lapply(atac.TF.cl.graphs, function(g.list) {
  lapply(g.list, function(g) {
    w <- edge_betweenness(g, directed = F)
    names(w) <- E(g)$edge.label
    w
  })
})


## Correlate gene centralities
h.type <- "HGE_fetal"
h.ths.cl.edges <- ths.cl.edges[[h.type]]
h.atac.cl.edges <- atac.cl.edges[[h.type]]

h.edges.corr.graphs <- lapply(names(h.ths.cl.edges), function(cl) {
  ths.edges <- h.ths.cl.edges[[cl]]
  atac.edges <- h.atac.cl.edges[[cl]]
  
  edges.use <- intersect(names(ths.edges), names(atac.edges))
  edges.use <- edges.use[!(ths.edges[edges.use] == 0 & atac.edges[edges.use] == 0)]
  ths.edges <- ths.edges[edges.use]
  atac.edges <- atac.edges[edges.use]
  
  pts.label <- names(c(head(sort(ths.edges, decreasing = T), n = 4), 
                       head(sort(atac.edges, decreasing = T), n = 4)))
  pts.label <- edges.use %in% pts.label
  # pts.label <- ths.edges[edges.use] > 0.25 | atac.edges[edges.use] > 0.25
  PlotCorrelation(ths.edges, atac.edges, box = F,
                  use.label = T, show.corr = T, pts.label = pts.label, title = cl,
                  label.font.size = 3, pt.size = 1) + 
    theme(title = element_text(size = 11))
})
names(h.edges.corr.graphs) <- names(h.ths.cl.edges)

pdf("fCTX-combined_ATAC_THS_Fetal_HGE_graph_edge_betweenness_corr.pdf", width = 6, height = 3)
plot_grid(plotlist = h.edges.corr.graphs[c("RG", "ipEx")], ncol = 2)
dev.off()



#### Validate ATAC HAR/HGE-linked genes using human/chimp gene logFC in organoids ####

## Load human/chimp DEG data
load("../Data/human_chimp_RNA_logFC.RData")
colnames(human_chimp_rna_logFC) <- plyr::revalue(colnames(human_chimp_rna_logFC), 
                                                 replace = c("IPC/eN"="ipEx", "eN" = "ebEx"))

## Run GSEA on HAR-linked genes
gsea.list <- lapply(names(har.list), function(h.type) {
  TF.cl.target.genes <- lapply(c(fetal.TF.cl.genes[[h.type]], atac.TF.cl.genes[[h.type]]), function(df) {
    # df <- subset(df, abs(cor) > 0.01)
    df <- subset(df, gene %in% rownames(human_chimp_rna_logFC))
    unique(df$gene)
  })
  names(TF.cl.target.genes) <- c(paste0(names(fetal.TF.cl.genes[[h.type]]), "-ths"),
                                 paste0(names(atac.TF.cl.genes[[h.type]]), "-atac"))
  sapply(TF.cl.target.genes, length)
  
  res <- lapply(colnames(human_chimp_rna_logFC), function(cl) {
    df <- bulk.gsea(human_chimp_rna_logFC[,cl], set.list = TF.cl.target.genes,
                    mc.cores = 8, rank = F)
    df$org.cluster <- cl
    df$graph.cluster <- rownames(df)
    df
  })
  gsea.df <- do.call(rbind, res); rownames(gsea.df) <- NULL;
  gsea.df$lp <- -log(gsea.df$p.val) * sign(gsea.df$sscore)
  # gsea.df$lp[gsea.df$lp < 0] <- 0
  gsea.df$dataset <- sapply(gsea.df$graph.cluster, function(x) strsplit(x, split = "-")[[1]][[2]])
  gsea.df$graph.cluster <- sapply(gsea.df$graph.cluster, function(x) strsplit(x, split = "-")[[1]][[1]])
  subset(gsea.df, org.cluster == graph.cluster)
})
names(gsea.list) <- names(har.list)

fetal.hge.plot <- ggplot(data = gsea.list$HGE_fetal, aes(x = graph.cluster, y = lp, fill = dataset)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, color = "black", angle = 90),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none") + ylim(-4, 10)

encode.ctrl.plot <- ggplot(data = gsea.list$ENCODE_DNAse, aes(x = graph.cluster, y = lp, fill = dataset)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, color = "black", angle = 90),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_blank()) + ylim(-4, 10)



pdf("fCTX-combined_ATAC_HAR_TF_human_chimp_gsea.pdf", width = 4, height = 3)
plot_grid(fetal.hge.plot, encode.ctrl.plot, rel_widths = c(0.4125, 0.6))
dev.off()

## Save results
save.image(output.file)



#### Run pathway analysis on NFIC targets ####

## Look at number of ATAC & THS binding sites
tf <- "NFIC"
h.types <- c("HGE_fetal")

## Count TF binding sites
tf.nbind <- subset(atac.TF.regions[[h.type]], TF == tf)
nrow(tf.nbind)

## Look at intersection of ATAC & THS targets
tf <- "NFIC"
h.types <- c("HGE_fetal")

tf.nGenes.plots <- lapply(h.types, function(h.type) {
  tf.ths.genes <- lapply(fetal.TF.cl.genes[[h.type]], function(df) subset(df, TF == tf)$gene)
  tf.atac.genes <- lapply(atac.TF.cl.genes[[h.type]], function(df) subset(df, TF == tf)$gene)
  
  tf.nGenes <- lapply(c("RG", "ipEx"), function(cl) {
    ths.genes <- tf.ths.genes[[cl]]
    atac.genes <- tf.atac.genes[[cl]]
    genes.ix <- intersect(ths.genes, atac.genes)
    data.frame(nGenes = c(nTHS = length(ths.genes), nATAC = length(atac.genes), nIntersect = length(genes.ix)),
               dataset = c("scTHS", "sciATAC", "Intersect"), cl = cl)
  })
  tf.nGenes <- do.call(rbind, tf.nGenes)
  
  ggplot(data = tf.nGenes, aes(x = cl, y = nGenes, fill = dataset)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.title = element_blank()) + ylim(0,500)
})
# tf.nGenes.plots[[1]] <- tf.nGenes.plots[[1]] + theme(legend.position = "none")

pdf("fCTX-combined_NFIC_target_nGenes.pdf", width = 3.5, height = 1.75)
# plot_grid(plotlist = tf.nGenes.plots, rel_widths = c(0.4, 0.6))
tf.nGenes.plots[[1]]
dev.off()


## Run GSEA on TF targets

## Load genesets
go.bp.genesets <- ReadGenesets("../Data/Genesets/go.bp.v7.2.gmt")
pathway.genesets <- ReadGenesets("../Data/Genesets/pathway.v7.2.gmt")
hallmark.genesets <- ReadGenesets("../Data/Genesets/hallmark.v7.2.gmt")

## Merge genesets
ctx.genes <- unique(c(rownames(atac.rna.cl.mean), rownames(fbcx)))
genesets <- c(hallmark.genesets, pathway.genesets)
genesets <- FilterGenesets(genesets, ctx.genes, min.size = 10, max.size = 1000)
names(genesets) <- sapply(names(genesets), function(x) gsub("KEGG_|REACTOME_|HALLMARK_", "", x))
names(genesets) <- sapply(names(genesets), function(x) gsub("_", " ", x))
names(genesets) <- sapply(names(genesets), function(x) simpleCap(tolower(x)))
length(genesets)

h.type <- "HGE_fetal"
tf.ths.genes <- lapply(fetal.TF.cl.genes[[h.type]], function(df) subset(df, TF == tf)$gene)
tf.atac.genes <- lapply(atac.TF.cl.genes[[h.type]], function(df) subset(df, TF == tf)$gene)

tf.ths.gsea <- lapply(tf.ths.genes, function(genes) {
  gsea.df <- FisherEnrich(genes, genesets, n.background = length(ctx.genes), correct = T)
  gsea.df[order(gsea.df$p.val),]
})

head(tf.ths.gsea$RG, n = 8)
head(tf.ths.gsea$ipEx, n = 8)

tf.atac.gsea <- lapply(tf.atac.genes, function(genes) {
  gsea.df <- FisherEnrich(genes, genesets, n.background = length(ctx.genes), correct = T)
  gsea.df[order(gsea.df$p.val),]
})

head(tf.atac.gsea$RG, n = 8)
head(tf.atac.gsea$ipEx, n = 8)

## Plot heatmap of geneset enrichment
genesets.ix <- tf.ths.gsea$RG$genesets
gsea.mat <- lapply(c(tf.ths.gsea[c("RG", "ipEx")], tf.atac.gsea[c("RG", "ipEx")]), function(df) {
  lp <- -log10(df$p.val)
  names(lp) <- df$genesets
  lp[genesets.ix]
})
gsea.mat <- do.call(cbind, gsea.mat)

top.genesets <- unlist(lapply(1:ncol(gsea.mat), function(i) {
  lp <- gsea.mat[,i]
  names(head(sort(lp, decreasing = T), n = 5))
}), F, F)
top.genesets <- unique(top.genesets)

pdf("fCTX-combined_NFIC_target_gsea.pdf", width = 5, height = 5.25)
ggHeat(gsea.mat[top.genesets,])
dev.off()

# ## Plot correlation of geneset enrichment
# top.n.plot <- 8
# use.label <- T
# tf.gsea.cor.plots <- lapply(c("RG", "ipEx"), function(cl) {
#   ths.gsea <- tf.ths.gsea[[cl]]
#   atac.gsea <- tf.atac.gsea[[cl]]
#   
#   ths.lp <- -log(ths.gsea$p.val)
#   names(ths.lp) <- ths.gsea$genesets
#   atac.lp <- -log(atac.gsea$p.val)
#   names(atac.lp) <- atac.gsea$genesets
#   
#   ix <- intersect(names(ths.lp), names(atac.lp))
#   ths.lp <- ths.lp[ix]; atac.lp <- atac.lp[ix];
#   
#   top.genesets <- c(head(names(sort(ths.lp, decreasing = T)), n = top.n.plot),
#                     head(names(sort(atac.lp, decreasing = T)), n = top.n.plot))
#   PlotCorrelation(ths.lp, atac.lp, box = F, use.label = use.label, pt.size = 0.75, show.corr = T,
#                   pts.label = ix %in% top.genesets, label.font.size = 2, title = cl)  
# })
# 
# # pdf("fCTX-combined_NFIC_target_gsea_cor_nolabels.pdf", width = 6, height = 3)
# plot_grid(plotlist = tf.gsea.cor.plots)
# # dev.off()


## Look at Notch signaling
notch.genes <- unique(unlist(genesets[c("Signaling By Notch", "Signaling By Notch3")], F, F))
notch.genes <- intersect(notch.genes, c(tf.ths.genes$RG, tf.atac.genes$RG))
notch.genes

ggHeat(atac.rna.cl.scale[notch.genes, atac.cl.order], clustering = "row")

## Look at cell adhesion genes
adhesion.genesets <- c("Naba Matrisome", "Naba Core Matrisome", "Naba Ecm Glycoproteins")
adhesion.genes <- unique(unlist(genesets[adhesion.genesets], F, F))
adhesion.genes <- intersect(adhesion.genes, intersect(tf.ths.genes$RG, tf.atac.genes$RG))
adhesion.genes

ggHeat(atac.rna.cl.scale[adhesion.genes, atac.cl.order], clustering = "row")

## Save results
save.image(output.file)



#### Visualize NFIC Cell Adhesion Sub-Network ####

cls <- c("RG", "ipEx")
h.type <- "HGE_fetal"

genesets.merged.df <- lapply(cls, function(cl) {
  df1 <- fetal.TF.cl.genes[[h.type]][[cl]]
  rownames(df1) <- paste0(df1$TF, "_", df1$gene)
  
  df2 <- atac.TF.cl.genes[[h.type]][[cl]]
  rownames(df2) <- paste0(df2$TF, "_", df2$gene)
  
  df.common <- df1[intersect(rownames(df1), rownames(df2)), c("TF", "gene")]
  df.common$edge <- rownames(df.common)
  subset(df.common, gene %in% adhesion.genes)
  
})
# names(genesets.merged.df) <- cls
genesets.merged.df <- do.call(rbind, genesets.merged.df)
rownames(genesets.merged.df) <- NULL
genesets.merged.df <- dplyr::distinct(genesets.merged.df)
genesets.merged.df

## Make graph of cell adhesion target genes
adhesion.g <- graph_from_data_frame(genesets.merged.df, directed = TRUE, vertices = NULL)

## Set vertex attributes
V(adhesion.g)$color <- sapply(names(V(adhesion.g)), function(x) {
  if (x %in% genesets.merged.df$TF) return("tomato")
  else return("skyblue")
})
V(adhesion.g)$size <- sapply(names(V(adhesion.g)), function(x) {
  if (x %in% genesets.merged.df$TF) return(7)
  else return(3)
})

## Compute layout
adhesion.l <- layout_with_fr(adhesion.g)

pdf("fCTX-combined_NFIC_adhesion_subgraph.pdf", width = 5, height = 5)
plot.igraph(adhesion.g, layout = adhesion.l, vertex.size = V(adhesion.g)$size, 
            vertex.label = names(V(adhesion.g)), vertex.color = V(adhesion.g)$color, 
            edge.arrow.size = 0.4, edge.width = 1, 
            vertex.label.cex = 0.8)
dev.off()

pdf("fCTX-combined_NFIC_adhesion_subgraph_nolabels.pdf", width = 5, height = 5)
plot.igraph(adhesion.g, layout = adhesion.l, vertex.size = V(adhesion.g)$size, 
            vertex.label = NA, vertex.color = V(adhesion.g)$color, 
            edge.arrow.size = 0.4, edge.width = 1, 
            vertex.label.cex = 0.8)
dev.off()



#### Visualize LHX2 targets ####

cls <- c("RG", "ipEx")
h.type <- "HGE_fetal"

merged.df <- lapply(cls, function(cl) {
  df1 <- fetal.TF.cl.genes[[h.type]][[cl]]
  rownames(df1) <- paste0(df1$TF, "_", df1$gene)
  
  df2 <- atac.TF.cl.genes[[h.type]][[cl]]
  rownames(df2) <- paste0(df2$TF, "_", df2$gene)
  
  df.common <- df1[intersect(rownames(df1), rownames(df2)), c("TF", "gene")]
})
merged.df <- do.call(rbind, merged.df)
merged.df <- dplyr::distinct(merged.df)

## Get LHX2 targets
lhx2.targets <- subset(merged.df, TF == "LHX2")$gene
lhx2.merged.df <- subset(merged.df, gene %in% lhx2.targets)
lhx2.merged.df

## Merged graph
g <- graph_from_data_frame(lhx2.merged.df, directed = TRUE, vertices = NULL)

## Set vertex attributes
V(g)$color <- sapply(names(V(g)), function(x) {
  if (x %in% lhx2.merged.df$TF) return("tomato")
  else return("skyblue")
})
V(g)$size <- sapply(names(V(g)), function(x) {
  if (x %in% lhx2.merged.df$TF) return(7)
  else return(3)
})
V(g)$label <- names(V(g))

## Compute layout
l <- layout_with_fr(g)

pdf("fCTX-combined_LHX2_subgraph.pdf", width = 5, height = 5)
plot.igraph(g, layout = l, vertex.size = V(g)$size, vertex.label = V(g)$label,
            vertex.color = V(g)$color, edge.arrow.size = 0.4, edge.width = 1, 
            vertex.label.cex = 0.8)
dev.off()

pdf("fCTX-combined_LHX2_subgraph_nolabels.pdf", width = 5, height = 5)
plot.igraph(g, layout = l, vertex.size = V(g)$size, vertex.label = NA,
            vertex.color = V(g)$color, edge.arrow.size = 0.4, edge.width = 1, 
            vertex.label.cex = 0.8)
dev.off()


## Save results
save.image(output.file)



#### Plot accessibility tracks for NFIC/LHX2 HGEs ####

## Generate coaccessibility plot centered on HAR/HGE-linked gene
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(cicero)

## Get txDB
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## Load THS networks
fetal.TF.regions <- readRDS("../Results/scTHS_TF_region_dfs.Robj")
fetal.region.gene.contacts <- readRDS("../Results/scTHS_region_gene_dfs.Robj")

## Parameters
h.type <- "HGE_fetal"
cls <- c("RG", "ipEx")
genes <- c("C1orf61")
buffer <- 5e4

## Check which peaks are contacting the target genes
contacts <- do.call(rbind, fetal.region.gene.contacts[[h.type]][cls])
# contacts <- dplyr::distinct(contacts[,c("region", "gene")])
gene.regions <- subset(contacts, gene %in% genes)$region

## Check which TFs are linked to those peaks
HiC.cl.map <- c("oRG" = "RG", "vRG" = "RG", "ipEx" = "ipEx", "ebEx" = "ebEx")
expr.tfs <- t(apply(HAR.tf.z.expr[[h.type]][,names(HiC.cl.map)], 1, function(x) {
  x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
  x[x != 0] <- 1
  tapply(x, HiC.cl.map, max)
}))
cls.expr.tfs <- apply(expr.tfs[,cls], 1, max)
cls.expr.tfs <- names(cls.expr.tfs[cls.expr.tfs != 0])  

motif.regions <- fetal.TF.regions[[h.type]]
motif.regions <- subset(motif.regions, region %in% gene.regions & TF %in% cls.expr.tfs)
head(motif.regions)


## Get Hi-C contacts that overlap HAR/HGE regions
conns <- do.call(rbind, conns.list[cls])
conns.regions.ix <- findOverlaps(peak2granges(conns$Peak1), peak2granges(gene.regions))
conns <- conns[conns.regions.ix@from,]
head(conns)

# ## Load pre-computed gene promoters
# genes.proms <- readRDS("gene_promoters.Robj")
# genes.proms <- genes.proms[genes.proms$gene %in% genes]

## Get Hi-C contacts that overlap target gene promoters
target.peak.anno <- annotatePeak(peak2granges(conns$Peak2), tssRegion = c(-3000, 3000), level = "transcript",
                                 TxDb = txdb, annoDb = "org.Hs.eg.db")
target.peak.anno <- as.data.frame(target.peak.anno)

conns <- conns[target.peak.anno$SYMBOL %in% genes & grepl("Promoter", target.peak.anno$annotation),]
conns <- dplyr::distinct(conns)
conns$coaccess <- 1
head(conns)

## Get limits of visualization
viz.peaks <- peak2granges(unique(c(conns$Peak1, conns$Peak2)))
viz.chr <- levels(seqnames(viz.peaks))
viz.start <- start(viz.peaks)
viz.end <- end(viz.peaks)
viz.min <- min(c(viz.start, viz.end)) - buffer
viz.max <- max(c(viz.start, viz.end)) + buffer

## Make plot
plot_list <- plot_connections(conns, viz.chr, viz.min, viz.max,
                              coaccess_cutoff = 0,
                              alpha_by_coaccess = T,
                              include_axis_track = T,
                              return_as_list = T,
                              connection_width = 1)

## Load HAR track
plot_list[[2]] <- AnnotationTrack(har.list$HGE_fetal, chromosome = viz.chr)

## Load accessibility data
dTrack1 <- DataTrack(range = "../Data/scTHS_cluster_bigwigs/RG.merged.sorted.bw",
                     genome = "hg38", type = "l", chromosome = viz.chr,
                     name = "bigWig")
plot_list[[4]] <- dTrack1

dTrack2 <- DataTrack(range = "../Data/scTHS_cluster_bigwigs/ipEx.merged.sorted.bw",
                     genome = "hg38", type = "l", chromosome = viz.chr,
                     name = "bigWig")
plot_list[[5]] <- dTrack2

dTrack3 <- DataTrack(range = "../Data/scTHS_cluster_bigwigs/Ex.merged.sorted.bw",
                     genome = "hg38", type = "l", chromosome = viz.chr,
                     name = "bigWig")
plot_list[[6]] <- dTrack3

## Replace gene annotation track
plot_list[[7]] <- GeneRegionTrack(txdb, chromosome = viz.chr, start = viz.min, end = viz.max,
                                  collapseTranscripts = "meta", geneSymbols = T,
                                  transcriptAnnotation = "symbol")
symbol(plot_list[[7]]) <- mapIds(org.Hs.eg.db, gene(plot_list[[7]]), 'SYMBOL', 'ENTREZID')

pdf("fCTX-combined_LHX2_C1orf61_HGE_conns.pdf", width = 5, height = 4)
Gviz::plotTracks(plot_list,
                 sizes = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25),
                 from = viz.min, to = viz.max, chromosome = viz.chr,
                 transcriptAnnotation = "symbol",
                 col.axis = "black",
                 fontsize.group = 10,
                 fontcolor.legend = "black",
                 lwd = 0.3,
                 title.width = 0.75,
                 background.title = "transparent",
                 col.border.title = "transparent",
                 cex.main = 4)
dev.off()


## Get all HAR regions in gviz plot
viz.gr <- GRanges(seqnames = viz.chr, ranges = IRanges(start = viz.min, end = viz.max))
print(as.data.frame(sort(subsetByOverlaps(har.regions, viz.gr))))



#### Export supplementary tables ####

## Export HAR/HGE cell type enrichment
write.table(HAR.z, "../Results/HAR_HGE_cell_type_zscores.tsv", sep = "\t", quote = F)

## Export HAR/HGE TF motif enrichment
for (h in names(merged.HAR.tf.z.expr)) {
  h.tf.z <- merged.HAR.tf.z.expr[[h]]
  h.tf.z[is.na(h.tf.z)] <- 0
  
  fname <- paste0("../Results/", h, "_", "tf_expr_zscores.tsv")
  write.table(h.tf.z, file = fname, sep = "\t", quote = F)
}

## Map peaks to HARs/HGEs
atac.peak.region.mappings <- lapply(har.list, function(gr) {
  peak.region.list <- RegionOverlapList(atac.peaks, gr, gr1.name = "Peak")
  peak.region.list <- lapply(peak.region.list, function(gr) {
    paste0(granges2peak(gr), collapse = ",")
  })
  hash.map <- new.env(hash = T)
  for (p in names(peak.region.list)) {
    hash.map[[p]] <- peak.region.list[[p]]
  }
  hash.map
})

## Export ATAC HAR/HGE to gene links
for (h in names(atac.region.gene.contacts)) {
  region.mapping <- atac.peak.region.mappings[[h]]
  df.list <- atac.region.gene.contacts[[h]]
  
  df.list <- lapply(names(df.list), function(cl) {
    df <- df.list[[cl]]
    df$cell.type <- cl
    df$accelerated.region <- sapply(df$region, function(p) region.mapping[[p]])
    df
  })
  df <- do.call(rbind, df.list)
  df$atac.region <- df$region
  df <- df[,c("cell.type", "atac.region", "accelerated.region", "gene")]
  print(head(df))
  
  fname <- paste0("../Results/", h, "_", "ATAC_region_gene_network.tsv")
  write.table(df, file = fname, sep = "\t", row.names = F)
}

## Export ATAC TF to HAR/HGE links
for (h in names(atac.TF.regions)) {
  region.mapping <- atac.peak.region.mappings[[h]]
  df <- atac.TF.regions[[h]]
  df$accelerated.region <- sapply(df$region, function(p) region.mapping[[p]])
  df$atac.region <- df$region
  df <- df[,c("TF", "atac.region", "accelerated.region")]
  print(head(df))
  
  fname <- paste0("../Results/", h, "_", "ATAC_TF_region_network.tsv")
  write.table(df, file = fname, sep = "\t", row.names = F, quote = F)
}

## Export ATAC TF to gene links
for (h in names(atac.TF.cl.genes)) {
  region.mapping <- atac.peak.region.mappings[[h]]
  df.list <- atac.TF.cl.genes[[h]]
  df.list <- lapply(names(df.list), function(cl) {
    df <- df.list[[cl]]
    df$cell.type <- cl
    df
  })
  df <- do.call(rbind, df.list)
  df$atac.region <- df$region
  df <- df[,c("cell.type", "TF", "gene")]
  print(head(df))
  
  fname <- paste0("../Results/", h, "_", "ATAC_TF_gene_network.tsv")
  write.table(df, file = fname, sep = "\t", row.names = F, quote = F)
}

## Get network stats
h.type <- "HGE_fetal"

sapply(atac.TF.cl.genes[[h.type]], function(df) length(unique(df$TF)))
sapply(atac.TF.cl.genes[[h.type]], function(df) mean(table(df$TF)))
mean(table(atac.TF.regions[[h.type]]$TF))

sapply(atac.region.gene.contacts[[h.type]], function(df) length(unique(df$region)))
sapply(atac.region.gene.contacts[[h.type]], function(df) mean(table(df$region)))