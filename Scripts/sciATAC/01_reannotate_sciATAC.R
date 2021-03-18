library(swne)
library(Seurat)
library(Signac)
library(cisTopic)
library(chromfunks)
library(cellMapper)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)


## Cell type ordering for visualizations
dev.cl.order <- c("oRG", "vRG", "PgG2M", "PgS", "IP", 
                  "ExN", "ExM", "ExM-U", "ExDp1", "ExDp2", 
                  "InCGE", "InMGE", "OPC", "Mic", "End", "Per")
adult.cl.order <- c("Ast", "Ex", "In", "Opc", "Oli", "Mic", "End")
atac.cl.order <- c("RG", "ipEx", "RG/OPC", "ebEx", "Ex", "In", "End")


#### Load and normalize THS and ATAC datasets ####

## Load fetal cortex reference data
raw_counts_mat <- ReadData("../Data/ref_fetal_cortex_counts.tsv.gz")
raw_counts_mat <- FilterData(raw_counts_mat, min.samples.frac = 0.001, trim = 5, 
                             min.nonzero.features = 0, max.sample.sum = Inf)
dim(raw_counts_mat)

meta.data <- read.table("../Data/ref_fetal_cortex_metadata.csv.gz", header = T, row.names = 1, sep = ",")
ref.annot <- meta.data$Cluster; names(ref.annot) <- rownames(meta.data)
table(ref.annot)

dev <- CreateSeuratObject(raw_counts_mat)
dev <- NormalizeData(dev)
Idents(dev) <- ref.annot[colnames(dev)]

dev.markers <- FindAllMarkers(dev, logfc.threshold = 0.25, only.pos = T)
dev.markers.genes <- unique(as.character(dev.markers$gene))
length(dev.markers.genes)

VariableFeatures(dev) <- dev.markers.genes


## Load adult reference data
vctx <- readRDS("../Data/snDrop_RNA_seurat_object.Robj")
vctx <- FindVariableFeatures(vctx)

adult.rna.idents <- as.character(vctx$Identity)
names(adult.rna.idents) <- colnames(vctx)

adult.rna.idents[grepl("Ex", adult.rna.idents)] <- "Ex"
adult.rna.idents[grepl("In", adult.rna.idents)] <- "In"
adult.rna.idents[adult.rna.idents == "OPC"] <- "Opc"
adult.rna.idents <- factor(adult.rna.idents)
table(adult.rna.idents)

## Load ATAC data
atac <- readRDS("../Data/scTHS_chromatin_seurat_object.Robj")
atac.ident <- atac$cell_type

DefaultAssay(atac) <- "RNA"
atac <- NormalizeData(atac, assay = "RNA")
atac <- FindVariableFeatures(atac, nfeatures = 3000)

atac.umap.emb <- cbind(atac$tissue_umap_1, atac$tissue_umap_2)
colnames(atac.umap.emb) <- c("umap1", "umap2")
atac[["umap"]] <- CreateDimReducObject(atac.umap.emb, key = "UMAP_", assay = "RNA")
atac[["umap"]] <- CreateDimReducObject(atac.umap.emb, key = "UMAP_", assay = "peaks")




#### Map ATAC clusters to THS clusters ####

## Subcluster ATAC Excitatory neurons
table(atac$cell_type)
atac.ex <- atac[,atac$cell_type == "Excitatory neurons"]

DefaultAssay(atac.ex) <- "peaks"
# atac.ex <- RunTFIDF(atac.ex)
# atac.ex <- FindTopFeatures(atac.ex, min.cutoff = 'q0')
# atac.ex <- RunSVD(atac.ex)
# 
# DepthCor(atac.ex)
# dims.use <- 2:30
# 
# atac.ex <- FindNeighbors(atac.ex, reduction = 'lsi', dims = dims.use)
atac.ex <- FindNeighbors(atac.ex, reduction = 'umap', dims = 1:2)
atac.ex <- FindClusters(atac.ex, resolution = 0.1)

## Plot results
DimPlot(object = atac.ex, label = TRUE) + NoLegend() + theme_void()
table(Idents(atac.ex))

## Add finer resolution excitatory neuron clusters to annotations
atac.ident <- as.character(atac$cell_type); names(atac.ident) <- colnames(atac);
atac.ident[names(Idents(atac.ex))] <- paste0("Ex", Idents(atac.ex))
atac.ident <- factor(atac.ident)
table(atac.ident)

## Identify genes to use for comparison
genes.map <- unique(c(VariableFeatures(atac), VariableFeatures(dev), VariableFeatures(vctx)))
length(genes.map)

dev <- ScaleData(dev, features = genes.map)
vctx <- ScaleData(vctx, features = genes.map)
atac <- ScaleData(atac, assay = "RNA", features = genes.map)

atac.dev.cor <- MapClustersCor(query.data = GetAssayData(atac, assay = "RNA", slot = "scale.data"),
                               query.clusters = atac.ident,
                               train.data = GetAssayData(dev, slot = "scale.data"),
                               train.clusters = Idents(dev))

atac.adult.cor <- MapClustersCor(query.data = GetAssayData(atac, assay = "RNA", slot = "scale.data"),
                                 query.clusters = atac.ident,
                                 train.data = GetAssayData(vctx, slot = "scale.data"),
                                 train.clusters = adult.rna.idents)

atac.dev.cor.mat <- atac.dev.cor$cluster.correlations[dev.cl.order,]
atac.adult.cor.mat <- atac.adult.cor$cluster.correlations[adult.cl.order,]
merged.cor.mat <- rbind(atac.dev.cor.mat, atac.adult.cor.mat)

pdf("fCTX-combined_ATAC_orig_ident_dev_cortex_cor.pdf", width = 6.75, height = 4.75)
ggHeat(t(merged.cor.mat))
dev.off()


## Visualize key marker gene expression
marker.genes <- c("GLI3", "VIM", "NES", "FABP7", "SOX2", ## RG markers
                  "EOMES", "PPP1R17", "NEUROG2", "PAX6", ## ipEx markers,
                  "S100B", "NEU1" ## Ast markers
)

## Make heatmap
atac <- ScaleData(atac, features = marker.genes, assay = "RNA")
atac.marker.heat <- t(apply(GetAssayData(atac, assay = "RNA", slot = "scale.data"), 1, function(x) {
  tapply(x, atac.ident, mean)
}))

pdf("fCTX-combined_ATAC_orig_ident_marker_heatmap.pdf", width = 5, height = 4.75)
ggHeat(t(atac.marker.heat[marker.genes,]))
dev.off()

## Rename ATAC idents to match better with THS idents
atac.ths.map <- c("Astrocytes"="RG", 
                  "Astrocytes/Oligodendrocytes"="RG/OPC",
                  "Cerebrum_Unknown.3"="ipEx", 
                  "Ex0"="ebEx",
                  "Ex1"="Ex",
                  "Ex2"="Ex",
                  "Ex3"="Ex",
                  "Ex4"="Ex",
                  "Ex5"="Ex",
                  "Ex6"="ebEx",
                  "Ex7"="Ex",
                  "Ex8"="Ex",
                  "Ex9"="Ex",
                  "Inhibitory neurons"="In", 
                  "Limbic system neurons"="LimbicNeurons",
                  "SKOR2_NPSR1 positive cells"="SKOR2_NPSR1_Cells",
                  "Vascular endothelial cells"="End")

## Rename ATAC idents
table(atac.ident)
atac.mapped.ident <- plyr::revalue(atac.ident, replace = atac.ths.map)
table(atac.mapped.ident)

## Plot UMAP embedding
atac <- SetIdent(atac, value = atac.mapped.ident)
cells.plot <- names(atac.mapped.ident[atac.mapped.ident %in% atac.cl.order])

mapped.colors.df <- read.table("../Data/harmonized_cluster_colors.txt", header = T, sep = "\t",
                               stringsAsFactors = F, comment.char = "")
mapped.colors <- mapped.colors.df$color
names(mapped.colors) <- mapped.colors.df$cluster
mapped.colors[["Ex"]] <- mapped.colors[["ExL2/3"]]
mapped.colors[["In"]] <- mapped.colors[["InCGE"]]
mapped.colors[["RG"]] <- mapped.colors[["oRG"]]
mapped.colors[["RG/OPC"]] <- mapped.colors[["OPC"]]

pdf("fCTX-combined_ATAC_renamed_ident_umap.pdf", width = 4, height = 4)
cells.sample <- sample(cells.plot, size = 10000)
PlotDims(atac.umap.emb[cells.sample,], sample.groups = atac.mapped.ident[cells.sample], 
         pt.size = 0.3, alpha = 0.5, do.label = T, show.legend = F,
         colors.use = mapped.colors, show.axes = F, label.size = 5)
dev.off()


## Plot correlation using mapped idents
atac <- ScaleData(atac, assay = "RNA", features = genes.map)
mapped.dev.cor <- MapClustersCor(query.data = GetAssayData(atac, assay = "RNA", slot = "scale.data"),
                                 query.clusters = atac.mapped.ident,
                                 train.data = GetAssayData(dev, slot = "scale.data"),
                                 train.clusters = Idents(dev))

mapped.adult.cor <- MapClustersCor(query.data = GetAssayData(atac, assay = "RNA", slot = "scale.data"),
                                   query.clusters = atac.mapped.ident,
                                   train.data = GetAssayData(vctx, slot = "scale.data"),
                                   train.clusters = adult.rna.idents)

mapped.dev.cor.mat <- mapped.dev.cor$cluster.correlations[dev.cl.order,]
mapped.adult.cor.mat <- mapped.adult.cor$cluster.correlations[adult.cl.order,]
mapped.cor.mat <- rbind(mapped.dev.cor.mat, mapped.adult.cor.mat)

pdf("fCTX-combined_ATAC_GA_correlations.pdf", width = 6.75, height = 3)
ggHeat(t(mapped.cor.mat[, rev(atac.cl.order)]))
dev.off()


## Plot marker expression for mapped idents
atac <- ScaleData(atac, features = marker.genes, assay = "RNA")
atac.marker.heat <- t(apply(GetAssayData(atac, assay = "RNA", slot = "scale.data"), 1, function(x) {
  tapply(x, atac.mapped.ident, mean)
}))

pdf("fCTX-combined_ATAC_marker_heatmap.pdf", width = 4.75, height = 3)
ggHeat(t(atac.marker.heat[marker.genes, rev(atac.cl.order)]))
dev.off()

# ## Subset to only cell types overlapping THS cell types
# atac.mapped.ident <- droplevels(atac.mapped.ident[atac.mapped.ident %in% atac.cl.order])
# table(atac.mapped.ident)
# 
# atac <- atac[,names(atac.mapped.ident)]
# dim(atac)


#### LiftOver ATAC Data to hg38 ####

## Set assay
DefaultAssay(atac) <- "peaks"

## Filter peaks
atac.cell.counts <- GetAssayData(atac, assay = "peaks")
dim(atac.cell.counts)

## Get peak granges in hg19
atac.peaks.hg19 <- peak2granges(rownames(atac.cell.counts), delim = c("-", "-"))

## LiftOver to hg38
ch <- import.chain("GWAS/Chains/hg19ToHg38.over.chain")
atac.peaks.map <- liftOver(atac.peaks.hg19, chain = ch)
names(atac.peaks.map) <- granges2peak(atac.peaks.hg19, delim = c("-", "-"))
sum(lengths(atac.peaks.map) > 1)

## Only keep uniquely mapping peaks
atac.peaks.map <- atac.peaks.map[lengths(atac.peaks.map) == 1]
atac.peaks.map <- unlist(atac.peaks.map)
hist(width(atac.peaks.map))


## Switch over to hg38 annotations
atac.cell.counts <- atac.cell.counts[names(atac.peaks.map),]
rownames(atac.cell.counts) <- granges2peak(atac.peaks.map)
dim(atac.cell.counts)

atac.peaks <- atac.peaks.map
names(atac.peaks) <- atac.peaks$Peak <- rownames(atac.cell.counts)

## Save hg38 results
saveRDS(atac.cell.counts, file = "../Data/scTHS_counts.Robj")

## Build and export Seurat object
chrom_assay <- CreateChromatinAssay(
  counts = atac.cell.counts,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 0,
  min.features = 0
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

## Add UMAP
atac[["umap"]] <- CreateDimReducObject(embeddings = atac.umap.emb, key = "UMAP_", assay = "peaks")

## Add clusters
atac$mapped.ident <- atac.mapped.ident

## Export
saveRDS(atac, file = "../Data/sciATAC_seurat_object.Robj")
