library(swne)
library(Seurat)
library(Signac)
library(cisTopic)
library(chromfunks)
library(cellMapper)
library(ggplot2)


## Cell type ordering for visualizations
dev.cl.order <- c("oRG", "vRG", "PgG2M", "PgS", "IP", 
                  "ExN", "ExM", "ExM-U", "ExDp1", "ExDp2", 
                  "InCGE", "InMGE", "OPC", "Mic", "End", "Per")
adult.cl.order <- c("Ast", "Ex", "In", "Oli", "Opc", "Mic", "End")
atac.rna.cl.order <- c("RG", "ipEx", "ebEx", "Ex", "OPC", "In", "Mic", "End")


#### Load RNA datasets ####

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


## Load and downsample ATAC RNA data
atac.rna.counts <- readRDS("../Data/sciRNA_gene_counts.Robj")
atac.rna.metadata <- readRDS("../Data/sciRNA_orig_cell_metadata.Robj")
atac.rna.metadata <- subset(atac.rna.metadata, Organ == "Cerebrum")
table(atac.rna.metadata$Main_cluster_name)

## Rename features
gene.df <- readRDS("../Data/sciRNA_gene_metadata.Robj")
gene.df <- subset(gene.df, !is.na(gene_short_name))
atac.rna.counts <- atac.rna.counts[rownames(atac.rna.counts) %in% rownames(gene.df),]
rownames(atac.rna.counts) <- gene.df[rownames(atac.rna.counts), "gene_short_name"]

## Downsample cells to at most 10k per cluster
atac.rna.orig.ident <- atac.rna.metadata$Main_cluster_name
names(atac.rna.orig.ident) <- atac.rna.metadata$sample
atac.rna.orig.ident <- atac.rna.orig.ident[!is.na(atac.rna.orig.ident)]
atac.rna.orig.ident <- UnflattenGroups(atac.rna.orig.ident)

set.seed(457412535)
atac.rna.orig.ident <- lapply(atac.rna.orig.ident, function(x) sample(x, min(10000, length(x))))
atac.rna.orig.ident <- FlattenGroups(atac.rna.orig.ident)
table(atac.rna.orig.ident)

rownames(atac.rna.metadata) <- atac.rna.metadata$sample
atac.rna.metadata <- atac.rna.metadata[names(atac.rna.orig.ident),]
atac.rna.counts <- atac.rna.counts[,names(atac.rna.orig.ident)]
dim(atac.rna.counts)

## Create Seurat object
atac.rna <- CreateSeuratObject(atac.rna.counts, meta.data = atac.rna.metadata,
                               min.cells = 10)
atac.rna <- NormalizeData(atac.rna)
Idents(atac.rna) <- atac.rna.orig.ident

atac.umap.emb <- cbind(atac.rna$Main_cluster_umap_1, atac.rna$Main_cluster_umap_2)
colnames(atac.umap.emb) <- c("umap1", "umap2")
atac.rna[["umap"]] <- CreateDimReducObject(atac.umap.emb, key = "UMAP_", assay = "RNA")

DimPlot(atac.rna, reduction = "umap")

## Clean up unused objects
rm(atac.rna.counts, atac.rna.metadata); invisible(gc());



#### Map ATAC RNA clusters to THS RNA Clusters ####

## Subcluster ATAC Excitatory neurons
atac.rna.ex <- atac.rna[,Idents(atac.rna) == "Excitatory neurons"]
atac.rna.ex <- NormalizeData(atac.rna.ex)
atac.rna.ex <- FindVariableFeatures(atac.rna.ex, nfeatures = 2000)
atac.rna.ex <- ScaleData(atac.rna.ex)
atac.rna.ex <- RunPCA(atac.rna.ex, npcs = 20, verbose = F)
atac.rna.ex <- FindNeighbors(atac.rna.ex)
atac.rna.ex <- FindClusters(atac.rna.ex, resolution = 0.8)

## Plot results
DimPlot(atac.rna.ex, reduction = "umap")
table(Idents(atac.rna.ex))

## Add finer resolution excitatory neuron clusters to annotations
atac.rna.orig.ident[names(Idents(atac.rna.ex))] <- paste0("Ex", Idents(atac.rna.ex))
table(atac.rna.orig.ident)

## Merge variable genesets
genes.map <- union(VariableFeatures(dev), VariableFeatures(vctx))
length(genes.map)

## Scale datasets
atac.rna <- ScaleData(atac.rna, features = genes.map)
dev <- ScaleData(dev, features = genes.map)
vctx <- ScaleData(vctx, features = genes.map)

## Map clusters
atac.dev.cor <- MapClustersCor(query.data = GetAssayData(atac.rna, assay = "RNA", slot = "scale.data"),
                               query.clusters = atac.rna.orig.ident,
                               train.data = GetAssayData(dev, slot = "scale.data"),
                               train.clusters = Idents(dev))
atac.dev.cor.mat <- atac.dev.cor$cluster.correlations[dev.cl.order,]

atac.vctx.cor <- MapClustersCor(query.data = GetAssayData(atac.rna, assay = "RNA", slot = "scale.data"),
                                query.clusters = atac.rna.orig.ident,
                                train.data = GetAssayData(vctx, slot = "scale.data"),
                                train.clusters = adult.rna.idents)
atac.vctx.cor.mat <- atac.vctx.cor$cluster.correlations[adult.cl.order,]

## Visualize cluster mappings
atac.rna.cor.mat <- rbind(atac.dev.cor.mat, atac.vctx.cor.mat)

pdf("fCTX-combined_ATAC_rna_orig_ident_dev_cortex_cor.pdf", width = 7, height = 5)
ggHeat(t(atac.rna.cor.mat))
dev.off()

## Visualize key marker gene expression
marker.genes <- c("GLI3", "VIM", "NES", "FABP7", "SOX2", ## RG markers
                  "EOMES", "PPP1R17", "NEUROG2", "PAX6", ## ipEx markers,
                  "S100B", "NEU1" ## Ast markers
)

## Make heatmap
atac.rna <- ScaleData(atac.rna, features = marker.genes, assay = "RNA")
atac.marker.heat <- t(apply(GetAssayData(atac.rna, assay = "RNA", slot = "scale.data"), 1, function(x) {
  tapply(x, atac.rna.orig.ident, mean)
}))
atac.marker.heat <- atac.marker.heat[marker.genes,]

max.scale <- 2
atac.marker.heat[atac.marker.heat > max.scale] <- max.scale

pdf("fCTX-combined_ATAC_rna_orig_ident_markers.pdf", width = 5, height = 5)
ggHeat(t(atac.marker.heat))
dev.off()



#### Remap clusters and make final plots ####

## Rename idents
atac.rna.map <- c("Astrocytes"="RG", 
                  "Oligodendrocytes"="OPC",
                  "Ex0"="Ex",
                  "Ex1"="ebEx",
                  "Ex2"="Ex",
                  "Ex3"="RG",
                  "Ex4"="Ex",
                  "Ex5"="Ex",
                  "Ex6"="Ex",
                  "Ex7"="ipEx",
                  "Ex8"="Ex",
                  "Ex9"="Ex",
                  "Ex10"="Ex",
                  "Ex11"="Ex",
                  "Ex12"="Ex",
                  "Ex13"="ipEx",
                  "Inhibitory neurons"="In", 
                  "Limbic system neurons"="LimbicNeurons",
                  "SKOR2_NPSR1 positive cells"="SKOR2_NPSR1_Cells",
                  "Vascular endothelial cells"="End",
                  "Microglia"="Mic")
atac.rna.ident <- plyr::revalue(atac.rna.orig.ident, replace = atac.rna.map)
table(atac.rna.ident)

## Remake UMAP plot
Idents(atac.rna) <- atac.rna.ident

mapped.colors.df <- read.table("../Data/harmonized_cluster_colors.txt", header = T, sep = "\t",
                               stringsAsFactors = F, comment.char = "")
mapped.colors <- mapped.colors.df$color
names(mapped.colors) <- mapped.colors.df$cluster
mapped.colors[["Ex"]] <- mapped.colors[["ExL2/3"]]
mapped.colors[["In"]] <- mapped.colors[["InCGE"]]
mapped.colors[["RG"]] <- mapped.colors[["oRG"]]
mapped.colors[["RG/OPC"]] <- mapped.colors[["OPC"]]

pdf("fCTX-combined_ATAC_rna_renamed_ident_umap.pdf", width = 4, height = 4)
cells.plot <- names(atac.rna.ident[atac.rna.ident %in% atac.rna.cl.order])
cells.sample <- sample(cells.plot, size = 10000)
PlotDims(atac.umap.emb[cells.sample,], sample.groups = atac.rna.ident[cells.sample], 
         pt.size = 0.5, alpha = 0.5, do.label = T, show.legend = F,
         colors.use = mapped.colors, show.axes = F, label.size = 5)
dev.off()

## Remap ATAC RNA data to snDropSeq datasets
atac.rna <- ScaleData(atac.rna, features = genes.map)
atac.dev.mapped.cor <- MapClustersCor(query.data = GetAssayData(atac.rna, assay = "RNA", slot = "scale.data"),
                                      query.clusters = atac.rna.ident,
                                      train.data = GetAssayData(dev, slot = "scale.data"),
                                      train.clusters = Idents(dev))
atac.dev.mapped.cor.mat <- atac.dev.mapped.cor$cluster.correlations[dev.cl.order,]

atac.vctx.mapped.cor <- MapClustersCor(query.data = GetAssayData(atac.rna, assay = "RNA", slot = "scale.data"),
                                       query.clusters = atac.rna.ident,
                                       train.data = GetAssayData(vctx, slot = "scale.data"),
                                       train.clusters = adult.rna.idents)
atac.vctx.mapped.cor.mat <- atac.vctx.mapped.cor$cluster.correlations[adult.cl.order,]
atac.rna.mapped.cor.mat <- rbind(atac.dev.mapped.cor.mat, atac.vctx.mapped.cor.mat)

## Visualize cluster map
pdf("fCTX-combined_ATAC_rna_DropSeq_cor.pdf", width = 6.75, height = 3)
ggHeat(t(atac.rna.mapped.cor.mat[,rev(atac.rna.cl.order)]))
dev.off()

## Remake heatmap with remapped idents
atac.rna <- ScaleData(atac.rna, features = marker.genes, assay = "RNA")
atac.mapped.marker.heat <- t(apply(GetAssayData(atac.rna, assay = "RNA", slot = "scale.data"), 1, function(x) {
  tapply(x, atac.rna.ident, mean)
}))
atac.mapped.marker.heat <- atac.mapped.marker.heat[marker.genes,]

max.scale <- 2
atac.mapped.marker.heat[atac.mapped.marker.heat > max.scale] <- max.scale

pdf("fCTX-combined_ATAC_rna_marker_heatmap.pdf", width = 4.75, height = 3)
ggHeat(t(atac.mapped.marker.heat[,rev(atac.rna.cl.order)]))
dev.off()


## Export results
atac.rna$ident <- atac.rna.ident
atac.rna <- atac.rna[,Idents(atac.rna) %in% atac.rna.cl.order]
dim(atac.rna)

saveRDS(atac.rna, file = "../Data/sciATAC_seurat_object.Robj")
