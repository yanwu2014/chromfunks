library(swne)
library(perturbLM)
library(Seurat)
library(cicero)
library(cisTopic)
library(chromfunks)
library(Signac)


## Output RData environment
output.file <- "scTHS_annotations.RData"



#### Map THS-Seq Clusters to RNA-Seq Clusters ####

## Load reference data
fbcx <- readRDS("../Data/snDrop_RNA_seurat_object.Robj")
ref.annot <- fbcx$annot
table(ref.annot)

## Plot scRNA-seq clusters
rna.umap.emb <- Embeddings(fbcx, "umap")
rna.time.point <- fbcx$age

## Get colors
library(scales)
plot.seed <- 65463947
set.seed(plot.seed)
rna.cluster.colors <- sample(hue_pal()(nlevels(ref.annot)))
names(rna.cluster.colors) <- levels(ref.annot)

pdf("scRNA_fCTX-combined_umap_clusters.pdf", height = 6, width = 6)
PlotDims(rna.umap.emb[names(ref.annot),], factor(ref.annot), alpha.plot = 0.5,
         pt.size = 0.75, label.size = 4.5, show.legend = F, show.axes = F, 
         colors.use = rna.cluster.colors)
dev.off()

pdf("scRNA_fCTX-combined_umap_clusters_nolabel.pdf", height = 6, width = 6)
PlotDims(rna.umap.emb[names(ref.annot),], factor(ref.annot), alpha.plot = 0.5,
         pt.size = 0.75, label.size = 0, show.legend = F, show.axes = F, 
         colors.use = rna.cluster.colors)
dev.off()

levels(rna.time.point) <- c("wk16", "wk18")
pdf("scRNA_fCTX-combined_umap_timepoint.pdf", width = 7, height = 6)
PlotDims(rna.umap.emb[names(rna.time.point),], factor(rna.time.point), alpha.plot = 0.4,
         pt.size = 0.5, label.size = 4.5, show.legend = T, show.axes = F, seed = plot.seed,
         do.label = F, use.brewer.pal = F) + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 16))
dev.off()


## Load chromatin data
load("../Data/scTHS_cisTopic_results.RData")

# ## Fine tune the clusters
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
# clusters <- Idents(seurat_obj)
# table(clusters)

pdf("fCTX-combined_umap_clusters_cisTopic.pdf", width = 7.5, height = 6)
PlotDims(umap.emb, sample.groups = clusters, x.lab = "umap1", y.lab = "umap2",
         pt.size = 0.3, alpha.plot = 0.4, label.size = 4.5, show.legend = T, show.axes = F,
         seed = 35246435)
dev.off()

## Normalize chromatin/rna datasets
norm_ga <- ScaleCounts(unnorm_ga)
norm_rna <- ScaleCounts(GetAssayData(fbcx, "counts"))
clusters <- clusters[colnames(norm_ga)]

## Load scRNA-seq marker genes
fbcx <- SetIdent(fbcx, value = ref.annot)
# scRNA.markers.df <- FindAllMarkers(fbcx, logfc.threshold = 0.25, return.thresh = 1e-3, only.pos = T)
# write.table(scRNA.markers.df, file = "scRNA-fCTX-combined_marker_genes.tsv", sep = "\t")
scRNA.markers.df <- read.table("../Data/snDrop_marker_genes.tsv.gz", header = T)
table(scRNA.markers.df$cluster)

## Compute gene activities
ga.obj <- CreateSeuratObject(unnorm_ga)
ga.obj <- NormalizeData(ga.obj)
ga.obj <- FindVariableFeatures(ga.obj)

## Correlate gene expression and gene activity
var.genes <- unique(c(VariableFeatures(fbcx), VariableFeatures(ga.obj), as.character(scRNA.markers.df$gene)))
var.genes <- var.genes[var.genes %in% rownames(norm_ga) & var.genes %in% rownames(fbcx)]
length(var.genes)

scale_ga <- as.matrix(norm_ga[var.genes,] - Matrix::rowMeans(norm_ga[var.genes,]))
scale_rna <- as.matrix(norm_rna[var.genes,] - Matrix::rowMeans(norm_rna[var.genes,]))

scRNA_centroids <- t(apply(scale_rna, 1, function(x) tapply(x, ref.annot, mean)))
ga_centroids <- t(apply(scale_ga, 1, function(x) tapply(x, clusters, mean)))

rna.ga.cor <- sapply(colnames(ga_centroids), function(x) {
  sapply(colnames(scRNA_centroids), function(y) cor(ga_centroids[,x], scRNA_centroids[,y]))
})
rna.ga.cor[abs(rna.ga.cor) < 0.2] <- 0

pdf("fCTX-combined_ga_cluster_cor_cisTopic.pdf", width = 7, height = 5)
ggHeat(rna.ga.cor, clustering = "both", x.lab.size = 12, y.lab.size = 12)
dev.off()

## Save results
save.image(output.file)


#### Rename clusters and create umap plot with named clusters ####
cluster.mapping.df <- read.table("fCTX-combined_cluster_mapping.txt", header = T, stringsAsFactors = F)
cluster.mapping <- cluster.mapping.df$named_cluster; names(cluster.mapping) <- cluster.mapping.df$seurat_cluster;
named.clusters <- plyr::revalue(clusters, replace = cluster.mapping)
table(named.clusters)

## Create metadata
# metadata.df <- data.frame(cell = names(clusters), orig.ident = clusters,
#                           named.ident = named.clusters,
#                           time.point = time.point)
# write.table(metadata.df, file = "../Data/scTHS_metadata.tsv", sep = "\t")

metadata.df <- read.table("../Data/scTHS_metadata.tsv", header = T, sep = "\t")
named.clusters <- metadata.df$named.ident; names(named.clusters) <- metadata.df$cell;

# sub.clusters <- as.character(named.clusters); names(sub.clusters) <- names(named.clusters);
# rg.metadata.df <- read.table("../Data/scTHS_RG_subclustering_metadata.tsv", header = T, sep = "\t")
# rg.sublusters <- as.character(rg.metadata.df$named.ident); names(rg.sublusters) <- rg.metadata.df$cell;
# sub.clusters[names(rg.sublusters)] <- rg.sublusters
# 
# In.metadata.df <- read.table("../Data/scTHS_In_subclustering_metadata.tsv", header = T, sep = "\t")
# In.sublusters <- as.character(In.metadata.df$named.ident); names(In.sublusters) <- In.metadata.df$cell;
# sub.clusters[names(In.sublusters)] <- In.sublusters
# 
# sub.clusters <- factor(sub.clusters)
# table(sub.clusters)
# 
# metadata.df$sub.named.ident <- sub.clusters
# write.table(metadata.df, file = "../Data/scTHS_RG_subclustering_metadata.tsv", sep = "\t")

plot.seed <- 35235

## Plot named clusters
pdf("fCTX-combined_umap_named_subclusters_cisTopic.pdf", width = 7.5, height = 6)
PlotDims(umap.emb, sample.groups = sub.clusters, pt.size = 0.2, alpha.plot = 0.35, 
         label.size = 5, show.legend = T, show.axes = F, 
         colors.use = rna.cluster.colors)
dev.off()

pdf("fCTX-combined_umap_named_subclusters_cisTopic_nolabel.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = sub.clusters, pt.size = 0.2, alpha.plot = 0.35, 
         label.size = 0, show.legend = F, show.axes = F, 
         colors.use = rna.cluster.colors)
dev.off()

library(ggplot2)
levels(time.point) <- c("wk16", "wk18")
pdf("fCTX-combined_umap_timepoint_cisTopic.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = time.point, pt.size = 0.025, 
         alpha.plot = 0.5, do.label = F, show.legend = F, show.axes = F, 
         seed = 124142, use.brewer.pal = F) + 
  theme(legend.text = element_text(size = 16))
dev.off()



## Save results
save.image(output.file)

## Write cluster colors
rna.cluster.colors.df <- data.frame(cluster = names(rna.cluster.colors), color = rna.cluster.colors)
write.table(rna.cluster.colors.df, file = "../Data/harmonized_cluster_colors.txt", sep = "\t", col.names = T,
            row.names = F, quote = F)



#### Create marker dotplots for RNA/THS-Seq ####


## Load scRNA markers
rna.markers.df <- read.table("../Data/snDrop_marker_genes.tsv.gz", header = T, row.names = 1)
head(subset(rna.markers.df, cluster == "InCGE" & pct.1/pct.2 > 3), n = 10)

## Feature & Violin plots
gene <- "EOMES"
FeaturePlotDims(umap.emb, norm_ga[gene,], feature.name = paste0(gene, "\nAct"), show.axes = F,
                pt.size = 0.2, alpha = 0.4, quantiles = c(0.001, 0.999))

FeaturePlotDims(rna.umap.emb, fbcx@data[gene,], feature.name = paste0(gene, "\nExpr"), show.axes = F,
                pt.size = 0.5, alpha = 0.4, quantiles = c(0.001, 0.999))


VlnPlot(ga.obj, gene, point.size.use = 0.01, x.lab.rot = T)
VlnPlot(fbcx, gene, point.size.use = 0.01, x.lab.rot = T)

# ga.markers.df <- FindAllMarkers(ga.obj, logfc.threshold = 0.25, only.pos = T)
cluster <- "InCGE"
cl.ga.markers.df <- FindMarkers(ga.obj, ident.1 = "InCGE", logfc.threshold = 0.3, only.pos = T)
# head(subset(ga.markers.df, cluster == "CycProg"), n = 20)

## Final list of genes to plot
genes.plot <- c("GLI3", ## RG markers
                "HOPX", ## oRG marker
                "CRYAB", ## vRG marker
                "HMGB2", ## CycProg marker
                "EOMES", ## ipEx/CycProg marker
                "PPP1R17", ## ipEx markers
                "CNTNAP2", ## ebEx marker 
                "CUX2", ## ebEx/ExL2/3 marker
                "SATB2", ## ExL2/3/4 marker
                "FOXP1", "ADAMTSL3", ## ExL4 markers
                "SEMA3E", "KLHL32", ## ExL5/6 markers
                "CDH18", "KIAA1217", ## ExL6b markers
                "DLX1", "SST", ## In markers
                "PDGFRA", "OLIG1", ## OPC markers
                "SAMSN1", "P2RY12" ## Mic markers
                )
all(genes.plot %in% rownames(ga.obj))

sub.cluster.order <- c("oRG", "vRG", "CycProg", "ipEx", "ebEx", "ExL2/3",
                       "ExL4", "ExL5/6", "ExL6b", "InMGE", "InCGE", "OPC", "Mic")
ga.ord.ident <- factor(sub.clusters, levels = sub.cluster.order)
table(ga.ord.ident)

scale_ga_markers <- as.matrix(norm_ga[genes.plot, names(ga.ord.ident)] - Matrix::rowMeans(norm_ga[genes.plot, names(ga.ord.ident)]))
scale_ga_markers <- apply(scale_ga_markers, 1, function(x) tapply(x, ga.ord.ident, mean))

pdf("fCTX-combined_ga_heatmap.pdf", width = 7, height = 4.5)
ggHeat(scale_ga_markers, x.lab.size = 11, y.lab.size = 11)
dev.off()

ref.annot <- factor(fbcx$annot, levels = c(sub.cluster.order, "End"))
ref.annot <- ref.annot[!is.na(ref.annot)]
fbcx <- fbcx[,names(ref.annot)]
fbcx <- SetIdent(fbcx, value = ref.annot)
rna.genes.plot <- c(genes.plot[genes.plot %in% rownames(fbcx)], 
                    "IGFBP7", "FN1") ## End markers

scale_rna_markers <- as.matrix(norm_rna[rna.genes.plot, names(ref.annot)] - Matrix::rowMeans(norm_rna[rna.genes.plot, names(ref.annot)]))
scale_rna_markers <- t(apply(scale_rna_markers, 1, tapply, ref.annot, mean))

pdf("fCTX-combined_rna_heatmap.pdf", width = 7, height = 4.5)
ggHeat(t(rna.markers.heat), x.lab.size = 11, y.lab.size = 11)
dev.off()


## Correlation heatmap between RNA and named chromatin clusters
named_ga_centroids <- t(apply(scale_ga, 1, function(x) tapply(x, ga.ord.ident, mean)))

scRNA_centroids <- scRNA_centroids[,c(sub.cluster.order, "End")]
rna.named.ga.cor <- sapply(colnames(named_ga_centroids), function(x) {
  sapply(colnames(scRNA_centroids), function(y) cor(named_ga_centroids[,x], scRNA_centroids[,y]))
})

pdf("fCTX-combined_named_ga_cluster_cor_cisTopic.pdf", width = 7, height = 5)
ggHeat(rna.named.ga.cor, clustering = "none", x.lab.size = 12, y.lab.size = 12)
dev.off()

## Save results
save.image(output.file)


#### Format scTHS and snDrop datasets and metadata as Seurat objects ####

## Create scTHS seurat object
chrom_assay <- CreateChromatinAssay(
  counts = cisTopicObject@binary.count.matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 0,
  min.features = 0
)

ths <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata.df[cisTopicObject@cell.names,]
)

## Add UMAP
ths[["umap"]] <- CreateDimReducObject(embeddings = umap.emb, key = "UMAP_", assay = "peaks")

## Add cisTopic topics
cisTopicObject <- getRegionsScores(cisTopicObject)
cisTopicObject <- binarizecisTopics(cisTopicObject)
topic.emb <- modelMatSelection(cisTopicObject, "cell", "Probability")
# topic.load <- modelMatSelection(cisTopicObject, "region", "Probability")
ths[["topics"]] <- CreateDimReducObject(embeddings = umap.emb, key = "topic_", assay = "peaks")

ga_assay <- CreateAssayObject(unnorm_ga[,cisTopicObject@cell.names])
ths[["RNA"]] <- ga_assay

## Write scTHS data to file
saveRDS(ths, file = "../Data/scTHS_chromatin_seurat_object.Robj")
