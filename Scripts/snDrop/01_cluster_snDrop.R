library(Seurat)
library(pagoda2)
library(swne)

## Output file
output.file <- "../Data/snDrop_analysis.RData"

counts <- ReadData("../Data/snDrop_counts.tsv.gz")
meta.data <- read.table("../Data/snDrop_metadata.tsv.gz", header = T, row.names = 1)

annot <- as.character(meta.data$annot); names(annot) <- rownames(meta.data);
annot.mapping <- c("Ex-CUX2"="ExL2/3", "Ex-BHLHE22"="ExL4", "Ex-FOXP1" = "ExL4",
                   "Ex-ETV1"="ExL5", "Ex-ASIC2"="ExL5", "Ex-FOXP2"="ExL5/6",
                   "Ex-CDH18"="ExL6b", "ebEX-NRP1"="ebEx-NRP1")
annot <- plyr::revalue(annot, replace = annot.mapping)
counts <- counts[,names(annot)]
table(factor(annot))

time.point <- meta.data$age; names(time.point) <- rownames(meta.data);
table(time.point)

batch <- sapply(colnames(counts), ExtractField, field = 1, delim = "_")
table(batch)

counts <- counts[Matrix::rowSums(counts > 0) > 10,]
dim(counts)

rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts, n.cores = 4, trim = 10, batch = batch)
r$adjustVariance(plot = T, gam.k = 10)
r$calculatePcaReduction(nPcs = 150, n.odgenes = 2e3, maxit = 1000)

fbcx <- CreateSeuratObject(counts, meta.data = meta.data[colnames(counts),])
fbcx <- NormalizeData(fbcx)
VariableFeatures(fbcx) <- rownames(r$misc$PCA$v)

pc.scores <- r$reductions$PCA
pca_dr <- CreateDimReducObject(embeddings = pc.scores,
                               key = "PC_", assay = DefaultAssay(fbcx),
                               misc = list(genes.use = r$misc$odgenes,
                                           cells.use = colnames(fbcx)))
fbcx[["pca"]] <- pca_dr

# fbcx <- FindNeighbors(fbcx, dims = 1:50, k.param = 10)
# fbcx <- FindClusters(fbcx, resolution = 1)
# clusters <- Idents(fbcx)
# table(clusters)

fbcx <- RunUMAP(fbcx, dims = 1:50, min.dist = 0.2, n.neighbors = 30)

umap.emb <- Embeddings(fbcx, "umap")
PlotDims(umap.emb[names(annot),], sample.groups = annot, alpha = 0.4, pt.size = 0.5, 
         show.axes = F, seed = 12423)


## Map subclusters to RNA data
raw_counts_mat <- ReadData("../Data/ref_fetal_cortex_counts.tsv.gz")
raw_counts_mat <- FilterData(raw_counts_mat, min.samples.frac = 0.001, trim = 5, 
                             min.nonzero.features = 0, max.sample.sum = Inf)
dim(raw_counts_mat)

ref.meta.data <- read.table("../Data/ref_fetal_cortex_metadata.csv.gz", header = T, row.names = 1, sep = ",")
ref.annot <- ref.meta.data$Cluster; names(ref.annot) <- rownames(ref.meta.data)
table(ref.annot)

donor <- factor(ref.meta.data$Donor); names(donor) <- rownames(ref.meta.data);
ref.var.genes.list <- lapply(levels(donor), function(i) {
  donor.counts <- raw_counts_mat[,names(donor[donor == i])]
  donor.counts <- FilterData(donor.counts, min.samples.frac = 5e-4, trim = 3, 
                             min.nonzero.features = 0, max.sample.sum = Inf)
  SelectFeatures(donor.counts, n.features = 2000)
})
ref.var.genes <- unique(unlist(ref.var.genes.list, F, F))
length(ref.var.genes)

## Correlate gene expression
var.genes <- union(ref.var.genes, VariableFeatures(fbcx))
var.genes <- var.genes[var.genes %in% rownames(fbcx) & var.genes %in% rownames(raw_counts_mat)]

scale_rna <- ScaleCounts(counts, batch = batch)
scale_rna <- as.matrix(scale_rna[var.genes,] - Matrix::rowMeans(scale_rna[var.genes,]))

scale_ref <- ScaleCounts(raw_counts_mat)
scale_ref <- as.matrix(scale_ref[var.genes,] - Matrix::rowMeans(scale_ref[var.genes,]))

rna_centroids <- t(apply(scale_rna[,names(annot)], 1, function(x) tapply(x, annot, mean)))
ref_centroids <- t(apply(scale_ref[,names(ref.annot)], 1, function(x) tapply(x, ref.annot, mean)))

rna.ref.cor <- sapply(colnames(rna_centroids), function(x) {
  sapply(colnames(ref_centroids), function(y) cor(rna_centroids[,x], ref_centroids[,y]))
})

pdf("scRNA_dev_cortex_cor_heatmap.pdf", width = 6, height = 6)
ggHeat(rna.ref.cor, clustering = "both", x.lab.size = 11, y.lab.size = 11)
dev.off()

## Rename some of our clusters based on this mapping
dev.annot.mapping <- c("ipEx"="CycProg", "ebEx-EOMES"="ipEx", "ebEx-NRP1"="ebEx",
                       "In-GALNTL6"="InMGE", "In-ADARB2"="InCGE", "In-PTPRT"="InMGE")
dev.annot <- plyr::revalue(annot, replace = dev.annot.mapping)
table(dev.annot)

PlotDims(umap.emb, sample.groups = dev.annot, alpha = 0.4, pt.size = 0.5, show.axes = F,
         seed = 2132085)

## Now let's try and separate oRG and vRG
fbcx.rg <- fbcx[,names(dev.annot[dev.annot == "RG"])]
fbcx.rg <- FindNeighbors(fbcx.rg, k.param = 10, dims = 1:50)
fbcx.rg <- FindClusters(fbcx.rg, resolution = 0.5)
fbcx.rg <- FindVariableFeatures(fbcx.rg)
table(Idents(fbcx.rg))

rg.ref.ident <- droplevels(ref.annot[ref.annot %in% c("oRG", "vRG")])
rg.var.genes <- union(SelectFeatures(raw_counts_mat[,names(rg.ref.ident)]), VariableFeatures(fbcx.rg))
rg.var.genes <- rg.var.genes[rg.var.genes %in% rownames(fbcx.rg) & rg.var.genes %in% rownames(raw_counts_mat)]

rg.scale_rna <- ScaleCounts(counts[,colnames(fbcx.rg)], batch = batch)
rg.scale_rna <- as.matrix(rg.scale_rna[rg.var.genes,] - Matrix::rowMeans(rg.scale_rna[rg.var.genes,]))

rg.scale_ref <- ScaleCounts(raw_counts_mat[,names(rg.ref.ident)])
rg.scale_ref <- as.matrix(rg.scale_ref[rg.var.genes,] - Matrix::rowMeans(rg.scale_ref[rg.var.genes,]))

rg.rna_centroids <- t(apply(rg.scale_rna, 1, function(x) tapply(x, Idents(fbcx.rg), mean)))
rg.ref_centroids <- t(apply(rg.scale_ref[,names(rg.ref.ident)], 1, function(x) tapply(x, rg.ref.ident, mean)))

rg.rna.ref.cor <- sapply(colnames(rg.rna_centroids), function(x) {
  sapply(colnames(rg.ref_centroids), function(y) cor(rg.rna_centroids[,x], rg.ref_centroids[,y]))
})

pdf("scRNA_dev_cortex_cor_heatmap_RG_subclusters.pdf", width = 3, height = 1.5)
ggHeat(rg.rna.ref.cor, x.lab.size = 11, y.lab.size = 11)
dev.off()

rg.ident <- plyr::revalue(Idents(fbcx.rg), replace = c("0"="oRG", "1"="vRG", "2"="vRG", "3"="oRG"))
final.annot <- as.character(dev.annot); names(final.annot) <- names(dev.annot);
final.annot[names(rg.ident)] <- as.character(rg.ident)
final.annot <- factor(final.annot)
table(final.annot)

fbcx$annot <- final.annot
fbcx <- SetIdent(fbcx, value = fbcx$annot)
markers.df <- FindAllMarkers(fbcx, logfc.threshold = 0.25)

## Save results
save.image(output.file)
save(fbcx, file = "../Data/snDrop_RNA_seurat_object.Robj")
write.table(markers.df, file = "../Data/snDrop_marker_genes.tsv", sep = "\t")


## Correlate finalized RNA cell types with dev cortex datasets
library(cellMapper)

rna.ref.cor <- MapClustersCor(query.data = scale_rna[,colnames(fbcx)], query.cluster = fbcx$annot, 
                              train.data = scale_ref[,names(ref.annot)], train.clusters = ref.annot)
rna.ref.cor.mat <- rna.ref.cor$cluster.correlations

fetal.cl.order <- c("oRG", "vRG", "CycProg", "ipEx", "ebEx",
                    "ExL2/3", "ExL4", "ExL5", "ExL5/6", "ExL6b",
                    "InCGE", "InMGE", 
                    "OPC", "Mic", "End")
dev.cortex.cl.order <- c("oRG", "vRG", "PgG2M", "PgS", "IP", 
                         "ExN", "ExM", "ExM-U", "ExDp1", "ExDp2", 
                         "InCGE", "InMGE", "OPC", "Mic", "End", "Per")
rna.ref.cor.mat <- rna.ref.cor.mat[dev.cortex.cl.order, fetal.cl.order]
# rna.ref.cor.mat <- rna.ref.cor.mat[rev(rownames(rna.ref.cor.mat)),]

pdf("scRNA_final_annot_dev_cortex_cor_heatmap.pdf", width = 6, height = 5.25)
ggHeat(rna.ref.cor.mat, clustering = "none", x.lab.size = 11, y.lab.size = 11)
dev.off()
