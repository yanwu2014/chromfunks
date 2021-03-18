library(methods)
library(cisTopic)
library(Seurat)
library(swne)
library(cicero)
library(chromfunks)

## Output file name
output.file <- "../Data/scTHS_In_subclustering.RData"

## load RData
counts <- readRDS("../Data/scTHS_counts.Robj")

## Load fraction of unique reads in peaks
nReads <- Matrix::colSums(counts)

## Binarize matrix
counts@x[counts@x > 1] <- 1

## Load existing metadata and subset to only interneurons
meta.data <- read.table("../scTHS_metadata.tsv.gz", sep = "\t")
named.ident <- meta.data$named.ident; names(named.ident) <- meta.data$cell;

cells.keep <- as.character(subset(meta.data, named.ident == "In")$cell)
counts <- counts[,cells.keep]

## Filter for common peaks
min.cells <- 10
peaks.idx <- Matrix::rowSums(counts) > min.cells
paste0("Peaks to keep: ", sum(peaks.idx))
counts <- counts[peaks.idx,]

## Get cell metadata
file.names <- colnames(counts)
time.point <- factor(sapply(file.names, ExtractField, field = 3, delim = "#"))
names(time.point) <- file.names; levels(time.point) <- c("wk16", "wk18");

## initialize cisTopic object from count matrix
cisTopicObject <- createcisTopicObject(counts, project.name = 'fCTX', keepCountsMatrix = F)

## run LDA model
topics.range <- c(20,25,30)
cisTopicObject <- runModels(cisTopicObject, topic = topics.range, seed = 987, nCores = length(topics.range),
                            burnin = 120, iterations = 250, addModels = F)

## select for model
cisTopicObject <- cisTopic::selectModel(cisTopicObject)

## run UMAP
cfg <- umap::umap.defaults; cfg$metric <- "cosine"; cfg$min_dist <- 0.3;
cisTopicObject <- runUmap(cisTopicObject, target = 'cell', method = 'Z-score', config = cfg)

## pull out umap coordinates
umap.emb <- cisTopicObject@dr$cell[["Umap"]]

## pull out topic coordinates
topic.emb <- t(modelMatSelection(cisTopicObject, target = "cell", method = 'Z-score'))
n.dims <- ncol(topic.emb)

## create Seurat object for clustering
seurat_obj = CreateSeuratObject(counts, meta.data = data.frame(time.point))
seurat_obj[["topic"]] <- CreateDimReducObject(embeddings = topic.emb, key = "topic_", 
                                              assay = DefaultAssay(seurat_obj))

seurat_obj <- FindNeighbors(seurat_obj, k.param = 10, dims = 1:n.dims, reduction = "topic")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, verbose = F)
clusters <- Idents(seurat_obj)
table(clusters)

## Make gene activity matrix
pData <- new("AnnotatedDataFrame", data = data.frame(time.point, clusters))

chrom <- sapply(rownames(counts), ExtractField, field = 1, delim = ":")
loc <- sapply(rownames(counts), ExtractField, field = 2, delim = ":")
loc_start <- as.integer(sapply(loc, ExtractField, field = 1, delim = "-"))
loc_end <- as.integer(sapply(loc, ExtractField, field = 2, delim = "-"))

fData <- data.frame(site_name = rownames(counts), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
rownames(fData) <- rownames(counts)
fData <- new("AnnotatedDataFrame", data = fData)

input_cds <- newCellDataSet(counts, phenoData = pData, featureData = fData,
                            expressionFamily = VGAM::binomialff())

# cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.emb, k = 40)
# conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run
# print("Done finding connections")
load("../Data/scTHS_conns.RData")
conns <- subset(conns, Peak1 %in% rownames(counts) & Peak2 %in% rownames(counts))

## Format peak annotation dataframe
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

peaks.gr <- makeGRangesFromDataFrame(data.frame(seqnames = chrom, start = loc_start, end = loc_end))
peak_anno <- annotatePeak(peaks.gr, tssRegion = c(-5000, 5000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb = "org.Hs.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df[!grepl("Promoter|Exon|Intron", peak_anno_df$annotation), "SYMBOL"] <- NA
peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")
peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]
head(peak_anno_df)

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

## Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)
print(dim(unnorm_ga))

## Save results
save.image(output.file)

## Map to InCGE/InMGE
raw_counts_mat <- ReadData("../Data/ref_fetal_cortex_counts.tsv.gz")

meta.data <- read.table("../Data/ref_fetal_cortex_metadata.csv.gz", 
                        header = T, row.names = 1, sep = ",")
ref.annot <- meta.data$Cluster; names(ref.annot) <- rownames(meta.data)
ref.annot <- droplevels(ref.annot[ref.annot %in% c("InCGE", "InMGE")])
table(ref.annot)

raw_counts_mat <- raw_counts_mat[,names(ref.annot)]
ref.obj <- CreateSeuratObject(raw_counts_mat, min.cells = 10)
ref.obj <- NormalizeData(ref.obj)
ref.obj <- FindVariableFeatures(ref.obj)
ref.obj <- ScaleData(ref.obj, features = rownames(ref.obj))
ref.obj <- SetIdent(ref.obj, value = ref.annot)
dim(ref.obj)

## Normalize chromatin dataset
ga.obj <- CreateSeuratObject(unnorm_ga)
ga.obj <- NormalizeData(ga.obj)
ga.obj <- FindVariableFeatures(ga.obj)
ga.obj <- ScaleData(ga.obj, features = rownames(ga.obj))
clusters <- clusters[colnames(ga.obj)]

var.genes <- union(VariableFeatures(ga.obj), VariableFeatures(ref.obj))
var.genes <- var.genes[var.genes %in% rownames(ga.obj) & var.genes %in% rownames(ref.obj)]
length(var.genes)

## Correlate gene expression and gene activity
scale_ga <- GetAssayData(ga.obj, slot = "scale.data")[var.genes,]
scale_rna <- GetAssayData(ref.obj, slot = "scale.data")[var.genes,]
scRNA_centroids <- t(apply(scale_rna[,names(ref.annot)], 1, function(x) tapply(x, ref.annot, mean)))
ga_centroids <- t(apply(scale_ga, 1, function(x) tapply(x, clusters, mean)))

rna.ga.cor <- sapply(colnames(ga_centroids), function(x) {
  sapply(colnames(scRNA_centroids), function(y) cor(ga_centroids[,x], scRNA_centroids[,y], method = "spearman"))
})

pdf("fCTX-combined_In_subclustering_ga_cor.pdf", width = 4, height = 1.5)
ggHeat(rna.ga.cor, clustering = "both", x.lab.size = 12, y.lab.size = 12)
dev.off()

## make umap plot
PlotDims(umap.emb, sample.groups = clusters, x.lab = "umap1", y.lab = "umap2",
         pt.size = 0.5, alpha.plot = 0.4, label.size = 4.5, show.legend = F, show.axes = F,
         seed = 625)

ga.obj <- SetIdent(ga.obj, value = clusters)

# ga.markers.df <- FindAllMarkers(ga.obj, logfc.threshold = 0.1, only.pos = T)
# top.ga.markers.df <- do.call(rbind, by(ga.markers.df, ga.markers.df$cluster, head, n = 10))

marker.genes <- c("GLI3", "HOPX", "CRYAB", "TOP2A", "EOMES")

ga.cl.marker.expr <- apply(GetAssayData(ga.obj, slot = "scale.data")[marker.genes,], 1, function(x) {
  tapply(x, clusters[colnames(ga.obj)], mean)
})
ggHeat(ga.cl.marker.expr, clustering = "none")

DotPlot(ref.obj, features = marker.genes) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cl.mapping <- c("0"="InCGE", "1"="InMGE", "2"="InMGE", 
                "3"=NA, "4"=NA, "5"=NA, "6"=NA)
named.clusters <- plyr::revalue(clusters, replace = cl.mapping)

pdf("fCTX-combined_In_subclustering_named_umap.pdf", width = 5, height = 5)
PlotDims(umap.emb, sample.groups = named.clusters, x.lab = "umap1", y.lab = "umap2",
         pt.size = 0.5, alpha.plot = 0.4, label.size = 4.5, show.legend = F, show.axes = F,
         seed = 625)
dev.off()

pdf("fCTX-combined_In_subclustering_named_umap_nolabels.pdf", width = 5, height = 5)
PlotDims(umap.emb, sample.groups = named.clusters, x.lab = "umap1", y.lab = "umap2",
         pt.size = 0.5, alpha.plot = 0.4, label.size = 0, show.legend = F, show.axes = F,
         seed = 625)
dev.off()

metadata.df <- data.frame(cell = names(clusters), orig.ident = clusters,
                          named.ident = named.clusters,
                          time.point = time.point[names(clusters)])
write.table(metadata.df, file = "../Data/scTHS_In_subclustering_metadata.tsv", sep = "\t")


## Save results
rm(counts, raw_counts_mat, cisTopicObject)
ref.obj@assays$RNA@scale.data <- ga.obj@assays$RNA@scale.data <- matrix(0,0,0)
invisible(gc())
save.image(output.file)
