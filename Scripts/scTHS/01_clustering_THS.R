library(methods)
library(cisTopic)
library(Seurat)
library(swne)
library(cicero)
library(chromfunks)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


## load RData
counts <- readRDS("../Data/scTHS_counts.Robj")

## Load fraction of unique reads in peaks
nReads <- Matrix::colSums(counts)

## Binarize matrix
counts@x[counts@x > 1] <- 1

## Filter for cells with enough peaks and enough reads in peaks
min.peaks <- 500
max.peaks <- 2e4

nPeaks <- Matrix::colSums(counts)
cells.idx <- nPeaks > min.peaks & nPeaks < max.peaks
paste0("Cells to keep: ", sum(cells.idx))
counts <- counts[,cells.idx]

# ## Filter for common peaks
# min.cells <- round(0.001*ncol(counts))
# peaks.idx <- Matrix::rowSums(counts) > min.cells
# paste0("Peaks to keep: ", sum(peaks.idx))
# counts <- counts[peaks.idx,]

## Get cell metadata
file.names <- colnames(counts)
time.point <- factor(sapply(file.names, ExtractField, field = 3, delim = "#"))
names(time.point) <- file.names; levels(time.point) <- c("wk16", "wk18");

## initialize cisTopic object from count matrix
cisTopicObject <- createcisTopicObject(counts, project.name = 'fCTX', keepCountsMatrix = F)

## run LDA model
topics.range <- c(20,25,30,35)
cisTopicObject <- runCGSModels(cisTopicObject, topic = topics.range, seed = 987, nCores = length(topics.range),
                               burnin = 120, iterations = 250, addModels = F)

## select for model
cisTopicObject <- cisTopic::selectModel(cisTopicObject)

# ## check likelihood stablization
# logLikelihoodByIter(cisTopicObject)

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
seurat_obj <- FindClusters(seurat_obj, resolution = 1, verbose = F)
clusters <- Idents(seurat_obj)
table(clusters)

## make umap plot
# pdf("fCTX-combined_umap_clusters_cisTopic.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = clusters, x.lab = "umap1", y.lab = "umap2",
         pt.size = 0.3, alpha.plot = 0.4, label.size = 4.5, show.legend = F, show.axes = F,
         seed = 625)
# dev.off()

# pdf("fCTX-combined_umap_timepoint_cisTopic.pdf", width = 6.5, height = 6)
PlotDims(umap.emb, sample.groups = time.point, x.lab = "umap1", y.lab = "umap2",
         pt.size = 0.1, alpha.plot = 0.3, label.size = 4.5, show.legend = T, show.axes = F,
         seed = 225, do.label = F)
# dev.off()


## Run cicero
hg38.chr.lengths <- read.table("../Data/hg38.chr.lengths.txt", header = F, sep = "\t")

pData <- data.frame(time.point, clusters)
pData <- new("AnnotatedDataFrame", data = pData)

chrom <- sapply(rownames(counts), ExtractField, field = 1, delim = ":")
loc <- sapply(rownames(counts), ExtractField, field = 2, delim = ":")
loc_start <- as.integer(sapply(loc, ExtractField, field = 1, delim = "-"))
loc_end <- as.integer(sapply(loc, ExtractField, field = 2, delim = "-"))
fData <- data.frame(site_name = rownames(counts), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
rownames(fData) <- rownames(counts)
fData <- new("AnnotatedDataFrame", data = fData)

input_cds <- newCellDataSet(counts, phenoData = pData, featureData = fData,
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.emb, k = 80)
print("Done making cicero cds")

conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run
print("Done finding connections")


## Format peak annotation dataframe
peaks.gr <- makeGRangesFromDataFrame(data.frame(seqnames = chrom, start = loc_start, end = loc_end))
peak_anno <- annotatePeak(peaks.gr, tssRegion = c(-3000, 3000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb = "org.Hs.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df[!grepl("Promoter|Exon|Intron", peak_anno_df$annotation), "SYMBOL"] <- NA

peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]
head(peak_anno_df)

# ## Write peak annotations to bed file
# write.table(peak_anno_df, file = "../Data/scTHS_peaks_filtered.bed", sep = "\t",
#             row.names = F, col.names = F)

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.1)
print(dim(unnorm_ga))

## Export gene activity matrix to assign cell types
seurat_obj@assays$RNA@counts <- seurat_obj@assays$RNA@data <- matrix(0)
save(unnorm_ga, seurat_obj, cisTopicObject, clusters, time.point, umap.emb,
     file = "../Data/scTHS_cisTopic_results.RData")

## Export cicero conns for visualization
save(conns, file = "../Data/scTHS_conns.RData")

