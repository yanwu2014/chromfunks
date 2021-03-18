library(methods)
library(cisTopic)
library(cicero)
library(Seurat)
library(Signac)
library(swne)
library(Matrix)
library(chromfunks)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

## Output file
output.file <- "../Data/adult_scTHS_cisTopic_results.RData"

## Load accessibility counts
# load("adult-vCTX_SnapATAC_counts.RData")
# colnames(counts) <- gsub("@", "", colnames(counts))
# colnames(counts) <- gsub("#SORTED", "", colnames(counts))
# saveRDS(counts, file = "adult_THS_counts.Robj")
counts <- readRDS("../Data/adult_scTHS_counts.Robj")

## Load metadata
metadata.df <- read.table("../Data/adult_scTHS_metadata.tsv.gz", 
                          header = T, stringsAsFactors = F)
metadata.df$cell <- sapply(metadata.df$cell, function(x) {
  fields <- strsplit(x, split = "_")[[1]]
  paste0(fields[[1]], fields[[2]], fields[[3]])
})
orig.ident <- factor(metadata.df$named.ident)
names(orig.ident) <- metadata.df$cell
table(orig.ident)

## Load fraction of unique reads in peaks
nReads <- Matrix::colSums(counts)

## Binarize matrix
counts@x[counts@x > 1] <- 1

## Filter for cells with enough peaks and enough reads in peaks
counts <- counts[,colnames(counts) %in% names(orig.ident)]

## reformatting matrix rownames
peaks <- peak2granges(rownames(counts))
chrom <- seqnames(peaks)
loc_start <- start(peaks)
loc_end <- end(peaks)

## initialize cisTopic object from count matrix
cisTopicObject <- createcisTopicObject(counts, project.name = 'fCTX', keepCountsMatrix = F)

## run LDA model
topics.range <- seq(20,30,5)
isTopicObject <- runCGSModels(cisTopicObject, topic = topics.range, seed = 987, nCores = length(topics.range),
                              burnin = 120, iterations = 250, addModels = F)

## select for model
cisTopicObject <- cisTopic::selectModel(cisTopicObject)
# logLikelihoodByIter(cisTopicObject, select = topics.range)

## run UMAP
cfg <- umap::umap.defaults; cfg$metric <- "cosine"; cfg$min_dist <- 0.2;
cisTopicObject <- runUmap(cisTopicObject, target = 'cell', method = 'Z-score', config = cfg)

# Save image
save.image(output.file)

## pull out umap coordinates
umap.emb <- cisTopicObject@dr$cell[["Umap"]]

## pull out topic coordinates
topic.emb <- t(modelMatSelection(cisTopicObject, target = "cell", method = 'Z-score'))

## pull out pca coordinates
pc.emb <- cisTopicObject@dr$cell[["PCA"]]$ind.coord
orig.ident <- orig.ident[names(orig.ident) %in% colnames(counts)]

## make umap plot
# pdf("adult-vCTX_clusters_orig_ident.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = orig.ident, pt.size = 0.3, alpha.plot = 0.4, label.size = 4.5,
         show.legend = F, show.axes = F, seed = 625)
# dev.off()

non.missing.ident <- orig.ident[!is.na(orig.ident)]
non.missing.ident <- non.missing.ident[names(non.missing.ident) %in% colnames(cisTopicObject@binary.count.matrix)]
# diff_peaks_list <- find_all_diff(cisTopicObject@binary.count.matrix[,names(non.missing.ident)], non.missing.ident)
# save(diff_peaks_list, file = "adult-vCTX_diff_peaks_cisTopic.RData")

## Run cicero
hg38.chr.lengths <- read.table("../Data/hg38.chr.lengths.txt", header = F, sep = "\t")
hg38.chr.lengths[[2]] <- as.numeric(hg38.chr.lengths[[2]])

orig.ident <- orig.ident[colnames(cisTopicObject@binary.count.matrix)]
pData <- data.frame(orig.ident)
rownames(pData) <- colnames(cisTopicObject@binary.count.matrix)

peaks.use <- rownames(cisTopicObject@binary.count.matrix)
fData <- data.frame(site_name = peaks.use, chromosome = chrom[peaks.use], 
                    bp1 = loc_start[peaks.use], bp2 = loc_end[peaks.use])
rownames(fData) <- peaks.use
input_cds <- newCellDataSet(cisTopicObject@binary.count.matrix[,names(orig.ident)], 
                            phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.emb, k = 50)
print("Done making cicero cds")

conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run
print("Done finding connections")
save.image(output.file)

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
# write.table(peak_anno_df, file = "../Data/adult_scTHS_peaks_filtered.bed", sep = "\t", quote = F,
#             row.names = F, col.names = F)

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

## Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.1)
print(dim(unnorm_ga))


## Build and export Seurat object
chrom_assay <- CreateChromatinAssay(
  counts = cisTopicObject@binary.count.matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 0,
  min.features = 0
)

adult.ths <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

## Add UMAP
adult.ths[["umap"]] <- CreateDimReducObject(embeddings = umap.emb, key = "UMAP_", assay = "peaks")

## Add clusters
adult.ths$ident <- orig.ident[colnames(adult.ths)]
table(adult.ths$ident)

## Add cisTopic topics
cisTopicObject <- getRegionsScores(cisTopicObject)
cisTopicObject <- binarizecisTopics(cisTopicObject)
topic.emb <- modelMatSelection(cisTopicObject, "cell", "Probability")
rownames(topic.emb) <- paste0("topic_", 1:nrow(topic.emb))
# topic.load <- modelMatSelection(cisTopicObject, "region", "Probability")
adult.ths[["topics"]] <- CreateDimReducObject(embeddings = t(topic.emb), key = "topic_", assay = "peaks")

ga_assay <- CreateAssayObject(unnorm_ga[,cisTopicObject@cell.names])
adult.ths[["RNA"]] <- ga_assay

## Write scTHS data to file
saveRDS(adult.ths, file = "../Data/adult_scTHS_chromatin_seurat_object.Robj")

# Save image
save.image(output.file)
