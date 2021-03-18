library(Seurat)
library(chromfunks)
library(Matrix)
library(TFBSTools)
library(SummarizedExperiment)
library(chromVAR)
library(motifmatchr)
library(ChIPseeker)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Repitools)
library(swne)
library(liger)

## Output file
output.file <- "../Data/human_chimp_dA_motif_enrich.RData"

## Get human/chimp dA peaks
hc_peaks <- readPeakFile("../Data/human_chimp_dA_peaks_hg38.bed")
names(hc_peaks) <- paste0(seqnames(hc_peaks), ":",
                          start(hc_peaks), "-",
                          end(hc_peaks))

human_peaks <- hc_peaks[hc_peaks$V4 == "human"]
chimp_peaks <- hc_peaks[hc_peaks$V4 == "chimp"]


## Get background peaks
bg.peaks <- c(human_peaks, chimp_peaks)
names(bg.peaks) <- paste0(seqnames(bg.peaks), ":",
                          start(bg.peaks), "-",
                          end(bg.peaks))

gc.frac <- gcContentCalc(bg.peaks, organism = BSgenome.Hsapiens.UCSC.hg38, verbose = T)
gc.frac.df <- data.frame(GC.percent = gc.frac, Peak = names(bg.peaks))
rownames(gc.frac.df) <- gc.frac.df$Peak

## Load JASPAR motifs
motifs <- getJasparMotifs()

## Get motif matrix
motif_ix <- matchMotifs(motifs, bg.peaks, genome = BSgenome.Hsapiens.UCSC.hg38)
motif_mat <- assay(motif_ix)
rownames(motif_mat) <- names(bg.peaks)
colnames(motif_mat) <- sapply(colnames(motif_mat), ExtractField, field = 2)

## Motif enrichment
human.motif <- MotifEnrich(gc.frac.df, names(human_peaks), motif_mat, n.background = 10000)
chimp.motif <- MotifEnrich(gc.frac.df, names(chimp_peaks), motif_mat, n.background = 10000)
head(human.motif)
head(chimp.motif)

## Save results
save(human.motif, chimp.motif, file = output.file)
