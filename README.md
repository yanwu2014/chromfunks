# chromfunks
A R package designed to help with chromatin accessibility analysis by building regulatory networks centered around noncoding accessible regions

Installation instructions:

```
BiocManager::install()
BiocManager::install(c("limma", "GenomicRanges", "biomaRt", "ChIPseeker", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "motifmatchr", "chromVAR", "SummarizedExperiment", 
"BSgenome.Hsapiens.UCSC.hg38"))
devtools::install_github("yanwu2014/chromfunks")
```

## Tutorials

### Creating transcription factor networks centered around noncoding regulatory regions ####
A [walkthrough](https://yanwu2014.github.io/chromfunks/Tutorials/TF_networks_example.html) on how to build transcription factor (TF) networks centered on noncoding regulatory regions. In these networks, TFs regulate genes via binding to one or more noncoding regulatory regions, thus adding biological context to the functions of these regions. Here we look at the regulatory context around Human Accelerated Regions (HARs) using single cell chromatin accessibility and gene expression datasets from the fetal human cortex and a Hi-C dataset from human Radial Glia.


### Using CausalDB and g-chromVAR to link GWAS traits to cell types using single cell chromatin accessibility data
A [walkthrough](https://yanwu2014.github.io/chromfunks/Tutorials/gChromVAR_example.html) on how to use [CausalDB](http://mulinlab.org/causaldb/index.html) and [gChromVAR](https://github.com/caleblareau/gchromVAR) to link GWAS traits to cell types using single-cell or bulk chromatin accessibility data.
