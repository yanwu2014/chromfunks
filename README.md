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


## Reproducing Analysis

This repository also contains the scripts necessary to reproduce the analysis in [Wu et al]().

1. Clone this repo: `git clone https://github.com/yanwu2014/chromfunks.git`

2. Enter the scripts directory: `cd chromfunks/Scripts/`

3. Download the raw/processed data from [our FTP server](ftp://genome-miner.ucsd.edu/fetal_brain_data/Data/): `wget -m ftp://genome-miner.ucsd.edu/fetal_brain_data/Data/`

4. If you're just interested in running the higher level transcription factor analysis set your working directory to the `TF_Analysis/` sub-directory and run those scripts. We include the annotated Seurat objects for the scTHS, snDrop, sciATAC, and sciRNA datasets necessary to run these scripts.

5. If you're interested in generating these Seurat objects from their respective counts matrices, run the R scripts within each sub-directory to reproduce the analysis. The scripts need to be run in a certain order reflected in their numbering. Be sure the set the working directory to the same sub-directory that the script is in before running it.
