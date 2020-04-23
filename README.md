# chromfunks
Helper functions for chromatin accessibility analysis

Installation instructions:

```
BiocManager::install()
BiocManager::install(c("limma", "GenomicRanges", "biomaRt", "ChIPseeker", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"))
devtools::install_github("yanwu2014/chromfunks")
```

## Links to tutorials
A [walkthrough](https://yanwu2014.github.io/chromfunks/Tutorials/gChromVAR_example.html) on how to use [CausalDB](http://mulinlab.org/causaldb/index.html) and [gChromVAR](https://github.com/caleblareau/gchromVAR) to link GWAS traits to cell types using single-cell or bulk chromatin accessibility data.
