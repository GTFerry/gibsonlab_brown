library(ArchR)

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
addArchRGenome("hg38")

library(BSgenome.Hsapiens.UCSC.hg38)