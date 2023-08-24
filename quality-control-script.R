# This is a scriptified version of the notebook for use on the server

library(ArchR)
library(Seurat)
library(Dune)

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("Dune")
# library("BSgenome.Hsapiens.UCSC.hg38")

addArchRGenome("hg38")

inputFiles <- "data/pbmc_unsorted_10k_atac_fragments.tsv.gz"
names(inputFiles) <- "PBMC_10k"

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  minTSS = 0,
  minFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
  #force = TRUE
)

proj <- ArchRProject(ArrowFiles)

seRNA <- import10xFeatureMatrix(
  input = c("data/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5"),
  names = c("PBMC_10k")
)

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

minFeatures <- 600
maxFeatures <- 4000
maxMT <- 0.019

# Filter out cells with NA values for the relevant metrics
proj <- proj[!is.na(proj$Gex_nGenes) & 
               !is.na(proj$Gex_MitoRatio)]

proj <- proj[proj$Gex_nGenes > minFeatures
             & proj$Gex_nGenes < maxFeatures
             & proj$Gex_MitoRatio < maxMT]

minTSS <- 12
minFrags <- 1200

proj <- proj[proj$TSSEnrichment > minTSS & proj$nFrags > minFrags]

saveArchRProject(ArchRProj = proj, outputDirectory = "saves/atac-rna-quality-controlled", load = FALSE)