#!/usr/bin/env Rscript

library(ArchR)
library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(Matrix)
library(Seurat)
library(stringr)

# 11-21-2022
# Get matrices for predicting GEX data:
# 1.) Peak matrix
# 2.) LogNorm GEX data

addArchRThreads(threads = 2)
addArchRGenome("hg38")

projIBD <- loadArchRProject("Singletons-Final-projIBD-ATACandGEX-11-2022")

# 1.) Get Peak Matrix:

raw_peak_mat <- getMatrixFromProject(
  ArchRProj = projIBD,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

pmat <- raw_peak_mat@assays@data$PeakMatrix
colnames(pmat) <- colnames(raw_peak_mat)

# chr
chr <- raw_peak_mat@rowRanges@seqnames
# region:
start <- data.frame(raw_peak_mat@rowRanges@ranges)$start
end <- data.frame(raw_peak_mat@rowRanges@ranges)$end

rownames(pmat) <- paste(chr, paste(start, end, sep = "-"), sep = ":")

write(colnames(pmat), file = "/storage/home/mfisher42/scProjects/Predict_GEX/input_data/sparse_peak_matrix_colnames.txt")
write(rownames(pmat), file = "/storage/home/mfisher42/scProjects/Predict_GEX/input_data/sparse_peak_matrix_rownames.txt")
writeMM(pmat, file = "/storage/home/mfisher42/scProjects/Predict_GEX/input_data/sparse_peak_matrix.txt")

# 1.) Get Gex Matrix:

raw_gex_mat <- getMatrixFromProject(
  ArchRProj = projIBD,
  useMatrix = "GeneExpressionMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

gmat <- raw_gex_mat@assays@data$GeneExpressionMatrix
colnames(gmat) <- colnames(raw_gex_mat)

# gene names
rownames(gmat) <- raw_gex_mat@elementMetadata$name

# gene regions
chr <- raw_gex_mat@elementMetadata$seqnames
start <- raw_gex_mat@elementMetadata$start
end <- raw_gex_mat@elementMetadata$end
region <- paste(chr, paste(start, end, sep = "-"), sep = ":")
gene_region <- data.frame(raw_gex_mat@elementMetadata$name, region)

write.table(gene_region, file = "saves/peak-calling-output/script-output/parse_gex_matrix_regions.txt")
write(colnames(gmat), file = "saves/peak-calling-output/script-output/sparse_gex_matrix_colnames.txt")
write(rownames(gmat), file = "saves/peak-calling-output/script-output/sparse_gex_matrix_rownames.txt")
writeMM(gmat, file = "saves/peak-calling-output/script-output/parse_gex_matrix.txt")



