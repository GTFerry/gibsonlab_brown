library(ArchR)
library(Dune)

### Load ArchR Object ###

args = commandArgs(trailingOnly=TRUE)
archR_project_path = args[1]

archR_project <- loadArchRProject(archR_project_path)

print("Loaded archr project")

addArchRThreads(threads = 6)
addArchRGenome("hg38")

raw_peak_mat <- getMatrixFromProject(
  ArchRProj = archR_project,
  useMatrix = "TileMatrix"
)

pmat <- raw_peak_mat@assays@data$PeakMatrix
colnames(pmat) <- colnames(raw_peak_mat)

chr <- raw_peak_mat@rowRanges@seqnames
start <- data.frame(raw_peak_mat@rowRanges@ranges)$start
end <- data.frame(raw_peak_mat@rowRanges@ranges)$end

rownames(pmat) <- paste(chr, paste(start, end, sep = "-"), sep = ":")

# get RNA matrix

raw_gex_mat <- getMatrixFromProject(
  ArchRProj = archR_project,
  useMatrix = "GeneExpressionMatrix",
)

gmat <- raw_gex_mat@assays@data$GeneExpressionMatrix
colnames(gmat) <- colnames(raw_gex_mat)
rownames(gmat) <- raw_gex_mat@elementMetadata$name
