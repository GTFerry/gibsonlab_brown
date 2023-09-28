library(ArchR)

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
addArchRGenome("hg38")

library(BSgenome.Hsapiens.UCSC.hg38)

archR_project <- loadArchRProject("saves/peak-called")


markersPeaks <- getMarkerFeatures(
  ArchRProj = archR_project, 
  useMatrix = "PeakMatrix", 
  groupBy = "SeuratCluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = archR_project, addDOC = FALSE)

saveArchRProject("saves/identified", load = FALSE, overwrite = TRUE)