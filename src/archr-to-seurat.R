library(ArchR)
library(Seurat)
library(dplyr)

### Load ArchR Object ###

# Generic
known_markers <- c(
  "VCAN",
  "CD14",
  "CD16",
  "LYZ",
  "FCGR3A",
  "MS4A7",
  "CD86",
  "ITGAX",
  "FCER1A",
  "CST3",
  "CD1C",
  "CLEC4C",
  "MS4A1",
  "TCL1A",
  "PAX5",
  "CD24",
  "CD83",
  "CD27",
  "PRDM1",
  "IRF4",
  "XBP1",
  "CD3D",
  "CD4",
  "CD8A",
  "IL7R",
  "CCR7",
  "SELL",
  "ANK3",
  "KLRC2",
  "ZEB1",
  "STAT3",
  "IL21R",
  "STAT1",
  "ISG15",
  "IL12RB2",
  "GZMH",
  "GZMB",
  "CCL5",
  "CD69",
  "KLRB1",
  "CCR5",
  "CCR6",
  "CXCR6",
  "IL18R1",
  "GNLY",
  "NKG7",
  "TGFB1",
  "ABCB1",
  "PPP1R16B",
  "MKI67",
  "PPBP"
)

# B Cell
# known_markers <- c(known_markers, c(
# "FOXP3",
# "STAT5A",
# "CTLA4",
# "PECAM1",
# "CCR7",
# "IL2RA",
# "CCR5",
# "TBX21",
# "GZMA",
# "CD69",
# "ITGAE",
# "IL7",
# "IL7R",
# "SELL",
# "CD27",
# "CD3D",
# "CCR7",
# "ANK3",
# "ZEB1"
# #"HLA-DR",
# ))

# T Cell
# known_markers <- c(known_markers,
#   "CD4",
#   "CD8A",
#   "CD8B",
#   "TXNIP",
#   "BTG1",
#   "ZFP36L2",
#   "ANK3",
#   "NABP1",
#   "IL7R",
#   "CTLA4",
#   "IL21",
#   "CD200",
#   "CXCL13",
#   "PDCD1",
#   "FKBP5",
#   "CXCR5",
#   "BCL6",
#   "SELL",
#   "KLF2",
#   "ITGA6",
#   "LEF1",
#   "CCR7",
#   "S100A4"
# )

# Monocytes

# known_markers <- c(
#   "S100A8",
#   "CD163",
#   "LYZ",
#   "CD14",
#   "NKG7",
#   "GNLY",
#   "GZMH",
#   "IL7R",
#   "CD69",
#   "MAL",
#   "MS4A1",
#   "CD79A",
#   "CD79B",
#   "FCGR3A",
#   "HES1",
#   "HES4",
#   "PF4",
#   "PPBP",
#   "HLAs",
#   "CD1C",
#   "PTGDS",
#   "ITM2C",
#   "CCDC50",
#   "NEATTAOK1"
# )

known_markers <- unique(known_markers)

args = commandArgs(trailingOnly=TRUE)
archR_project_path = args[1]

archR_project <- loadArchRProject("saves/clustered")
gex_data <- Read10X_h5("data/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5")
seurat_project <- CreateSeuratObject(counts = gex_data$`Gene Expression`)

qc_cells <- archR_project@cellColData@rownames
qc_cells_seurat <- substring(qc_cells, 10)
seurat_cell_ids <- gex_data[["Gene Expression"]]@Dimnames[[2]]
common_cell_ids <- intersect(seurat_cell_ids, qc_cells_seurat)

gex_data_subset <- subset(seurat_project, cells = common_cell_ids)

gex_data_subset <- RenameCells(gex_data_subset, new.names = paste0("PBMC_10k#", colnames(gex_data_subset)))

gex_data_subset <- NormalizeData(gex_data_subset, normalization.method = "LogNormalize", scale.factor = 10000)
gex_data_subset <- ScaleData(gex_data_subset, features = known_markers)

rm(gex_data)

addArchRThreads(threads = 6)
addArchRGenome("hg38")

# raw_peak_mat <- getMatrixFromProject(
#   ArchRProj = archR_project,
#   useMatrix = "TileMatrix"
# )
# 
# pmat <- raw_peak_mat@assays@data$PeakMatrix
# colnames(pmat) <- colnames(raw_peak_mat)
# 
# chr <- raw_peak_mat@elementMetadata@listData[["seqnames"]]
# start <- data.frame(raw_peak_mat@rowRanges@ranges)$start
# end <- data.frame(raw_peak_mat@rowRanges@ranges)$end
# 
# rownames(pmat) <- paste(chr, paste(start, end, sep = "-"), sep = ":")

# get RNA matrix

# print("Getting GeneExpressionMatrix")
# 
# raw_gex_mat <- getMatrixFromProject(
#   ArchRProj = archR_project,
#   useMatrix = "GeneExpressionMatrix",
# )
# 
# gmat <- raw_gex_mat@assays@data$GeneExpressionMatrix
# colnames(gmat) <- colnames(raw_gex_mat)
# rownames(gmat) <- raw_gex_mat@elementMetadata$name
# 
# print("Getting Metadata")
metadata <- getCellColData(
  ArchRProj = archR_project
)

metadata$Barcodes <- rownames(metadata)
sorted_metadata <- metadata[match(colnames(gex_data_subset), metadata$Barcodes),]

cluster_columns <- grep("Cluster", colnames(sorted_metadata), value = TRUE)

for (col_name in cluster_columns) {
  sorted_metadata[[col_name]] <- as.numeric(gsub("C", "", sorted_metadata[[col_name]]))
}

print("Adding Metadata")

for (i in cluster_columns) {
  gex_data_subset <- AddMetaData(gex_data_subset, eval(parse(text=paste("sorted_metadata$",i))), i)
}

available_embeddings <- names(archR_project@embeddings)

# Filter to only include embeddings with "UMAP_Combined" in the name
umap_combined_embeddings <- grep("UMAP_Combined", available_embeddings, value = TRUE)



# Loop through each UMAP_Combined embedding and add it to the Seurat object
for (embedding_name in umap_combined_embeddings) {

  # Get the UMAP coordinates from the ArchR project
  umap_coords <- getEmbedding(ArchRProj = archR_project,
                              embedding = embedding_name,
                              returnDF = TRUE)

  umap_coords <- umap_coords[match(colnames(gex_data_subset), rownames(umap_coords)),]

  # Convert to a matrix
  mat_umap_coords <- as.matrix(umap_coords)

  umap <- CreateDimReducObject(embeddings = mat_umap_coords)

  # Add the DimReduc object to the Seurat object under the appropriate name
  gex_data_subset[[embedding_name]] <- umap
  
  # gex_data_subset@reductions[[embedding_name]]@key <- "UMAP_"

}


counter <- 1

pdf("dotplots.pdf", width = 16, height = 8)
for (i in cluster_columns) {
  print(paste("Saving dotplot: ", i))
  gex_data_subset <- SetIdent(gex_data_subset, value = i)

  levels(gex_data_subset) <- as.character(sort(as.numeric(levels(gex_data_subset))))
  clusters_to_plot <- c(1:4, 18)
  
  #### To Rename Clusters (18)
  new.cluster.ids <- c(
   "Non-classical Monocytes", # 1
   "Unknown", # 2
   "Conventional Dendritic Cells", # 3
   "FCGR3A+ Monocytes", # 4
   "CD14+ Monocytes", # 5
   "CD14+ Monocytes", # 6
   "CD8 Naive T Cells", # 7
   "CD4 Naive T Cells", # 8
   "CD8 Memory T Cells", # 9
   "CD4 Memory T Cells", # 10
   "MAIT Cells", # 11
   "NK Cells", # 12
   "CD8+ Cytotoxic T Cells", # 13
   "Plasmacytoid Dendritic Cells", # 14
   "Plasma B Cells", # 15
   "Activated B Cells", # 16
   "Resting Naive B Cells", # 17
   "Transitional B Cells" # 18
   )

  names(new.cluster.ids) <- levels(gex_data_subset)
  # names(new.cluster.ids) <- rev(c(
  #   5, 6, 4, 1, 2, 3, 14,
  #   18, 17, 16, 15,
  #   7, 8, 9, 10, 13,
  #   11, 12
  # ))
  
  gex_data_subset <- RenameIdents(gex_data_subset, new.cluster.ids)
  
  monocyteComparison <- FindMarkers(gex_data_subset,
                           ident.1 = "Unknown",
                           ident.2 = "Non-classical Monocytes", ident.3 = "FCGR3A+ Monocytes", ident.4 = "CD14+ Monocytes")
  monocyteComparison <- monocyteComparison[monocyteComparison$p_val_adj < 0.05,]
  monocyteComparison <- monocyteComparison[monocyteComparison$avg_log2FC > 0,]
  
  monocyteComparison <- monocyteComparison %>% arrange(desc(avg_log2FC))# %>% slice_head(n = 15)
  
  bcellComparison <- FindMarkers(gex_data_subset,
                           ident.1 = "Unknown",
                           ident.2 = "Plasma B Cells", ident.3 = "Activated B Cells", ident.4 = "Resting Naive B Cells", ident.4 = "Transitional B Cells")
  bcellComparison <- bcellComparison[bcellComparison$p_val_adj < 0.05,]
  bcellComparison <- bcellComparison[bcellComparison$avg_log2FC > 0,]
  
  bcellComparison <- bcellComparison %>% arrange(desc(avg_log2FC))#  %>% slice_head(n = 15)
  
  # Subset gex_data_subset to only include the clusters you want to plot
  # gex_data_subset <- subset(gex_data_subset, idents = new.cluster.ids[clusters_to_plot])
  ###
  
  print(DotPlot(gex_data_subset, features = known_markers, scale = FALSE) +
    ggtitle(i) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
  print(DimPlot(gex_data_subset, reduction = paste0("UMAP_Combined_", counter), label = TRUE, label.size = 5)) + ggtitle(paste0(i, "  - UMAP"))
  counter <- counter + 1
  
  # diff_markers <- FindMarkers(gex_data_subset, ident.1 = "Trm", ident.2 = "Tcm", ident.3 = "GdT Cells", ident.4 = "Th1/Th17 Cells", min.pct = 0.25)
  # write.csv(diff_markers, "diff.csv")
  
  known_markers_length <- length(known_markers)
  group_size <- 9
  
  for (i in seq(1, known_markers_length, by = group_size)) {
    end_index <- min(i + group_size - 1, known_markers_length)
    print(FeaturePlot(gex_data_subset, features = known_markers[c(i:end_index)]))
  }
  
  print(VlnPlot(gex_data_subset, features = "nFeature_RNA", idents = new.cluster.ids))  
  
  print("Monocyte comparison, p < 0.05, top 15")
  for (i in rownames(monocyteComparison)) {cat(c(i, "\n"))}
  print("B Cell comparison, p < 0.05, top 15")
  for (i in rownames(bcellComparison)) {cat(c(i, "\n"))}
}
dev.off()

seurat_clusters <- data.frame(gex_data_subset@active.ident)

sorted_metadata <- seurat_clusters[match(archR_project$cellNames, rownames(seurat_clusters)),]

sorted_metadata <- as.character(sorted_metadata)

# Add the cluster labels to the ArchR object
archR_project <- addCellColData(
  ArchRProj = archR_project, 
  data = sorted_metadata,
  cells = archR_project$cellNames,
  name = "SeuratCluster",
  force = TRUE
)


plotEmbedding(
  ArchRProj = archR_project, 
  colorBy = "cellColData",
  name = "SeuratCluster",
  embedding = "UMAP_Combined_1"
)

saveArchRProject(archR_project, "saves/labeled", load = FALSE, overwrite = TRUE)