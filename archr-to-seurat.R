library(ArchR)
library(Seurat)
library(dplyr)

### Load ArchR Object ###

known_markers <- c(
  "VCAN",
  "CD14",
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
  "PPP1R16B"
)

args = commandArgs(trailingOnly=TRUE)
archR_project_path = args[1]

archR_project <- loadArchRProject("saves/clustered_manual")
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
# 
# sorted_metadata <- metadata[match(colnames(gex_data_subset$nFeature_RNA), metadata$Barcodes),]
# 
# 
# # Identify columns in the metadata that contain the word "Cluster"
cluster_columns <- grep("Cluster", colnames(sorted_metadata), value = TRUE)
# 
# 
# # Remove the 'C' from the beginning of cluster IDs and convert to numeric
for (col_name in cluster_columns) {
  sorted_metadata[[col_name]] <- as.numeric(gsub("C", "", sorted_metadata[[col_name]]))
}
# 
# print("Creating RNA Assay")
# RNA_assay <- CreateAssayObject(counts = gmat)
# 
# print("Creating Seurat Object")
# seurat_obj <- CreateSeuratObject(
#   counts = RNA_assay,
#   assay = "RNA"
# )
# 
# # Initialize a vector to hold mock column names
# mock_names <- vector("character", length = ncol(sorted_metadata))

# Counter to keep track of repetitions
counter <- 1

# Generate mock names
# for (i in 1:26) {
#   for (j in 1:26) {
#     if (counter > ncol(sorted_metadata)) {
#       break
#     }
#     mock_names[counter] <- paste0(letters[i], letters[j], letters[i], letters[j])
#     counter <- counter + 1
#   }
#   if (counter > ncol(sorted_metadata)) {
#     break
#   }
# }

# Assign mock names to columns
# colnames(sorted_metadata) <- mock_names
# 
# all(rownames(sorted_metadata) %in% colnames(seurat_obj))
# sum(is.na(sorted_metadata))


print("Adding Metadata")
# print(dim(gmat))
# print(dim(sorted_metadata))
# print(sapply(sorted_metadata, class))

for (i in cluster_columns) {
  gex_data_subset <- AddMetaData(gex_data_subset, eval(parse(text=paste("sorted_metadata$",i))), i)
}

# seurat_obj <- AddMetaData(object = seurat_obj, metadata = sorted_metadata)
# Make sure the rows of metadata match with the columns of the Seurat object
# if (all(rownames(sorted_metadata) %in% colnames(seurat_obj))) {
#   
#   # Re-order the metadata to match the Seurat object
#   sorted_metadata <- sorted_metadata[match(colnames(seurat_obj), rownames(sorted_metadata)),]
#   
#   # Assign metadata directly
#   seurat_obj@meta.data <- cbind(seurat_obj@meta.data, sorted_metadata)
#   
# } else {
#   stop("Row names in metadata do not match with Seurat object column names.")
# }

# Get names of all available embeddings in the ArchR project
available_embeddings <- names(archR_project@embeddings)

# Filter to only include embeddings with "UMAP_Combined" in the name
umap_combined_embeddings <- grep("UMAP_Combined", available_embeddings, value = TRUE)



# Loop through each UMAP_Combined embedding and add it to the Seurat object
for (embedding_name in umap_combined_embeddings) {

  # Get the UMAP coordinates from the ArchR project
  umap_coords <- getEmbedding(ArchRProj = archR_project,
                              embedding = embedding_name,
                              returnDF = TRUE)

  # Match the order of the UMAP coordinates to the order of the cells in the Seurat object
  # print(match(colnames(gex_data_subset), rownames(umap_coords)),)
  # print(colnames(gex_data_subset$nFeature_RNA))
  #print(rownames(umap_coords))
  umap_coords <- umap_coords[match(colnames(gex_data_subset), rownames(umap_coords)),]

  # Convert to a matrix
  mat_umap_coords <- as.matrix(umap_coords)

  # Create a Seurat DimReduc object
  umap <- CreateDimReducObject(embeddings = mat_umap_coords)

  # Add the DimReduc object to the Seurat object under the appropriate name
  gex_data_subset[[embedding_name]] <- umap

}


# available_embeddings <- names(archR_project@embeddings)
# 
# # Filter to only include embeddings with "UMAP_Combined" in the name
# umap_combined_embeddings <- grep("UMAP_Combined", available_embeddings, value = TRUE)
# 
# # Get cell names from Seurat object
# seurat_cell_names <- colnames(gex_data_subset)
# 
# # Loop through each UMAP_Combined embedding and add it to the Seurat object
# for (embedding_name in umap_combined_embeddings) {
#   
#   # Get the UMAP coordinates from the ArchR project
#   umap_coords <- getEmbedding(ArchRProj = archR_project, 
#                               embedding = embedding_name, 
#                               returnDF = TRUE)
#   
#   # Match the order of the UMAP coordinates to the order of the cells in the Seurat object
#   matching_indices <- match(seurat_cell_names, rownames(umap_coords))
#   
#   if (any(is.na(matching_indices))) {
#     warning(paste("There are", sum(is.na(matching_indices)), "cells in Seurat object not found in the embedding", embedding_name))
#   }
#   
#   # Get the matched UMAP coordinates
#   umap_coords_matched <- umap_coords[matching_indices, , drop = FALSE]
#   
#   # Add the UMAP coordinates to the metadata of the Seurat object
#   gex_data_subset[["umap_coords1"]] <- umap_coords_matched[, 1]
#   gex_data_subset[["umap_coords2"]] <- umap_coords_matched[, 2]
# }


pdf("dotplots.pdf", width = 16, height = 8)
for (i in cluster_columns) {
  print(paste("Saving dotplot: ", i))
  gex_data_subset <- SetIdent(gex_data_subset, value = i)
  
  # DimPlot(gex_data_subset, reduction = <imagine that this is right>)

  levels(gex_data_subset) <- as.character(sort(as.numeric(levels(gex_data_subset))))
  
  print(DotPlot(gex_data_subset, features = known_markers, scale = FALSE) +
    ggtitle(i) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
}
dev.off()

