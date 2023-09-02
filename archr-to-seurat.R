library(ArchR)
library(Seurat)

### Load ArchR Object ###

args = commandArgs(trailingOnly=TRUE)
archR_project_path = args[1]

archR_project <- loadArchRProject(archR_project_path)

print("Loaded archr project")

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

print("Getting GeneExpressionMatrix")

raw_gex_mat <- getMatrixFromProject(
                                    ArchRProj = archR_project,
                                    useMatrix = "GeneExpressionMatrix",
)

gmat <- raw_gex_mat@assays@data$GeneExpressionMatrix
colnames(gmat) <- colnames(raw_gex_mat)
rownames(gmat) <- raw_gex_mat@elementMetadata$name

print("Getting Metadata")
metadata <- getCellColData(
                           ArchRProj = archR_project
)

metadata$Barcodes <- rownames(metadata)

sorted_metadata <- metadata[match(colnames(gmat), metadata$Barcodes),]


# Identify columns in the metadata that contain the word "Cluster"
cluster_columns <- grep("Cluster", colnames(sorted_metadata), value = TRUE)


# Remove the 'C' from the beginning of cluster IDs and convert to numeric
for (col_name in cluster_columns) {
  sorted_metadata[[col_name]] <- as.numeric(gsub("C", "", sorted_metadata[[col_name]]))
}

print("Creating RNA Assay")
RNA_assay <- CreateAssayObject(counts = gmat)

print("Creating Seurat Object")
seurat_obj <- CreateSeuratObject(
                                 counts = RNA_assay,
                                 assay = "RNA"
)

# Initialize a vector to hold mock column names
# mock_names <- vector("character", length = ncol(sorted_metadata))
# 
# # Counter to keep track of repetitions
# counter <- 1
# 
# # Generate mock names
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
# 
# # Assign mock names to columns
# colnames(sorted_metadata) <- mock_names

all(rownames(sorted_metadata) %in% colnames(seurat_obj))
sum(is.na(sorted_metadata))


print("Adding Metadata")
print(dim(gmat))
print(dim(sorted_metadata))
print(sapply(sorted_metadata, class))



# seurat_obj <- AddMetaData(object = seurat_obj, metadata = sorted_metadata)
# Make sure the rows of metadata match with the columns of the Seurat object
if (all(rownames(sorted_metadata) %in% colnames(seurat_obj))) {
  
  # Re-order the metadata to match the Seurat object
  sorted_metadata <- sorted_metadata[match(colnames(seurat_obj), rownames(sorted_metadata)),]
  
  # Assign metadata directly
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, sorted_metadata)
  
} else {
  stop("Row names in metadata do not match with Seurat object column names.")
}

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
  umap_coords <- umap_coords[match(colnames(gmat), rownames(umap_coords)),]
  
  # Convert to a matrix
  mat_umap_coords <- as.matrix(umap_coords)
  
  # Create a Seurat DimReduc object
  umap <- CreateDimReducObject(embeddings = mat_umap_coords)
  
  # Add the DimReduc object to the Seurat object under the appropriate name
  seurat_obj[[embedding_name]] <- umap
  
}


# save seurat object
saveRDS(seurat_obj, file = "seurat_obj.rds")
