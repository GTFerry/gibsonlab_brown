library(ArchR)
library(Seurat)
library(Dune)

proj <- loadArchRProject(path = "saves/atac-rna-quality-controlled")
results_dir <- "results/"

# Add ATAC LSI
proj <- addIterativeLSI(
  ArchRProj = proj,
  iterations = 1,
  clusterParams = list(
    resolution = 0.75,
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  dimsToUse = 1:3,
  name = "manual_0.75_atac",
  force = TRUE
)

# Add RNA LSI
proj <- addIterativeLSI(
  ArchRProj = proj,
  iterations = 1,
  clusterParams = list(
    resolution = 0.75,
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  dimsToUse = 1:3,
  name = "manual_0.75_rna",
  force = TRUE
)

# Combined Dims
proj <- addCombinedDims(
  proj,
  reducedDims = c("manual_0.75_atac",
                  "manual_0.75_rna"),
  name = "combined_lsi_manual"
)

# Add Combined UMAP
proj <- addUMAP(
  proj,
  reducedDims = "combined_lsi_manual",
  name = "umap_manual",
  minDist = 0.8,
  force = TRUE
)

# Add Clusters
proj <- addClusters(
  proj,
  reducedDims = "combined_lsi_manual",
  name = "manual_clusters",
  resolution = 0.75,
  force = TRUE
)

### Dimplot
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

# Load Data & Create Obj
gex_data <- Read10X_h5("data/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5")
seurat_project <- CreateSeuratObject(counts = gex_data$`Gene Expression`)

# Infer what cells were quality controlled by ArchR
qc_cells <- proj@cellColData@rownames
qc_cells_seurat <- substring(qc_cells, 10)
seurat_cell_ids <- colnames(gex_data$`Gene Expression`)
common_cell_ids <- intersect(seurat_cell_ids, qc_cells_seurat)

# Slice subset & Rename Cells to follow same naming standard
gex_data_subset <- subset(seurat_project, cells = common_cell_ids)

gex_data_subset <- RenameCells(gex_data_subset, new.names = paste0("PBMC_10k#", colnames(gex_data_subset)))

gex_data_subset <- NormalizeData(gex_data_subset, normalization.method = "LogNormalize", scale.factor = 10000)
gex_data_subset <- ScaleData(gex_data_subset, features = known_markers)

rm(gex_data)


metadata <- getCellColData(
  ArchRProj = proj
)

metadata$Barcodes <- rownames(metadata)
sorted_metadata <- metadata[match(colnames(gex_data_subset), metadata$Barcodes),]

sorted_metadata[["manual_clusters"]] <- as.numeric(gsub("C", "", sorted_metadata[["manual_clusters"]]))

gex_data_subset <- AddMetaData(gex_data_subset, sorted_metadata$manual_clusters, "manual_clusters")


# Get the UMAP coordinates from the ArchR project
umap_coords <- getEmbedding(ArchRProj = proj,
                            embedding = "umap_manual",
                            returnDF = TRUE)

#### THIS IS THE PROBLEM, NOTHING MATCHES
umap_coords <- umap_coords[match(colnames(gex_data_subset), rownames(umap_coords)),]

# Convert to a matrix
mat_umap_coords <- as.matrix(umap_coords)

# Create a Seurat DimReduc object
umap <- CreateDimReducObject(embeddings = mat_umap_coords)

# Add the DimReduc object to the Seurat object under the appropriate name
gex_data_subset[["umap_manual"]] <- umap

gex_data_subset <- SetIdent(gex_data_subset, value = "manual_clusters")
DimPlot(gex_data_subset, reduction = "umap_manual")
FeaturePlot(gex_data_subset, feature = "MS4A1")
