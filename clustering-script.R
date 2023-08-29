library(ArchR)
library(Seurat)
library(Dune)

proj <- loadArchRProject(path = "saves/atac-rna-quality-controlled")

results_dir <- "results/"

elapsedTime <- function(difference) {
  total_secs <- as.numeric(difference, units = "secs")
  hours <- as.integer(total_secs %/% 3600)
  minutes <- as.integer((total_secs %% 3600) %/% 60)
  seconds <- as.integer(total_secs %% 60)
  
  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

appendToCSV <- function(newColumnDF, filename) {
  # Check if the input is a one-dimensional dataframe
  if(ncol(newColumnDF) != 1) {
    stop("Input dataframe must be one-dimensional!")
  }
  
  # Check if the file exists
  if(!file.exists(filename)) {
    stop(paste("File", filename, "does not exist!"))
  }
  
  # Read the existing CSV file
  existingDF <- read.csv(filename, stringsAsFactors = FALSE)
  
  # Check if the number of rows matches
  if(nrow(existingDF) != nrow(newColumnDF)) {
    stop("Number of rows in the new column does not match the existing CSV!")
  }
  
  # Add the new column to the existing dataframe
  colName <- colnames(newColumnDF)
  existingDF[[colName]] <- newColumnDF[[colName]]
  
  # Write the updated dataframe back to the CSV file
  write.csv(existingDF, file = filename, row.names = FALSE)
}

testClusters <- function(proj, resolutions_ATAC, resolutions_RNA,
                         iterationsList, dimsList, verbosity = 1,
                         writeLogs = TRUE,
                         logFile = "automatedClusteringLogfile.log") {
  
  if(verbosity >= 1) {
    start_time <- Sys.time()
    cat(
      strftime(start_time, format="%Y-%m-%d %H:%M:%S"),
      "- Starting Clustering Loop."
    )
  }
  
  # Check if all input lists are of the same length
  lengths <- c(length(resolutions_ATAC), length(resolutions_RNA), length(iterationsList), length(dimsList))
  
  if(length(unique(lengths)) != 1){
    stop("All input lists must be of the same length!")
  }
  
  first_col <- data.frame(1:nCells(proj))
  colnames(first_col) <- "Cell Index"
  
  # Write the dataframe to a CSV file
  write.csv(first_col, file = paste(results_dir, "ClusterAssignments.csv", sep = ""), row.names = FALSE)
  
  plots <- list()
  clusterAssignments <- data.frame(cellBarcode = getCellNames(proj))
  
  for(i in 1:length(resolutions_ATAC)){
    
    if(verbosity >= 1){
      cat(
        strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"),
        "- Starting iteration", i, 
        "of", length(resolutions_ATAC), "\n")
    }
    
    res_ATAC <- resolutions_ATAC[i]
    res_RNA <- resolutions_RNA[i]
    iters <- iterationsList[i]
    dims <- dimsList[i]
    
    ATAC_LSI_name <- paste0("LSI_ATAC_", i, "_Res_", gsub("\\.", "-", as.character(res_ATAC)), "_Iters_", iters, "_Dims_1to", dims)
    
    if(verbosity >= 2){
      cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"),
          "- Performing dimensionality reduction and clustering for ATAC with \n",
          "Resolution =", res_ATAC, "\n",
          "Iterations = ", iters, "\n",
          "Dimensions = 1 :", dims, "\n",
          "Name = ", ATAC_LSI_name, "\n")
    }
    
    if (!writeLogs) {
      sink(logFile, append=TRUE)
    }
    
    # Dimensionality reduction and clustering for ATAC
    proj <- addIterativeLSI(
      ArchRProj = proj,
      iterations = iters,
      clusterParams = list(resolution = res_ATAC, sampleCells = 10000, n.start = 10),
      saveIterations = FALSE,
      useMatrix = "TileMatrix",
      depthCol = "nFrags",
      dimsToUse = 1:dims,
      name = ATAC_LSI_name,
      force = TRUE
    )
    
    if (!writeLogs) {
      sink()
    }
    
    RNA_LSI_name <- paste0("LSI_RNA_", i, "_Res_", gsub("\\.", "-", as.character(res_RNA)), "_Iters_", iters, "_Dims_1to", dims)
    
    if(verbosity >= 2){
      cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"),
          "- Performing dimensionality reduction and clustering for RNA with:",
          "Resolution =", res_RNA,
          "Iterations = ", iters, "\n",
          "Dimensions = 1 :", dims, "\n",
          "Name = ", RNA_LSI_name, "\n")
    }
    
    if (!writeLogs) {
      sink(logFile, append=TRUE)
    }
    
    # Dimensionality reduction and clustering for RNA
    proj <- addIterativeLSI(
      ArchRProj = proj,
      iterations = iters,
      clusterParams = list(resolution = res_RNA, sampleCells = 10000, n.start = 10),
      saveIterations = FALSE,
      useMatrix = "GeneExpressionMatrix", 
      depthCol = "Gex_nUMI",
      varFeatures = 2500,
      firstSelection = "variable",
      binarize = FALSE,
      dimsToUse = 1:dims,
      name = RNA_LSI_name,
      force = TRUE
    )
    
    if (!writeLogs) {
      sink()
    }
    
    # Combined dimensions
    if(verbosity >= 2){
      cat(
        strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"),
        "- Combining dimensions for ATAC and RNA\n"
      )
    }
    
    COMB_LSI_name <- paste0("LSI_Combined_", i, "_ResRNA_", gsub("\\.", "-", as.character(res_RNA)), "_ResATAC_", gsub("\\.", "-", as.character(res_ATAC)), "_Iters_", iters, "_Dims_1to", dims)
    
    if (!writeLogs) {
      sink(logFile, append=TRUE)
    }
    
    proj <- addCombinedDims(
      proj, 
      reducedDims = c(
        ATAC_LSI_name,
        RNA_LSI_name
      ), 
      name = COMB_LSI_name
    )
    
    if (!writeLogs) {
      sink()
    }
    
    # UMAP and clustering for each dimensionality reduction
    if(verbosity >= 2){
      cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"), "- Applying UMAP and clustering \n")
    }
    
    if (!writeLogs) {
      sink(logFile, append=TRUE)
    }

    proj <- addUMAP(
      proj,
      reducedDims = ATAC_LSI_name,
      name = paste0("UMAP_ATAC_", i),
      minDist = 0.8,
      force = TRUE
    )

    proj <- addUMAP(
      proj,
      reducedDims = RNA_LSI_name,
      name = paste0("UMAP_RNA_", i),
      minDist = 0.8,
      force = TRUE
    )

    proj <- addUMAP(
      proj,
      reducedDims = COMB_LSI_name,
      name = paste0("UMAP_Combined_", i),
      minDist = 0.8,
      force = TRUE
    )

    if (!writeLogs) {
      sink()
    }

    cluster_name <- paste0("Cluster_", i, "_ResRNA_", gsub("\\.", "-", as.character(res_RNA)), "_ResATAC_", gsub("\\.", "-", as.character(res_ATAC)), "_Iters_", iters, "_Dims_1to", dims)

    if(verbosity >= 2){
      cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"), "- Saving clusters with name:", cluster_name, "\n")
    }

    if (!writeLogs) {
      sink(logFile, append = TRUE)
    }

    # TODO: check assumption of resolution here
    proj <- addClusters(
      proj,
      reducedDims = COMB_LSI_name,
      name = cluster_name,
      resolution = (res_ATAC + res_RNA)/2,
      force = TRUE
    )

    if (!writeLogs) {
      sink()
    }

    if(verbosity >= 1){
      cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"), "- Saving cluster assignments to CSV\n")
    }

    currentClusteringLabels <- data.frame(proj@cellColData@listData[[cluster_name]])

    colnames(currentClusteringLabels) <- c(cluster_name)

    # View(currentClusteringLabels)

    appendToCSV(currentClusteringLabels, paste(results_dir, "ClusterAssignments.csv", sep = ""))

    if(verbosity >= 1){
      cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"), "- Finished iteration", i, "\n")
    }
  }

  if(verbosity >= 1) {
    difference <- Sys.time() - start_time
    cat(strftime(Sys.time(), format="%Y-%m-%d %H:%M:%S"), "- Finished Loops, total time: ", elapsedTime(difference), "\n")
  }

  return(proj)
}

input_data <- read.csv("input.csv")

resolutions_ATAC_csv <- input_data$resolutions_ATAC
resolutions_RNA_csv  <- input_data$resolutions_RNA
iterationsList_csv   <- input_data$iterationsList
dimsList_csv         <- input_data$dimsList

proj <- testClusters(
  proj,
  resolutions_ATAC = resolutions_ATAC_csv,
  resolutions_RNA  = resolutions_RNA_csv,
  iterationsList   = iterationsList_csv,
  dimsList         = dimsList_csv,
  verbosity        = 2
)

saveArchRProject(ArchRProj = proj, outputDirectory = "saves/clustered", load = FALSE)

print("Starting Dune Merging")

merger <- read.csv(paste(results_dir, "ClusterAssignments.csv", sep = ""), row.names=1)

print(head(merger))

print("Read File")

### ALL Plotted Together ###
# plot ARIs as barplot
getARIs <- ARIs(merger)

print("Got ARIs")
sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_report.pdf", height = 15, width = 15)
plotARIs(clusMat = merger) + RotatedAxis()
dev.off()

# print("Plotting Before ARI")
# print("Merging")
merger <- Dune(clusMat = merger, verbose = FALSE)
# print("Merged")
# pdf(file = "RandIndex_heatmap_all_after.pdf", height = 15, width = 15)
plotARIs(clusMat = merger$currentMat) + RotatedAxis()
# dev.off()

# print("Plotting After ARI")

# pdf(file = "RandIndex_before_after_clusters.pdf", height = 15, width = 15)
plotPrePost(merger)
dev.off()
