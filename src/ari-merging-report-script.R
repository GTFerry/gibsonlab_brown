library(ArchR)
library(Seurat)
library(Dune)
library(pheatmap)

results_dir <- "results/"

proj <- loadArchRProject(path = "saves/clustered")

print("Starting Dune Merging")

merger <- read.csv(paste(results_dir, "ClusterAssignments.csv", sep = ""), row.names=1)
cluster_names <- colnames(read.csv(paste(results_dir, "ClusterAssignments.csv", sep = "")))
cluster_names <- cluster_names[-1]

print("Read File")
print(cluster_names)

### ALL Plotted Together ###
# plot ARIs as barplot
getARIs <- ARIs(merger)

print("Got ARIs")
sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

print("Starting Plots")
# save plots
pdf(file = "RandIndex_heatmap_report.pdf", height = 15, width = 15)
plot.new()
plotARIs(clusMat = merger) + RotatedAxis() + title("Plot Before ARI Merging (Dune)")

print("Merging")
merger <- Dune(clusMat = merger, verbose = FALSE)
ari_matrix <- ARIs(merger$currentMat)
save(ari_matrix, file="saves/merged_data.Rdata")
print("Finished Merging")
title("Plot After ARI Merging (Dune)")
plot.new()
plotARIs(clusMat = merger$currentMat) + RotatedAxis()

plot.new()
title("Cluster Change in each clustering Method")
plotPrePost(merger)

print(cluster_names)

plot.new()

for (i in 1:(length(cluster_names))) {
    print(paste("Plotting Embedding:", paste("UMAP_Combined_", i, sep = "")))
#    print("Plotting Confusion Matrix")
#    cM <- confusionMatrix(paste0(proj$cluster_names[i]), paste0(proj$Sample))
#    cM <- cM / Matrix::rowSums(cM)
#    confmatrix <- pheatmap::pheatmap(
#        mat = as.matrix(cM),
#        color = paletteContinuous("whiteBlue"),
#        border_color = "black"
#    )
#
#    print(confmatrix)
#
#    print("Plotting UMAP")
    print(plotEmbedding(proj, name = cluster_names[i], embedding = paste("UMAP_Combined_", i, sep = "")))
}

dev.off()

