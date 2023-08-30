library(ArchR)
library(Seurat)
library(Dune)

results_dir <- "results/"

proj <- loadArchRProject(path = "saves/clustered")

print("Starting Dune Merging")

merger <- read.csv(paste(results_dir, "ClusterAssignments.csv", sep = ""), row.names=1)
cluster_names <- colnames(read.csv(paste(results_dir, "ClusterAssignments.csv", sep = "")))
cluster_names <- cluster_names[-1]

print("Read File")

### ALL Plotted Together ###
# plot ARIs as barplot
getARIs <- ARIs(merger)

print("Got ARIs")
sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_report.pdf", height = 15, width = 15)
plot.new()
plot.new()
title("Plot Before ARI Merging (Dune)")
plotARIs(clusMat = merger) + RotatedAxis()

merger <- Dune(clusMat = merger, verbose = FALSE)
title("Plot After ARI Merging (Dune)")
plotARIs(clusMat = merger$currentMat) + RotatedAxis()

title("Cluster Change in each clustering Method")
plotPrePost(merger)

print(cluster_names)

plot.new()

for (i in 1:(length(cluster_names))) {
    print(paste("Plotting Embedding:", paste("UMAP_Combined_", i, sep = "")))
    print(plotEmbedding(proj, name = cluster_names[i], embedding = paste("UMAP_Combined_", i, sep = "")))
}

dev.off()
