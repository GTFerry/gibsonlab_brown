library(ArchR)
library(Seurat)
library(Dune)

results_dir <- "results/"

proj <- loadArchRProject(ArchRProj = proj)

print("Starting Dune Merging")

merger <- read.csv(paste(results_dir, "ClusterAssignments.csv", sep = ""), row.names=1)
cluster_names <- colnames(read.csv(paste(results_dir, "ClusterAssignments.csv", setp = "")))

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
title("Plot Before ARI Merging (Dune)")

merger <- Dune(clusMat = merger, verbose = FALSE)
plotARIs(clusMat = merger$currentMat) + RotatedAxis()
title("Plot After ARI Merging (Dune)")

plotPrePost(merger)
title("Cluster Change in each clustering Method")
dev.off()

for (i in 1:length(cluster_name)) {
    plotEmbedding(proj, name = cluster_names[i], embedding = paste("UMAP_Combined_", i, sep = ""), size = 1.5, labelAsFactors=F, labelMeans=F)
    title(paste("UMAP Visualization for: ", cluster_names[i], sep = ""))
}

