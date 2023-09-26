library(ArchR)

archR_project <- loadArchRProject("saves/labeled")

print("Loaded Project. Starting Pseudo Bulk...")

archR_project <- addGroupCoverages(ArchRProj = archR_project, groupBy = "SeuratCluster")

print("Pseudo Bulk completed. Saving...")

saveArchRProject(archR_project, "saves/pseudobulked", load = FALSE, overwrite = TRUE)

print("Starting Peak calling...")

pathToMacs2 <- findMacs2()

archR_project <- addReproduciblePeakSet(
  ArchRProj = archR_project, 
  groupBy = "SeuratCluster", 
  pathToMacs2 = pathToMacs2
)

getPeakSet(archR_project)

print("Peak Calling Finished. Saving...")

saveArchRProject(archR_project, "saves/peak-called", load = FALSE, overwrite = TRUE)