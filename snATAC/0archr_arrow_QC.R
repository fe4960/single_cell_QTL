library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 

##########inputFiles=read.table("sc_human_retina/scripts/single_cell/macLobe_ATAC20ppl_list",row.names=1)

inputFiles=read.table("sc_human_retina/scripts/single_cell/macular_ATAC20ppl_list",row.names=1)
inputFiles1=as.vector(inputFiles$V2)
names(inputFiles1)=row.names(inputFiles)
inputFiles=inputFiles1
inputFiles=as.list(inputFiles)
addArchRGenome("hg19")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)



doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Save-Proj20ppl1",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

paste0("Memory Size = ", round(object.size(proj1) / 10^6, 3), " MB")

#getAvailableMatrices(projHeme1)

getAvailableMatrices(proj1)
proj1 <- filterDoublets(ArchRProj = proj1)

p1 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

p2 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p3 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

p4 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj1, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj1)
p2 <- plotTSSEnrichment(ArchRProj = proj1)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)



saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj20ppl1", load = FALSE)


