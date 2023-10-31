library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
library(Seurat)
#proj1=loadArchRProject(path = "ArchRSubset_20ppl_1", force = FALSE, showLogo = TRUE)

addArchRGenome("hg19")

proj1=loadArchRProject(path = "Save-Proj20ppl11", force = FALSE, showLogo = TRUE)

proj1 <- addIterativeLSI(
    ArchRProj = proj1,
    useMatrix = "TileMatrix",
    name = "IterativeLSI4",
    iterations = 5,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1,0.2,0.4,0.6,0.8),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 30000,
    dimsToUse = 1:40,
    LSIMethod = 2,
    selectionMethod = "var",
    filterBias = TRUE
)

proj1 <- addHarmony(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI4",
    name = "Harmony4",
    groupBy = "Sample",
    corCutOff = 0.25,
    lambda = 0.75,
    sigma = 0.2
)

proj1 <- addClusters(input = proj1, reducedDims = "IterativeLSI4",method = "Seurat",    name = "Clusters20",    resolution = 1.5)


cM <- confusionMatrix(paste0(proj1$Clusters20), paste0(proj1$Sample))


proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "IterativeLSI4", name = "UMAP4",     nNeighbors = 80, minDist = 0.45, metric = "cosine")
p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP4")
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters20", embedding = "UMAP4")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters20.pdf",
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)



proj1 <- addUMAP(
    ArchRProj = proj1, 
    reducedDims = "Harmony4", 
    name = "UMAPHarmony4", 
    nNeighbors = 80, 
    minDist = 0.45, 
    metric = "cosine"
)


p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony4")
p4 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters20", embedding = "UMAPHarmony4")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters20.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj20ppl31", load = FALSE)

markerGenes <- c("GAD1","GAD2","GFAP","ARR3","ESAM","ONECUT2","APOE","RGR","CD74","GRIK1","GRM6","NEFM","PDE6G")
p <- plotEmbedding(
    ArchRProj = proj1,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
#    embedding = "TSNEHarmony", #    embedding = "UMAP",
    embedding = "UMAP4",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)




plotPDF(plotList = p,
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
    ArchRProj = proj1,
    addDOC = FALSE, width = 5, height = 5)

proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "IterativeLSI4", name = "UMAP5",     nNeighbors = 30, minDist = 0.5, metric = "cosine")
p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP5")
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters20", embedding = "UMAP5")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters21.pdf",
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

proj1 <- addClusters(input = proj1, reducedDims = "Harmony4",method = "Seurat",    name = "Clusters22",    resolution = 1.5)

