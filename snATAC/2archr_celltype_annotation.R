library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
library(Seurat)
#proj1=loadArchRProject(path = "ArchRSubset_20ppl_1", force = FALSE, showLogo = TRUE)

addArchRGenome("hg19")


#saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj20ppl62", load = FALSE)
proj1=loadArchRProject(path = "Save-Proj20ppl62", force = FALSE, showLogo = TRUE)

seRNA0=readRDS("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/macular20ppl_merged_afterScPred_count.rds")
seRNA0@active.assay="RNA"

seRNA=subset(x=seRNA0, subset=scpred_max >=0.9)
seRNA@active.assay="RNA"


proj2 <- addGeneIntegrationMatrix(
    ArchRProj = proj1,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI4",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "scpred_prediction",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)


saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj20ppl62", load = FALSE)

cM <- as.matrix(confusionMatrix(proj2$Clusters20, proj2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

cM

unique(unique(proj2$predictedGroup_Un))

cPR <- paste0(c("Rod","Cone"), collapse="|")
cIN <- paste0(c("OFFBC", "ONBC", "RGC", "HC", "AC"), collapse="|")

cNN <- paste0(c("Mic", "MG", "Astro","Endo"), collapse="|")

clustPR <- rownames(cM)[grep(cPR, preClust)]
clustIN <- rownames(cM)[grep(cIN, preClust)]

clustNN <- rownames(cM)[grep(cNN, preClust)]

rnaPR <- colnames(seRNA)[grep(cPR, seRNA@meta.data$scpred_prediction)]
rnaIN <- colnames(seRNA)[grep(cIN, seRNA@meta.data$scpred_prediction)]

rnaNN <- colnames(seRNA)[grep(cNN, seRNA@meta.data$scpred_prediction)]

groupList <- SimpleList(
    PR = SimpleList(
        ATAC = proj2$cellNames[proj2$Clusters20 %in% clustPR],
        RNA = rnaPR
    ),
    IN = SimpleList(
        ATAC = proj2$cellNames[proj2$Clusters20 %in% clustIN],
        RNA = rnaIN
    ),
    NN = SimpleList(
        ATAC = proj2$cellNames[proj2$Clusters20 %in% clustNN],
        RNA = rnaNN
    )

)

proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI4",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupList = groupList,
    groupRNA = "scpred_prediction",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)


pal <- paletteDiscrete(values = seRNA@meta.data$scpred_prediction)

p1 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup_Un",
    pal = pal
)

p2 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup_Co",
    pal = pal
)

plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj20ppl62", load = FALSE)

proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj3,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI4",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "scpred_prediction",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

proj3 <- addImputeWeights(proj3)

saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj20ppl63", load = FALSE)
###########

cells=rownames(proj3[!(proj3@cellColData$Clusters20%in%c("C18","C15")),])

proj4=subsetArchRProject(  
  ArchRProj = proj3,
  cells = cells,
  outputDirectory = "Save-Proj20ppl63_rm_C15_C18",
  dropCells = TRUE)

cM <- confusionMatrix(proj4$Clusters20, proj4$predictedGroup)


labelOld <- rownames(cM)
labelOld



labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew



proj4$Clusters40 <- mapLabels(proj4$Clusters20, newLabels = labelNew, oldLabels = labelOld)

proj4@cellColData[proj4$Clusters20=="C16"&proj4$predictedGroup=="OFFBC",]$Clusters40 = "OFFBC"

proj4@cellColData[proj4$Clusters20=="C13"&proj4$predictedGroup=="ONBC",]$Clusters40 = "ONBC"

cM <- confusionMatrix(proj4$Clusters40, proj4$predictedGroup)




p1 <- plotEmbedding(proj4, colorBy = "cellColData", name = "Clusters40",embedding = "UMAP4")

p2 <- plotEmbedding(proj3, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP4")
p3 <- plotEmbedding(proj4, colorBy = "cellColData", name = "Sample", embedding = "UMAP4")


plotPDF(p1,p2,p3, name = "Plot-UMAP-Remap-Clusters20_cluster40.pdf", ArchRProj = proj4, addDOC = FALSE, width = 5, height = 5)



p4 <- plotEmbedding(proj4[proj4$predictedScore>=0.8,], colorBy = "cellColData", name = "Clusters40",embedding = "UMAP4")
p5 <- plotEmbedding(proj4[proj4$predictedScore>=0.5,], colorBy = "cellColData", name = "Clusters40",embedding = "UMAP4")


plotPDF(p4,p5, name = "Plot-UMAP-Remap-Clusters20_cluster40_ps_0.8_0.5.pdf", ArchRProj = proj4, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj20ppl63", load = FALSE)
saveArchRProject(ArchRProj = proj4, outputDirectory = "Save-Proj20ppl63_rm_C15_C18", load = FALSE)

write.table(proj4@cellColData,file="sc_human_retina/data/ASE/ATAC_cell_archr_04062021",sep="\t",quote=F)

cM1=confusionMatrix(proj4$Clusters40, proj4$Sample)
data.frame(cM1)
write.table(data.frame(cM1),file="sc_human_retina/data/ASE/ATAC_cell_archr_04062021_cell_count_sample",sep="\t",quote=F)

proj4$Sample_main=proj4$Sample
proj4$Sample_main=gsub("_Lobe|_Macular", "",proj4$Sample_main)
cM2=confusionMatrix(proj4$Clusters40, proj4$Sample_main)
data.frame(cM2)
write.table(data.frame(cM2),file="sc_human_retina/data/ASE/ATAC_cell_archr_04062021_cell_count_donor",sep="\t",quote=F)

cM3=confusionMatrix(proj4[proj4$predictedScore>=0.5,]$Clusters40, proj4[proj4$predictedScore>=0.5,]$Sample)
data.frame(cM3)
write.table(data.frame(cM3),file="sc_human_retina/data/ASE/ATAC_cell_archr_04062021_cell_count_sample_predScore_more_eq_0.5",sep="\t",quote=F)

p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "Clusters40", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)

plotPDF(plotList = p,
    name = "Plot-Tracks-Marker-Genes_cluster40_rm_C15_C18.pdf",
    ArchRProj = proj3,
    addDOC = FALSE, width = 5, height = 5)



##########

