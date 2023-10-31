library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 

inputFiles=read.table("sc_human_retina/scripts/single_cell/macLobe_ATAC20ppl_list_2",row.names=1)

#inputFiles=read.table("sc_human_retina/scripts/single_cell/macular_ATAC20ppl_list",row.names=1)
inputFiles1=as.vector(inputFiles$V2)
names(inputFiles1)=row.names(inputFiles)
inputFiles=inputFiles1
addArchRGenome("hg19")
#Dont set this too high because you can always increase later

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
#ArrowFiles=read.table("sc_human_retina/data/single_cell/arrowFile_list")
#ArrowFiles=ArrowFiles$V1
#addArchRGenome("hg19")

#ArrowFiles=read.table("sc_human_retina/data/single_cell/snATAC/20ppl_mac_lobe_list")
#ArrowFiles=ArrowFiles$V1



