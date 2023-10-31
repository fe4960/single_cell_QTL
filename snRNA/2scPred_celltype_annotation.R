library(scPred)
library(Seurat)
library(magrittr)
set.seed(12345)

#reference=readRDS("sc_human_retina/data/single_cell/scPred/Human_retina_combined_reference.rds")
reference=readRDS("sc_human_retina/data/single_cell/scPred/Human_retina_combined_reference_new.rds")

#query=readRDS("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/macular20ppl.rds")
sample_list=read.table("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/sample_list_5_0")

#sample_list<-c("D028_13", "D027_13", "D026_13", "D021_13", "D019_13") #, "D018_13", "D017_13", "D013_13", "D009_13","D005_13",
#               "D030_13", "19_D019", "19_D011", "19_D010","19_D009", "19_D008",  "19_D007", "19_D006","19_D005", "19_D003")
#head<-"/storage/novaseq/Data/200318_A00431_0173_AHNGHTDSXX/Data/Intensities/BaseCalls/10x/10x3_"
#tail<-"/outs/filtered_feature_bc_matrix/"
#setwd("/storage/chen/data_share_folder/human_20_snRNAseq/output_1/")
#dir_save="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/macular20ppl_nFeature500_mt15/"
dir_save="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/scPred_new/"

in_dir="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/"
for (sample in sample_list$V1){
#dir=paste0("/storage/chen/data_share_folder/human20_10X3/post_QC_nFeature500_mt15/","10x3_",sample)
dir=paste0(in_dir,"/10x3_",sample)

#dir=paste0("/storage/chen/data_share_folder/human20_10X3/post_QC/","10x3_",sample)
#query <- Read10X(data.dir ="/storage/chen/data_share_folder/human20_10X3/post_QC/10x3_19_D003/")
query <- Read10X(data.dir =dir)
query <- CreateSeuratObject(counts = query)
######D003_19$region <- "macular"

query$sample <- sample

query <- NormalizeData(query)

query <- scPredict(query, reference)



#crossTab(query, "cell_type", "scpred_prediction")


file1=paste0(dir_save, sample, "_scPred.pdf")
#pdf("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/19_D003_scPred.pdf")
pdf(file1)
print(DimPlot(query, group.by = "scpred_prediction", reduction = "scpred"))
dev.off()



query <- SCTransform(query, verbose = FALSE)
query <- RunPCA(query, verbose = FALSE)
#query<-RunUMAP(query, reduction = "pca", dims = 1:15, umap.method='umap-learn',metric = 'correlation')
query<-RunUMAP(query, reduction = "pca", dims = 1:15)
query <- RunTSNE(query, dims = 1:15, verbose = FALSE)
query <- FindNeighbors(query, dims = 1:15, verbose = FALSE)
query <- FindClusters(query, verbose = FALSE, resolution = 0.5)

#query<-RunUMAP(query, reduction = "pca", dims = 1:15)

#pdf("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/19_D003_scPred_umap.pdf")
file1=paste0(dir_save, sample, "_umap.pdf")
pdf(file1)
print(DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE))
dev.off()


markerGenes <- c("ARR3","APOE","GLUL","PROX1","PCP2","PAX6","NRL","OTX2","ONECUT2","GFAP","SOD3","GAD1","GAD2","SLC6A9","TPBG","TBR1","VSX1")

file1=paste0(dir_save, sample, "_dotPlot.pdf")
pdf(file1,width=12)
#pdf("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/19_D003_scPred_dotPlot.pdf")
print(DotPlot(query, features = markerGenes, cols = "red", group.by = "scpred_prediction"))
dev.off()

file1=paste0(dir_save, sample, "_scPred.rds")
#saveRDS(query, "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/19_D003_scPred_out.rds")
saveRDS(query, file1)
mt=query@meta.data
table(mt$scpred_prediction)
table(mt[mt$scpred_max>=0.85,]$scpred_prediction)
  ms_add<-paste0(dir_save, sample, "_maxScore.pdf")
  pdf(ms_add, width = 8, height = 6)
  print(FeaturePlot(object = query, features = "scpred_max"))
  dev.off()
  gc_add<-paste0(dir_save, sample, "_genecount.pdf")
  pdf(gc_add, width = 8, height = 6)
  print(FeaturePlot(object = query, features = "nFeature_RNA", max.cutoff = "q99"))
  dev.off()
  Conf_test_scPred<- table(query@active.ident, query@meta.data$scpred_prediction)
  cm_add<-paste0(dir_save, sample, "_conf.csv")
  write.csv(Conf_test_scPred, cm_add)
  mt_add<-paste0(dir_save, sample, "_metadata.csv")
  write.csv(query@meta.data, mt_add)

#query.markers <- FindAllMarkers(query, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#query.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#mt_add<-paste0(dir_save, sample, "_marker.csv")
#  write.csv(query.markers, mt_add)
}

#pdf("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/19_D003_scPred_featurePlot.pdf")
#FeaturePlot(query, features =markerGenes)
#dev.off()
