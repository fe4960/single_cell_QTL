args<-commandArgs(TRUE)
library(Seurat)
library(cowplot)
library(dplyr)
#devtools::install_github('chris-mcginnis-ucsf/DoubletFinder',force=TRUE)
library(DoubletFinder)
library(DropletUtils)
#file=read.table("/storage/chen/data_share_folder/human20_10X3/file_list")
file=read.table("/storage/chen/home/jw29/sc_human_retina/scripts/single_cell/macular_10x3_7_0")
setwd("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/")
#sam=c("Macular_19D13","Macular_19D14","Macular_19D15","Macular_19D16")
#ifnb.list=list()
for (i in 1:length(file$V1)){
#cur_dir=paste("/storage/novaseq/Data/191016_A00431_0119_AHLCT2DSXX/Data/Intensities/BaseCalls/10x/10x_",sam[i],"_Nu/outs/filtered_feature_bc_matrix",sep="")
#cur_dir=paste("/storage/chen/data_share_folder/human20_10X3/backup/",file$V1[i],"/outs/filtered_feature_bc_matrix",sep="")
cur_dir=file$V1[i]
#cur_dir=file$V1[i]
data=Read10X(data.dir=cur_dir)
#ifnb.list[[file$V1[i]]]=CreateSeuratObject(count=ifnb.list[[file$V1[i]]], project=file$V1[i], min.cells = 10, min.features = 500)
data=CreateSeuratObject(count=data, min.cells = 5)
print(dim(data))
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-")
#ifnb.list[[file$V1[i]]] <- subset(ifnb.list[[file$V1[i]]], subset = percent.mt <= 15)
data <- subset(data, subset = nFeature_RNA >= 500 & percent.mt <= 15) 
print(dim(data))
#QC1=paste("/storage/chen/data_share_folder/human20_10X3/QC_nFeature500_mt15/mac_",file$V1[i],"_QC1_flt.pdf",sep="")
QC1=paste0("QC_nFeature500_mt15/mac_",file$V1[i],"_QC1_flt.pdf")

pdf(QC1,width=11)
print(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#QC2=paste("/storage/chen/data_share_folder/human20_10X3/QC_nFeature500_mt15/mac_",file$V1[i],"_QC2_flt.pdf",sep="")
QC2=paste0("QC_nFeature500_mt15/mac_",file$V1[i],"_QC2_flt.pdf")
pdf(QC2,width=11)
print(CombinePlots(plots = list(plot1, plot2)))
dev.off()
data <- SCTransform(data, verbose = FALSE)
data <- RunPCA(data, verbose = FALSE)

# t-SNE and Clustering
data <- RunUMAP(data, reduction = "pca", dims = 1:10)
########ifnb.list_new[[sam[i]]] <- NormalizeData(ifnb.list_new[[sam[i]]] )
#######ifnb.list_new[[sam[i]]] <- ScaleData(ifnb.list_new[[sam[i]]])
########ifnb.list_new[[sam[i]]] <- FindVariableFeatures(ifnb.list_new[[sam[i]]], selection.method = "vst", nfeatures = 2000)
########ifnb.list_new[[sam[i]]] <- RunPCA(ifnb.list_new[[sam[i]]])
#######ifnb.list_new[[sam[i]]] <- RunUMAP(ifnb.list_new[[sam[i]]], dims = 1:10)
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data, resolution = 0.5)

######sweep.res.list_mac <- paramSweep_v3(ifnb.list_new[[sam[i]]], PCs = 1:20, sct = FALSE)
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = TRUE)

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
#sweep.stats_mac <- summarizeSweep(x, GT = FALSE)

bcmvn_mac <- find.pK(sweep.stats)

pK=(bcmvn_mac[bcmvn_mac$BCmetric==max(bcmvn_mac$BCmetric),'pK'])[1]
pK=as.numeric(as.character(pK))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
cutoffDroublet <- file$V6[i]/1000 * 0.01 # assign doublets for 7.5% of total cells
annotations    <- data@meta.data$seurat_cluster
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi       <- round(cutoffDroublet*length(colnames(data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
#nExp_poi.adj 

data   <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
tmp_pANN = paste("pANN_0.25_",pK,"_",nExp_poi,sep="")
#ifnb.list[[file$V1[i]]]   <- doubletFinder_v3(ifnb.list[[file$V1[i]]], PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = tmp_pANN, sct = TRUE)
#tmp_class1=paste("DF.classifications","_0.25_",pK,"_",nExp_poi.adj,sep="")

tmp_class=paste("DF.classifications","_0.25_",pK,"_",nExp_poi,sep="")
propDoublet <- prop.table(table(data@meta.data[[tmp_class]]))
print(propDoublet)
umap_doublet=paste0("QC_nFeature500_mt15/mac_",file$V1[i],"_umap_doublet_sct_flt.pdf")

#umap_doublet=paste("/storage/chen/data_share_folder/human20_10X3/QC_nFeature500_mt15/mac_",file$V1[i],"_umap_doublet_sct_flt.pdf",sep="")
#pdf("5.1.2_umap_integrated.pdf")
pdf(umap_doublet,width=11)
print(DimPlot(object = data, reduction = "umap", label = T,label.size = 5))
# umap - w/o doublet
# need to replace the right col name for 'DF.classifications_0.25_0.09_723' in the meta-data
print(DimPlot(object = data, reduction = "umap", group.by = tmp_class,label = T,label.size = 5))
#print(DimPlot(object = data, reduction = "umap", group.by = tmp_class1,label = T,label.size = 5))

dev.off()


data@meta.data$dblt=data@meta.data[[tmp_class]]
data=subset(data, subset=dblt =="Singlet")
print(dim(data))
#data_output_dir_1="/storage/chen/data_share_folder/human20_10X3/post_QC_nFeature500_mt15/"
data_output_dir=paste0("QC_nFeature500_mt15/",file$V1[i])
write10xCounts(x=data@assays$RNA@counts, version="3", path=data_output_dir)
# umap
#setwd(dirOut)
#head(ifnb.list_new[[file$V1[i]]]@meta.data)

# need to replace col name for 'pANN_0.25_0.09_723' in the meta data
#QC3=paste("/storage/chen/data_share_folder/human20_10X3/QC_nFeature500_mt15/mac_",file$V1[i],"_QC_afterRm_doublet_sct_flt.pdf",sep="")
QC3=paste0("QC_nFeature500_mt15/mac_",file$V1[i],"_QC_afterRm_doublet_sct_flt.pdf")

print(pdf(QC3,width=11))
print(FeaturePlot(data, features = "percent.mt", pt.size = 0.1, ncol = 1))
print(FeaturePlot(data, features = "nCount_RNA", pt.size = 0.1, ncol = 1))
print(FeaturePlot(data, features = "nFeature_RNA", pt.size = 0.1, ncol = 1))
print(FeaturePlot(data, features = tmp_pANN, pt.size = 0.1, ncol = 1))
dev.off()

cell_single=rownames(data@meta.data)
#write(cell_single,file=paste0("/storage/chen/data_share_folder/human20_10X3/QC_nFeature500_mt15/mac_",file$V1[i],"_doublet_sct_flt"))
write(cell_single,file=paste0("QC_nFeature500_mt15/mac_",file$V1[i],"_doublet_sct_flt"))

}


