library(SoupX)
library(Seurat)
library(cowplot)
library(dplyr)
#devtools::install_github('chris-mcginnis-ucsf/DoubletFinder',force=TRUE)
library(DoubletFinder)
library(DropletUtils)


file=read.table("sc_human_retina/scripts/single_cell/macular_10x3_7_0")
#sam=c("Macular_19D13","Macular_19D14","Macular_19D15","Macular_19D16")
#ifnb.list=list()
for (i in 1:length(file$V1)){
#cur_dir=paste("/storage/novaseq/Data/191016_A00431_0119_AHLCT2DSXX/Data/Intensities/BaseCalls/10x/10x_",sam[i],"_Nu/outs/filtered_feature_bc_matrix",sep="")
cur_dir=paste0("/storage/chen/data_share_folder/human20_10X3/backup/",file$V1[i],"/outs/")
#sc = load10X("/storage/chen/data_share_folder/human20_10X3/backup/10x3_19_D005/outs/")
print(file$V1[i])

sc = load10X(cur_dir)
sc = autoEstCont(sc)
out = adjustCounts(sc)

new_dir=paste0("sc_human_retina/data/single_cell/snRNA/afterSoupX/",file$V1[i])
write10xCounts(x=out, version="3", path=new_dir)

gg = plotMarkerMap(sc, "NRL")
pdf(paste0("sc_human_retina/data/single_cell/snRNA/afterSoupX/",file$V1[i],"_MarkerMap_NRL.pdf"))
print(plot(gg))
dev.off()

pdf(paste0("sc_human_retina/data/single_cell/snRNA/afterSoupX/",file$V1[i],"_SoupX_NRL.pdf"))
print(plotChangeMap(sc, out, "NRL"))
dev.off()

data1=CreateSeuratObject(counts = sc$toc,  min.cells = 5)
data2=CreateSeuratObject(counts = out,  min.cells = 5)
data1[["percent.mt"]] <- PercentageFeatureSet(data1, pattern = "^MT-")
data2[["percent.mt"]] <- PercentageFeatureSet(data2, pattern = "^MT-")

pdf(paste0("sc_human_retina/data/single_cell/snRNA/afterSoupX/",file$V1[i],"_QC1.pdf"),width=20)
p1=VlnPlot(data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2=VlnPlot(data2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(list(p1,p2),ncol=2)
dev.off()
#pdf()
#pdf("QC2_Scatter_19_D005.pdf",width=15)
pdf(paste0("sc_human_retina/data/single_cell/snRNA/afterSoupX/",file$V1[i],"_QC2.pdf"),width=20)

p1 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p4 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(list(p1,p2,p3,p4),ncol=4)
dev.off()


}


