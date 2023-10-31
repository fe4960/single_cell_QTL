library(Seurat)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)

#scp<-readRDS("/storage/singlecell/qingnanl/human_10x_process/25_scPred/1_full_ref/scp.rds")
#sample_list<-c("D028_13", "D027_13", "D026_13", "D021_13", "D019_13", "D018_13", "D017_13", "D013_13", "D009_13","D005_13",
#               "D030_13", "19_D019", "19_D011", "19_D010","19_D009", "19_D008",  "19_D007", "19_D006","19_D005", "19_D003")
#sample_list=read.table("sc_human_retina/data/single_cell/snRNA/scPred/sample_list")
sample_list=read.table("sc_human_retina/data/single_cell/snRNA/scPred/sample_list")
#head<-"/storage/novaseq/Data/200318_A00431_0173_AHNGHTDSXX/Data/Intensities/BaseCalls/10x/10x3_"
#tail<-"/outs/filtered_feature_bc_matrix/"
##############!setwd("/storage/chen/data_share_folder/human20_10X3/post_QC_nFeature500_mt15/")
setwd("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/")

#setwd("/storage/chen/data_share_folder/human_20_snRNAseq/output_1/")
#celltype=c("Rod", "BC",    "AC", "MG", "Astro","HC","RGC", "Cone")
celltype=c("BC", "Astro",   "Rod",     "MG",      "Cone",    "RGC",     "HC", "AC") # "Mic")
for (sample in sample_list$V1){
#  ct_info=read.csv(paste0("10x3_",sample,"_metadata.csv"),header=T)
#  address<-paste0(head, sample, tail)
##########!ct_info=read.csv(paste0("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/scPred/macular20ppl_nFeature500_mt15/", sample, "_metadata.csv"), header=T)
ct_info=read.csv(paste0("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/scPred_new/",sample, "_metadata.csv"), header=T)
  address = paste0("10x3_", sample)
  query <- Read10X(data.dir = address)
#  query <- CreateSeuratObject(counts = query,  min.cells = 5, min.features = 500)
  query <- CreateSeuratObject(counts = query)

  DefaultAssay(query)<-"RNA"
  query.sce <- as.SingleCellExperiment(query)
  query_counts <- counts(query.sce)
#  query_cpm  <- apply(query_counts, 2, function(x) (x/sum(x))*1000000)
  dim_num = dim(query_counts)
  exp_matrix=matrix(NA, nrow=dim_num[1],ncol=length(celltype))
  colnames(exp_matrix) = celltype
  rownames(exp_matrix) = rownames(query_counts)
  print(sample)
  for(ct in celltype){
  print(ct)
  cols=ct_info[ct_info$scpred_max>=0.9&ct_info$scpred_prediction==ct,]$X

  #cols=ct_info[ct_info$scpred_max>=0.9&ct_info$scpred_prediction==ct&ct_info$sample==address,]$X
#  pattern=paste0(address,"_")
#  cols1=gsub(pattern,"",cols)
#  if(length(cols)>1){
  query_cpm_tmp = rowSums(query_counts[,cols])

#  query_cpm_tmp = rowSums(query_counts[,paste0("10x3_",cols1)])
  total=sum(query_cpm_tmp)
#  exp_matrix[,ct]=rowMeans(query_cpm[,cols])
  exp_matrix[,ct]=query_cpm_tmp/total*1000000
 #  }
  }

#write.table()
######write.table(exp_matrix,file=paste0("/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq/",sample,"_macula_fovea_snRNA_pool_ave_new"),quote=F,sep="\t")
write.table(exp_matrix,file=paste0("/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq/",sample,"_macula_fovea_snRNA_pool_ave_new_afterSoupX_new"),quote=F,sep="\t")

}
