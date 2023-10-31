library(EWCE)
library(MAGMA.Celltyping)
library(ggplot2)
library(dplyr)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)


#rds="/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_cleansample_downsample5000_group.rds"
#obs=read.table("/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_cleansample_downsample5000_group.obs.gz",header=T,comment.char="")
dirIn="/storage/chenlab/Users/junwang/human_meta/data/GWAS/MAGMA"
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

rds="/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_v1_downsample5000.rds"
obs=read.table("/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_v1_downsample5000.obs.gz",header=T,comment.char="",sep="\t")
exp=readRDS(rds)

sce=as.SingleCellExperiment(exp)

#scpred=readRDS("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/macular20ppl_merged_afterScPred_umap_flt_scpred_2fpkm.rds")
#annotLevels <- list(l1=obs$celltype) #list(l1 = l1, l2 = l2)
annotLevels <- list(l1=obs$majorclass) #list(l1 = l1, l2 = l2)

fNames_ALLCELLS <- EWCE::generate_celltype_data(
    exp=sce,
    annotLevels = annotLevels,
    groupName = "retina"
)

retina=load(file=fNames_ALLCELLS)
#file=read.table("/storage/chen/home/jw29/sc_human_retina/scripts/GWAS_new/3MAGMA_list_all")
#file=read.table("/storage/chen/home/jw29/file_t")
file=read.table("/storage/chen/home/jw29/human_meta/scripts/GWAS/MAGMA_list_100")
#file=read.table(paste0("human_meta/scripts/GWAS/",args[1]))
for(i in 1:length(file$V1)){
MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
# magma_dirs= "/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Springelkamp2017_28073927/CA/MAGMA_Files/Meta_CA_caucasians_reform_MungeSumstats.txt.35UP.10DOWN/", 
magma_dirs = file$V1[i],
# magma_dirs="/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Chen_32888493/Chen_32888493cell/MAGMA_Files/GCST90002374_WBC_MungeSumstats_rmFRQ.txt_reform.35UP.10DOWN/",
 ctd = ctd,
  ctd_species = "human", 
  ctd_name = "retina", 

  run_linear = TRUE, 
  run_top10 = TRUE)

merged_results <- MAGMA.Celltyping::merge_results(
MAGMA_results = MAGMA_results)

##data=as.data.frame(MAGMA_results)

#if(dim(data)[1]>0){
#write.table(merged_results,file=paste0("MAGMA_result_atlas"),quote=F,sep="\t")

write.table(merged_results,file=paste0(file$V1[i],"MAGMA_result_atlas"),quote=F,sep="\t")
#}
}


