library(rasqualTools)
library(ArchR)
library(icesTAF)
#cell=c("Rod","Cone","ONBC","OFFBC","HC","AC","MG","RGC")
cell=c("Rod","Cone","ONBC","OFFBC","HC","AC","MG")
#region=c("mac","merged")
region=c("merged")
#region=c("lobe","mac","merged")
proj10=loadArchRProject(path = "Save-Proj20ppl64", force = FALSE, showLogo = TRUE)
peakset=getPeakSet(proj10)
peakset_new=data.frame(peakset)
peakset_new$id=paste0(peakset_new$seqnames,":", peakset_new$start,"-",   peakset_new$end)
gene_metadata=data.frame(gene_id=peakset_new$id,percentage_gc_content=peakset_new$GC)
gene_size = data.frame(gene_id=peakset_new$id, len=peakset_new$end-peakset_new$start+1)
#write.table(gene_metadata,file="/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/peak_count/peakset_new")
#write.table(gene_metadata,file="/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/peak_count/peakset_new",quote=F,sep="\t",row.name=F)
write.table(gene_metadata,file="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/peak_count/peakset_new",quote=F,sep="\t",row.name=F)
#setwd("/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/peak_count/")
setwd("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/peak_count/")
randomize <- function(x,g=NULL){
  if(is.null(g)){
    n=ncol(x);
    t(apply(x,1,function(xx){xx[order(runif(n))]}))
  }else{
    for(i in unique(g)){
      x[,g==i]=randomize(x[,g==i,drop=F])
    }
    x
  }
}

rasqualMakeCovariates <- function(counts, size_factors) {
### rasqualMakeCovariates <- function(counts, size_factors, out) {

  #Map parameters to Natsuhiko's variables
  Y=counts
  K=size_factors
  n=ncol(Y)

  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9

  # Singular value decomposition
  fpkm.svd   = svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
  fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))

  # Covariate selection
  sf=log(apply(Y,2,sum))
  covs=fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
#  if(cor(sf,covs[,1])^2<0.9){covs=cbind(sf, covs)}
  if(cor(sf,covs)^2<0.9){covs=cbind(sf, covs)}

  # Write covariates
  return(covs)
####  write.table(covs,file=out,quote=F,sep="\t")
}



for(reg in region){
for(c in cell){
#counts_matrix=read.table(paste0(c,"_snATAC_peak_count"),header=T)
counts_matrix=read.table(paste0(c,"_",reg,"_snATAC_peak_count"),header=T)

mkdir(paste0("rasqualTool_new/",c,"_",reg))
saveRasqualMatrices(list(cellType = counts_matrix), paste0("rasqualTool_new/",c,"_",reg), file_suffix = "expression")
size_factors = rasqualCalculateSampleOffsets(counts_matrix, gene_metadata, gc_correct = TRUE)
saveRasqualMatrices(list(cellType = size_factors), paste0("rasqualTool_new/",c,"_",reg), file_suffix = "size_factors_gc")
m=match(rownames(counts_matrix), gene_size$gene_id)
#cov=rasqualMakeCovariates(counts_matrix, gene_size$len[m] )
cov=rasqualMakeCovariates(counts_matrix,size_factors)
saveRasqualMatrices(list(cellType = cov), paste0("rasqualTool_new/",c,"_",reg), file_suffix = "covariate")
}
}
#####gene_data=data.frame(gene_id=peakset_new$id,chr=peakset_new$seqnames,strand=peakset_new$strand,exon_starts=as.character(peakset_new$start),exon_ends=as.character(peakset_new$end))
#####snp_id=read.table("/storage/chen/home/jw29/sc_human_retina/data/ca_eQTL/1000GP_Phase3_20ppl_all.phased.wRef.merge.vcf.snp_id")
#####colnames(snp_id)=c("chr","pos","snp_id")
#####snp_counts = countSnpsOverlapingExons(gene_data, snp_coords, cis_window = 200)
#####dplyr::select(snp_counts, gene_id, feature_snp_count, cis_snp_count)
######saveRasqualMatrices(list(cellType = snp_counts), paste0("rasqualTool/"), file_suffix = "snp_counts")



