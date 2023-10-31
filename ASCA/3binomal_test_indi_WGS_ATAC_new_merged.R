#contig  position        variantID       refAllele       altAllele       refCount        altCount        totalCount      lowMAPQDepth    lowBaseQDepth   rawDepth       otherBases      improperPairs
#chr1    10583   chr1:10583      G       A       114     1       115     0       0       116     0       1

#/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam/out-cluster-ATAC_D026_13_BC.snp.count_higher0_binP_PASS_qval_WGS
#chr1    940390  chr1:875770     A       G       6       13      19      0       0       20      0       1       0.167068481445313       1       19      18
#sample=read.table("sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl_snp_list_30_0")
#sample=read.table("sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl_snp_list_30_0")

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sample=read.table(args[1])
for(i in 1:length(sample$V1)){

#file=paste0("/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-ATAC_",sample$V1[i],".filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio")
file=paste0("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/",sample$V1[i],".snATAC.merged.filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio")

if(file.size(file)>0){
M=read.table(file)

fisher_test=function(x){
TeaTasting = matrix(c(as.numeric(x[6]), as.numeric(x[7]), as.numeric(x[14]) , as.numeric(x[15])), nrow = 2, dimnames = list(Guess = c("Milk", "Tea"), Truth = c("Milk", "Tea")))
p=fisher.test(TeaTasting,alternative = "greater")$p.value
return(p)
}

fisher_test1=function(x){
TeaTasting = matrix(c(as.numeric(x[6]), as.numeric(x[7]), as.numeric(x[14]) , as.numeric(x[15])), nrow = 2, dimnames = list(Guess = c("Milk", "Tea"), Truth = c("Milk", "Tea")))
p=fisher.test(TeaTasting,alternative = "less")$p.value
return(p)
}
pval=apply(M,1,FUN=fisher_test)
M$greater=pval
pval1=apply(M,1,FUN=fisher_test1)
M$less=pval1

#file_new=paste0("/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-ATAC_",sample$V1[i],".filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP")
file_new=paste0("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/",sample$V1[i],".snATAC.merged.filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP")
write.table(M,file=file_new,sep="\t",quote=F)
}
}
