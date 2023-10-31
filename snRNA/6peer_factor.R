library(peer)
#ct_name=c("ONBC","OFFBC","Astro","Rod","MG","Cone","RGC","HC","AC")
ct_name=c("BC","Astro","Rod","MG","Cone","RGC","HC","AC")
#ct_name=c("BC","Astro","Rod","MG","Cone","RGC","HC","AC","Mic")
for(ct in ct_name){
####ct_file = paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct)
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular_24pp_mean5cpm")
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_pool_macular_20pp")
#####ct_file = paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_pool_macular_20pp")
ct_file = paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_pool_macular_20pp_afterSoupX_new")
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/normalized_avg10cmp_",ct,"_pool_macular_20pp_afterSoupX")

exp=read.table(ct_file,header=T)

#exp=read.table("sc_human_retina20/data/snRNA_seq/normalized_cmp_AC",header=T)
dim(exp)
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(exp)))
dim(PEER_getPhenoMean(model))
#PEER_setNk(model,5) #5 factor
#PEER_setNk(model,1) #1 factor
#PEER_setNk(model,5) #2 factor
PEER_setNk(model,3) #2 factor
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)
dim(factors)
weights = PEER_getW(model)
dim(weights)
precision = PEER_getAlpha(model)
dim(precision)
residuals = PEER_getResiduals(model)
dim(residuals)
plot(precision)
PEER_setAdd_mean(model, TRUE)
rownames(factors)=colnames(exp)
#####write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_5peer_factor"),sep="\t",quote=F)
######write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular_24pp_mean5cpm_5peer_factor"),sep="\t",quote=F)
######write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular_20pp_mean5cpm_3peer_factor"),sep="\t",quote=F)
#######write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_macular_20pp_mean5cpm_3peer_factor"),sep="\t",quote=F)
write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_macular_20pp_mean5cpm_3peer_factor_afterSoupX_new"),sep="\t",quote=F)
#write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg10cmp_",ct,"_macular_20pp_mean10cpm_3peer_factor_afterSoupX"),sep="\t",quote=F)
#write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg10cmp_",ct,"_macular_20pp_mean10cpm_5peer_factor_afterSoupX"),sep="\t",quote=F)

#write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular_24pp_1peer_factor"),sep="\t",quote=F)
#write.table(factors,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular_24pp_2peer_factor"),sep="\t",quote=F)

}
