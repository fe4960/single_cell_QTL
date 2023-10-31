library(preprocessCore)
my.invnorm = function(x)
{
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
}

transform_standard_normal = function(df)
{
#    data_valid_expressed_full_qn = normalize.quantiles(as.matrix(df), copy=FALSE)
    data_valid_expressed_full_qn = normalize.quantiles(as.matrix(df), copy=FALSE)

    input_mat = as.data.frame(t(apply(t(data_valid_expressed_full_qn), 2, my.invnorm)))
    
    return(input_mat)
}

#normalize_tpm = function(tpm, min_tpm = 2, min_samples = 0.1,sample) # add PEER factors?
normalize_tpm = function(tpm, min_tpm = 1, min_samples = 0.1,sample) # add PEER factors?
{
    message("Normalizing expression...")
    message(paste("Total genes/peaks", nrow(tpm), sep = " = "))
    message(paste("Total samples"    , ncol(tpm), sep = " = "))
#    tpm=as.data.frame(tpm)    
    expressed = as.matrix(tpm)

#    expressed[as.matrix(tpm) <  min_tpm] = 0
#    expressed[as.matrix(tpm) >= min_tpm] = 1

#    tpm_f          = tpm[rowSums(expressed) >= (min_samples * ncol(tpm)),]
    tpm_f          = tpm[rowMeans(expressed) >= 5,]

    message(paste("Total expressed genes/peaks"    , nrow(tpm_f), sep = " = "))
    tpm_f_std_norm = transform_standard_normal(tpm_f)
    
#    return(list(tpm_f = tpm_f, tpm_f_std_norm = tpm_f_std_norm, gene_ids = rownames(tpm_f), sample_ids = colnames(tpm_f)))
return(tpm_f_std_norm)
}
#ct_name=c("ONBC","OFFBC","Astro","Rod","MG","Cone","RGC","HC","AC")
ct_name=c("BC","Astro","Rod","MG","Cone","RGC","HC","AC")
#ct_name = c("RGC","HC","AC")
#ct_name=c("BC","Astro","Rod","MG","Cone","RGC","HC","AC","Mic")
for(ct in ct_name){
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/",ct,"_snRNA_tmp_ave")
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/",ct,"_snRNA_tmp_ave_macular")
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/",ct,"_snRNA_tmp_ave_macular_24pp")
########ct_file = paste0("sc_human_retina20/data/snRNA_seq/",ct,"_snRNA_pool_ave_macular_20pp")
#ct_file = paste0("sc_human_retina20/data/snRNA_seq/",ct,"_snRNA_pool_ave_macular_20pp_afterSoupX")
ct_file = paste0("sc_human_retina20/data/snRNA_seq/",ct,"_snRNA_pool_ave_macular_20pp_afterSoupX_new")

#ct_exp=read.table(ct_file, header=T, na.strings = 0)

ct_exp=read.table(ct_file,header=T)
tpm_f_std_norm=normalize_tpm(ct_exp,min_tpm = 1, min_samples = 0.1,sample=ct)
print(ct)
#write.table(tpm_f_std_norm,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct),quote=F,sep="\t")
#write.table(tpm_f_std_norm,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular"),quote=F,sep="\t")
#write.table(tpm_f_std_norm,paste0("sc_human_retina20/data/snRNA_seq/normalized_cmp_",ct,"_macular_24pp"),quote=F,sep="\t")
#####write.table(tpm_f_std_norm,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_pool_macular_20pp"),quote=F,sep="\t", row.names=T)
#write.table(tpm_f_std_norm,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_pool_macular_20pp_afterSoupX"),quote=F,sep="\t", row.names=T)
write.table(tpm_f_std_norm,paste0("sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_",ct,"_pool_macular_20pp_afterSoupX_new"),quote=F,sep="\t", row.names=T)

}

