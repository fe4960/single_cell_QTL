#!/bin/sh
command=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all.sh
########for file in clinical_c_Block_H30-H36 clinical_c_H33 clinical_c_H35 clinical_c_H36 selfReported_n_1301 
#######do
#nohup sh $command sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1000g_flt_20ppl > sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_${file}.out 2> sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_${file}.err &
#sbatch --mem=5000MB --output=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_${file}.out --error=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_${file}.err $command sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1000g_flt_20ppl
###########sbatch --mem=5000MB --output=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_${file}.out --error=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_${file}.err $command sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1fpm.1000g_flt_20ppl

########done


#nohup sh $command sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1000g_flt_20ppl > sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Gharahkhani.out 2> sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Gharahkhani.err &

#############sbatch --mem=5000MB --output=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Gharahkhani.out --error=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Gharahkhani.err $command sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1fpm.1000g_flt_20ppl 

#nohup sh $command sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1000g_flt_20ppl > sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Fritsche.out 2> sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Fritsche.err &
sbatch --mem=5000MB --output=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Fritsche_20ppl.out --error=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Fritsche_20ppl.err $command sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1fpm.1000g_flt_20ppl 

#nohup sh $command sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1000g_flt_20ppl > sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Hysi.out 2> sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Hysi.err &
sbatch --mem=5000MB --output=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Hysi_20ppl.out --error=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all_Hysi_20ppl.err $command sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1fpm.1000g_flt_20ppl 


