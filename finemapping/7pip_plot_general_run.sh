#!/bin/sh
#command=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all.sh
#command=sc_human_retina/scripts/GWAS_new/finemapping/3generate_r_plink_general
command=sc_human_retina/scripts/GWAS_new/finemapping/5pip_plot_general
err=${command}_err
mkdir $err
for file in clinical_c_Block_H30-H36 clinical_c_H33 clinical_c_H35 clinical_c_H36 selfReported_n_1301 
do
#nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1000g_flt_20ppl > ${command}_${file}.out 2> ${command}_${file}.err &
sbatch -p gpu --mem=5000MB ${command}.sh sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1fpm.1000g_flt_20ppl --output=${err}/${file}.out --error=${err}/${file}.err 

# software/R-4.0.0/bin/Rscript --vanilla  ${command}.R sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1fpm.1000g_flt_20ppl > ${command}_${file}.out 2> ${command}_${file}.err &

done

###nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R  sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1fpm.1000g_flt_20ppl > ${command}_Gharahkhani.out 2> ${command}_Gharahkhani.err &

file="sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1fpm.1000g_flt_20ppl" 
sbatch -p gpu --mem=5000MB  ${command}.sh $file  --output=${err}/${file}.out --error=${err}/${file}.err

#nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R  sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1000g_flt_20ppl > ${command}_Gharahkhani.out 2> ${command}_Gharahkhani.err &

#nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1000g_flt_20ppl > ${command}_Fritsche.out 2> ${command}_Fritsche.err &
file="sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1fpm.1000g_flt_20ppl"
sbatch -p gpu --mem=5000MB ${command}.sh $file  --output=${err}/${file}.out --error=${err}/${file}.err


#nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1fpm.1000g_flt_20ppl > ${command}_Fritsche.out 2> ${command}_Fritsche.err &

#nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1000g_flt_20ppl > ${command}_Hysi.out 2> ${command}_Hysi.err &

file="sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1fpm.1000g_flt_20ppl"
sbatch -p gpu --mem=5000MB ${command}.sh $file  --output=${err}/${file}.out --error=${err}/${file}.err

#nohup software/R-4.0.0/bin/Rscript --vanilla  ${command}.R sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1fpm.1000g_flt_20ppl > ${command}_Hysi.out 2> ${command}_Hysi.err &



