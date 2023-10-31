#!/bin/sh
#command=sc_human_retina/scripts/GWAS_new/finemapping/2run_torus_all.sh
command=sc_human_retina/scripts/GWAS_new/finemapping/3generate_r_plink_general
for file in clinical_c_Block_H30-H36 clinical_c_H33 clinical_c_H35 clinical_c_H36 selfReported_n_1301 
do
#nohup perl ${command}.pl sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1000g_flt_20ppl > ${command}_${file}.out 2> ${command}_${file}.err &
#sbatch --mem=5000MB --output=${command}_${file}.out --error=${command}_${file}.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1000g_flt_20ppl
sbatch --mem=5000MB --output=${command}_${file}.out --error=${command}_${file}.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format.vcf.1fpm.1000g_flt_20ppl

echo "$file"
done


#nohup perl ${command}.pl sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1000g_flt_20ppl > ${command}_Gharahkhani.out 2> ${command}_Gharahkhani.err &

#sbatch --mem=5000MB --output=${command}_Gharahkhani.out --error=${command}_Gharahkhani.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1000g_flt_20ppl
sbatch --mem=5000MB --output=${command}_Gharahkhani.out --error=${command}_Gharahkhani.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format.vcf.1fpm.1000g_flt_20ppl

#nohup perl ${command}.pl sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1000g_flt_20ppl > ${command}_Fritsche.out 2> ${command}_Fritsche.err &

#sbatch --mem=5000MB --output=${command}_Fritsche.out --error=${command}_Fritsche.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1000g_flt_20ppl
sbatch --mem=5000MB --output=${command}_Fritsche.out --error=${command}_Fritsche.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format.vcf.1fpm.1000g_flt_20ppl

#nohup perl ${command}.pl sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1000g_flt_20ppl > ${command}_Hysi.out 2> ${command}_Hysi.err &
#sbatch --mem=5000MB --output=${command}_Hysi.out --error=${command}_Hysi.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1000g_flt_20ppl
sbatch --mem=5000MB --output=${command}_Hysi.out --error=${command}_Hysi.err ${command}.sh sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format.vcf.1fpm.1000g_flt_20ppl



