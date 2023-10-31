#chr:6:29027255  a       g       0.4933  0.0068  232207.00       5.977   2.271e-09       ++++    0.0     1.056   3       0.7876
#!/bin/sh
plink="/storage/chen/Software/plink_linux_x86_64/plink"

for i in {1..22}
do
#chr_file=sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_chr${i}
chr_file=sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_chr${i}
#1       8436802 rs6678140       T       C       -0.0564 0.0101  2.13e-08
awk -v b=${i} '{if($1==b){ print $1"\t"$2"\t"$2"\t"$1":"$2":"$3}}' sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8 > ${chr_file}
out=${chr_file}_1000G

input=/storage/chen/home/jw29/sc_human_retina/data/GWAS/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}
$plink --bfile $input --extract range ${chr_file}  --make-bed --out ${out}
done
