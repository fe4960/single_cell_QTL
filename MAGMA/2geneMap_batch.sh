#!/bin/sh

#for file in 2ctd_list_4_0 2ctd_list_4_1 2ctd_list_4_2 2ctd_list_4_3 2ctd_list_4_4 2ctd_list_4_5
#for file in 2ctd_list_4_5 2ctd_list_4_4 
#for file in 2ctd_list_10_0 2ctd_list_10_1
#for file in 2ctd_list_10_1 2ctd_list_10_0_1 2ctd_list_10_0_2  2ctd_list_10_0_3  2ctd_list_10_0_4  2ctd_list_10_0_5 2ctd_list_10_0_6
for file in 2ctd_list_new1
do
sbatch -p interactive --mem=10000MB --nodelist=mhgcp-r01 human_meta/scripts/GWAS/2geneMap.sh $file
done 

#####for file in 2ctd_list_10_0_3  2ctd_list_10_0_4
#####do
#####sbatch --mem=10000MB --nodelist=mhgcp-c01 human_meta/scripts/GWAS/2geneMap.sh $file
#####done 

####for file in 2ctd_list_10_0_5 2ctd_list_10_0_6
####do
####sbatch --mem=10000MB --nodelist=mhgcp-c02 human_meta/scripts/GWAS/2geneMap.sh $file
#####done 



#sbatch --mem=5000MB -p short --nodelist=mhgcp-d00 human_meta/scripts/GWAS/2geneMap.sh 2ctd_list_4_3

#for file in 2ctd_list_4_2 
#do
#sbatch --mem=5000MB -p short --nodelist=mhgcp-d01 human_meta/scripts/GWAS/2geneMap.sh $file
#done

#sbatch --mem=5000MB -p short --nodelist=mhgcp-d03 human_meta/scripts/GWAS/2geneMap.sh 2ctd_list_4_0

#sbatch --mem=5000MB -p short --nodelist=mhgcp-d01 human_meta/scripts/GWAS/2geneMap.sh 2ctd_list_4_1

