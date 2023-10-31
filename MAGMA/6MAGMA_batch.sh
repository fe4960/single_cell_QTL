#!/bin/sh

for file in MAGMA_list_5_0 MAGMA_list_5_1 MAGMA_list_5_2 MAGMA_list_5_3 
do 
sbatch --mem=10000MB -p interactive human_meta/scripts/GWAS/6MAGMA.sh $file

#sbatch --mem=10000MB -p interactive --nodelist=mhgcp-r03 human_meta/scripts/GWAS/6MAGMA.sh $file
done 
