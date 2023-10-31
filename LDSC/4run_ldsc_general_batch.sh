#!/bin/sh

#err_dir="human_meta/scripts/GWAS/4run_ldsc_general_DAR"
#mkdir $err_dir


file="/storage/chenlab/Users/junwang/human_meta/data/cell_list"

#while read c
#do 

#sbatch --mem=5000MB -p interactive   --output=${err_dir}/${c}.out --error=${err_dir}/${c}.err human_meta/scripts/GWAS/4run_ldsc_general_DAR.sh $c

#done < $file

err_dir="human_meta/scripts/GWAS/4run_ldsc_general_all_peak"
mkdir $err_dir

while read c
do 

sbatch --mem=5000MB -p interactive   --output=${err_dir}/${c}.out --error=${err_dir}/${c}.err human_meta/scripts/GWAS/4run_ldsc_general_all_peak.sh $c

done < $file
