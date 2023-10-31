#!/bin/sh
#for i in {0..6}



#for i in 1
#for i in 0 1 
for i in 0_1
do
file="sc_human_retina/data/phase/sample_ceu_3_${i}"
race="ceu"
#sh sc_human_retina/scripts/ASE/WASP_remap_final_run.sh ${file} ${race} 
out=sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run_ceu_${i}.out
err=sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run_ceu_${i}.err
#sbatch --mem=40000MB --output=${out} --error=${err} sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh ${file} ${race} 
nohup sh sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh  ${file} ${race} > ${out} 2> ${err} &
done 

#file="sc_human_retina/data/phase/sample_oth"
#race="oth"
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_run.sh ${file} ${race}
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh ${file} ${race}

#file="sc_human_retina/data/phase/sample_eas"
#race="eas"
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_run.sh ${file} ${race}
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh ${file} ${race}

