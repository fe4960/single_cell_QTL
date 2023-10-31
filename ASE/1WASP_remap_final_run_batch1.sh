#!/bin/sh
#for i in {0..6}



#for i in 2
for i in 3 4
do
file="sc_human_retina/data/phase/sample_ceu_3_${i}"
race="ceu"
sh sc_human_retina/scripts/ASE/WASP_remap_final_run.sh ${file} ${race} 
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh ${file} ${race} 
#sh sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh  ${file} ${race}

done 

#file="sc_human_retina/data/phase/sample_oth"
#race="oth"
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_run.sh ${file} ${race}
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh ${file} ${race}

#file="sc_human_retina/data/phase/sample_eas"
#race="eas"
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_run.sh ${file} ${race}
#sbatch --mem=40000MB sc_human_retina/scripts/ASE/WASP_remap_final_ATAC_run.sh ${file} ${race}

