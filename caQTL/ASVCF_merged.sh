#!/bin/sh
export RASQUALDIR="/storage/chen/home/jw29/software/rasqual/"
chr=$1
cell=$2
bam_list=/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/bam_list_${cell}_merged
sh /storage/chen/home/jw29/software/rasqual/src/ASVCF/createASVCF.sh ${bam_list} sc_human_retina/data/ca_eQTL/20ppl_${cell}/1000GP_Phase3_20ppl_chr${chr}.phased.wRef.merge.vcf.correct_ref_new.flt20ppl.gz  sc_human_retina/scripts/ca_eQTL/tmp_RASQUAL/tmp_merged.${cell}.${chr}.sh  sc_human_retina/data/ca_eQTL/tmp/tmp_merged.${cell}.${chr} 
