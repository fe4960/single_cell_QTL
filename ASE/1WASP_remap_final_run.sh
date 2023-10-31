#!/bin/sh
dir="/storage/chen/home/jw29/sc_human_retina20/err_dir"
mkdir $dir
#for indi in  19_D006 19_D007 19_D008 19_D009
#for indi in 19_D010 19_D011 19_D019 D005_13 D009_13
#for indi in D013_13 D017_13 D018_13 D019_13 
#for indi  in D021_13 D026_13 
#for indi in D026_13 
file=$1
while read line
do 
indi=${line}
#D013_13 D017_13 D018_13 D019_13 D021_13 D026_13 D027_13 D028_13 
#for indi in 19_D003
#for indi in D013_13  D009_13 D005_13  D030_13 19_D019 19_D011 19_D010
#for indi in D028_13 D013_13

#for cell in Rod
for cell in ONBC OFFBC Astro Rod MG Cone RGC HC AC

#for cell in BC Astro Rod MG Cone RGC HC AC
#for cell in Astro
#for cell in BC Astro
#for cell in Astro
#for cell in Rod MG Cone RGC HC AC
do
race=$2
label=${indi}_${cell}
bam="/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam20ppl/out-cluster-snRNA_"${indi}"_"${cell}

#bam="/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam/out-cluster-snRNA_"${indi}"_"${cell}
#bam="/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam/out-cluster-ATAC_"${indi}"_"${cell}
#vcf="/storage/chen/home/jw29/sc_human_retina20/data/WGS20/"${indi}".GATK.HaplotypeCaller.mark_coor_hg19ToHg38"
vcf="/storage/chen/home/jw29/sc_human_retina20/data/WGS20/"${indi}".GATK.HaplotypeCaller.mark"
#echo "sh 2dedup_ASE_caller.sh $bam $vcf" | msub -q scavenger  -l nodes=1:ppn=1,pmem=30000MB -e ${dir}/${label}.err -o ${dir}/${label}.out -N ${label} -d $PWD
output=${dir}/${label}.out
err=${dir}/${label}.err
#sbatch --error=${err} --job-name=${label} --output=${output} --mem=40000MB /storage/chen/home/jw29/sc_human_retina/scripts/ASE/2dedup_ASE_caller.sh $bam $vcf 
#sh /storage/chen/home/jw29/sc_human_retina/scripts/ASE/2dedup_ASE_caller.sh $bam $vcf 
#sh /storage/chen/home/jw29/sc_human_retina/scripts/ASE/2dedup_ASE_caller_umi.sh $bam $vcf 
#sh sc_human_retina/scripts/ASE/WASP_remap_final.sh  $label $indi 
sh sc_human_retina/scripts/ASE/WASP_remap_final_phased20ppl.sh  $label $indi  $race
#sh sc_human_retina/scripts/ASE/WASP_remap_final_phased_rest.sh  $label $indi  $race

done
done < $file

