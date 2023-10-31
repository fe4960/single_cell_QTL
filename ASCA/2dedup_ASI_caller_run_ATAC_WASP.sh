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
#for indi in D028_13 D027_13 D026_13 D021_13 D019_13 D018_13

#for cell in ONBC OFFBC Astro Rod MG Cone RGC HC AC
for cell in ONBC OFFBC Rod MG Cone RGC HC AC
#for cell in BC Astro
#for cell in Astro
#for cell in Rod MG Cone RGC HC AC
#for cell in Rod
do
label=${cell}_${indi}
#bam="/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam/out-cluster-snRNA_"${indi}"_"${cell}
##########bam=sc_human_retina20/data/snRNA_seq_bam/WASP_corrected_phased/out-cluster-snRNA_${indi}_${cell}.keep.merged.sorted
#######bam=sc_human_retina20/data/snATAC_seq_bam/WASP_corrected_phased20ppl/out-cluster-snRNA_${indi}_${cell}.keep.merged.sorted
#sc_human_retina20/data/snATAC_seq_bam/WASP_corrected_phased20ppl/out-cluster-ATAC_D028_13_Rod.filter.keep.merged.sorted.bam
#bam=sc_human_retina20/data/snATAC_seq_bam/WASP_corrected_phased20ppl/out-cluster-ATAC_${indi}_${cell}.filter.keep.merged.sorted
bam=/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/${cell}_${indi}.snATAC.merged.filter.keep.merged.sorted
#bam="/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam/out-cluster-ATAC_"${indi}"_"${cell}
#vcf="/storage/chen/home/jw29/sc_human_retina20/data/WGS20/"${indi}".GATK.HaplotypeCaller.mark_coor_hg19ToHg38"
vcf="/storage/chen/home/jw29/sc_human_retina20/data/WGS20/"${indi}".GATK.HaplotypeCaller.mark"
#echo "sh 2dedup_ASE_caller.sh $bam $vcf" | msub -q scavenger  -l nodes=1:ppn=1,pmem=30000MB -e ${dir}/${label}.err -o ${dir}/${label}.out -N ${label} -d $PWD
output=${dir}/${label}.out
err=${dir}/${label}.err
#sbatch --error=${err} --job-name=${label} --output=${output} --mem=40000MB /storage/chen/home/jw29/sc_human_retina/scripts/ASE/2dedup_ASE_caller.sh $bam $vcf 
sh /storage/chen/home/jw29/sc_human_retina/scripts/ASE/2dedup_ASE_caller.sh $bam $vcf 
#sh /storage/chen/home/jw29/sc_human_retina/scripts/ASE/2dedup_ASE_caller_umi_WASP.sh $bam $vcf 

done
done < $file

