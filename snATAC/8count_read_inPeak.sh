#!/bin/sh
dir="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/snATAC_featureCounts_combi/"
mkdir ${dir}
file="sc_human_retina/data/sample_list20"
#for cell in AC Cone Rod BC Astro MG HC RGC
for cell in $1
do
while read line 
do
bam="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/${cell}_${line}.snATAC.merged.filter.keep.merged.sorted.rmdup.filtered.bam"

software/subread-2.0.1-Linux-x86_64/bin/featureCounts -O -F SAF -g GeneID -a /storage/chen/home/jw29/sc_human_retina/data/ATAC_peak/narrowPeak_combi/${cell}_narrowPeak_combi_reform_macular_lobe  -o ${dir}/${line}_${cell}_merged_featureCount  ${bam}

done < $file
done 
