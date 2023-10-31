#!/bin/sh
grep -v s37d5 /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed_flt_1fpm | grep -v GL > /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_flt_1fpm 


bedtools intersect -wao -a /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_flt_1fpm  -b sc_human_retina/data/ca_eQTL/wgEncodeHg19ConsensusSignalArtifactRegions.sorted.bed  > /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_overlap_Blacklist_flt_1fpm
awk '{if($NF==0){a=$3-$2; print $1"\t"$2"\t"$3"\t"a}}' /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_overlap_Blacklist_flt_1fpm > /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_rmBlacklist_flt_1fpm

grep -v "chrY" /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_rmBlacklist_flt_1fpm > /storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_rmBlacklist_rmY_flt_1fpm


