#!/bin/sh
for cell in AC ONBC OFFBC MG HC RGC Rod Cone
#for cell in AC
do
#rm  /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/*
#mv /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/*_${cell}_phaser.gene_ae.txt /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/
#for sam in `ls /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/`
#do
#gunzip /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/${sam}.gz
#file=/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/${sam}
file="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop_new.gw_phased.bed.gz"
file_new="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop_new_sort.gw_phased.bed"
(less $file | head -n 1   && less $file | tail -n +2  | grep -v -P "MT|X|Y" | sort -k 1n,1n -k 2n,2n) > ${file_new}
bgzip ${file_new}
tabix -p bed ${file_new}.gz
#bgzip /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/${sam}
#tabix -p bed /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}/${sam}.gz
#done
done
