#!/bin/sh
#main="sc_human_retina/data/single_cell/snATAC"
#main="sc_human_retina/data/single_cell/snATAC/lobe_macular/"
export PYTHONPATH="/storage/chen/home/jw29/software/anaconda3/lib/python3.8/site-packages:/storage/chen/home/jw29/software/anaconda3/lib/python3.8:/storage/chen/home/jw29/software/anaconda3/bin"
main="/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_bam_new_merged/"
dir="/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_peak"
cd $dir
#cd $main
#for cell in RGC ONBC OFFBC
#for cell in ONBC
#for cell in Astro Cone HC 
#for cell in Astro Cone
for cell in Astro
#MG OFFBC RGC 
do
#sortname=${main}/${cell}.snATAC.merged.nameSorted.bam
sortname=${main}/${cell}.snATAC.merged.bam

/storage/chen/home/jw29/software/anaconda3/bin/macs3Env/bin/macs3 callpeak -f BAMPE -t ${sortname} -g hs -n ${cell}_macs3 -B -q 0.01 
done 
