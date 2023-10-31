#!/bin/sh
export PYTHONPATH=/storage/chen/Software/miniconda2/bin/python2.7:$PYTHONPATH
export LD_LIBRARY_PATH=/storage/chen/Software/miniconda2/lib/python2.7/site-packages/:$LD_LIBRARY_PATH
dir="/storage/chen/home/jw29/software"
#python=software/anaconda3/bin/python
#python=/storage/chen/Software/miniconda3/bin/python
python=/storage/chen/Software/miniconda2/bin/python2.7 
#vcf="/storage/chen/home/jw29/sc_human_retina/data/phase/1000GP_Phase3_20ppl_All.phased.wRef.merge.vcf.correct_ref_new.flt.gz"
file_dir="/storage/chen/home/jw29/software/phaser/useful_files/"
#out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/AC_phaser_pop"
#main="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/AC"
vcf="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/20ppl_phaser_all.vcf.gz"
#out=
#for cell in Cone ONBC OFFBC MG 
#for cell in AC RGC HC Rod
#for cell in AC
for cell in $1
do
#bed="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop"
#main="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}"
map="sc_human_retina/scripts/ASE/1phaser_cis_map_file_new_${cell}"

pair="sc_human_retina/data/ca_eQTL/phaser_map_${cell}_ens_rmX"
#bed="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop.gw_phased.bed.gz"
#out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop.gw_phased.cis_var.new_out1"

bed="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop_format_sort.gw_phased.bed.gz"
out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop_format_sort.gw_phased.cis_var.new_out"

#mkdir ${main}
$python ${dir}/phaser/phaser_pop/phaser_cis_var.py --bed ${bed} --vcf ${vcf} --pair ${pair} --map  ${map}  --o ${out} --ignore_v 1 --t 4 --chr ${2}
done 
