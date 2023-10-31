#!/bin/sh
export PYTHONPATH=/storage/chen/Software/miniconda2/bin/python2.7:$PYTHONPATH
export LD_LIBRARY_PATH=/storage/chen/Software/miniconda2/lib/python2.7/site-packages/:$LD_LIBRARY_PATH
dir="/storage/chen/home/jw29/software"
python=/storage/chen/Software/miniconda2/bin/python2.7 
#python=software/anaconda3/bin/python
#python=/storage/chen/Software/miniconda3/bin/python
vcf="/storage/chen/home/jw29/sc_human_retina/data/phase/1000GP_Phase3_20ppl_All.phased.wRef.merge.vcf.correct_ref_new.flt.gz"
file_dir="/storage/chen/home/jw29/software/phaser/useful_files/"
#out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/AC_phaser_pop"
#main="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/AC"
for cell in Cone ONBC OFFBC MG AC RGC HC Rod

do
#out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop_new"
out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${cell}_phaser_pop_python2"
main="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/gene_ae/${cell}"

mkdir ${main}
$python ${dir}/phaser/phaser_pop/phaser_expr_matrix.py --gene_ae_dir ${main} --features ${file_dir}/gencode.v19.GRCh37.genes.bed  --o ${out}
done 
