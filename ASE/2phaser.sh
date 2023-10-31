#!/bin/sh
export PYTHONPATH=/storage/chen/Software/miniconda2/bin/python2.7:$PYTHONPATH
export LD_LIBRARY_PATH=/storage/chen/Software/miniconda2/lib/python2.7/site-packages/:$LD_LIBRARY_PATH
dir="/storage/chen/home/jw29/software/"
python=/storage/chen/Software/miniconda2/bin/python2.7 
sample=$1
sample1=$2
cell=$3
label=${sample1}_${cell}
#bam="sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_Rod.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_Cone.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_AC.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_ONBC.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_OFFBC.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_HC.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_MG.keep.merged.sorted.bam,sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_RGC.keep.merged.sorted.bam"
#vcf="/storage/chen/home/jw29/sc_human_retina/data/phase/1000GP_Phase3_20ppl_chr${chr}.phased.wRef.merge.vcf.correct_ref_new_flt_noChr.gz"
bam="sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-snRNA_${sample1}_${cell}.keep.merged.sorted.bam"
vcf="/storage/chen/home/jw29/sc_human_retina/data/phase/1000GP_Phase3_20ppl_All.phased.wRef.merge.vcf.correct_ref_new.flt.gz"
file_dir="/storage/chen/home/jw29/software/phaser/useful_files/"
out="/storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/${label}_phaser"
mkdir /storage/chen/home/jw29/sc_human_retina/data/ASE/phaser/
$python ${dir}/phaser/phaser/phaser.py --vcf ${vcf} --bam ${bam} --paired_end 0 --mapq 255 --baseq 10 --sample ${sample} --blacklist ${file_dir}/hg19_hla.bed --haplo_count_blacklist ${file_dir}/hg19_haplo_count_blacklist.bed --threads 4 --o ${out} --gw_phase_vcf 1 --gw_phase_method 0 
$python ${dir}/phaser/phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts ${out}.haplotypic_counts.txt --features ${file_dir}/gencode.v19.GRCh37.genes.bed --o ${out}.gene_ae.txt

