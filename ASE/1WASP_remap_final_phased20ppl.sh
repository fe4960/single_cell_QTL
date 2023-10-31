#!/bin/sh
# Set these environment vars to point to
# your local installation of WASP
s="out-cluster-snRNA_${1}"
#s=$1
indi=$2
#indi1=$(echo $indi | cut -d "-" | awk '{print $1$2}')
indi1=${indi/_/}
#echo $indi1
#exit
gt=genotypes_${3}_ref_phased
WASP=/storage/chen/home/jw29/software/WASP
DATA_DIR1=/storage/chen/home/jw29/sc_human_retina/data/ASE/new_analysis
DATA_DIR=/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl/
mkdir ${DATA_DIR}
bam_dir=/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam20ppl
/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/find_intersecting_snps.py \
       --output_dir $DATA_DIR  \
       --snp_index $DATA_DIR1/${gt}/snp_index.h5 \
       --snp_tab $DATA_DIR1/${gt}/snp_tab.h5 \
       --haplotype $DATA_DIR1/${gt}/haps.h5 \
       --samples $indi1 \
        ${bam_dir}/${s}.bam


#exit
#STAR   --runThreadN 4   --genomeDir /storage/chen/Reference/Ref_cellranger.hg19.premrna/star       --readNameSeparator space      --outStd SAM   --outSAMtype SAM      --outSAMunmapped Within   KeepPairs      --outSAMorder PairedKeepInputOrder --outSAMattrRGline ID:10x3_D026_13:0:1:HNGHTDSXX:1 SM:10x3_D026_13 LB:0.1 PU:10x3_D026_13:0:1:HNGHTDSXX:1 PL:ILLUMINA --readFilesIn ${DATA_DIR}/${s}.remap.single.fq.gz   | samtools view -S -b -q 10 - >  ${DATA_DIR}/${s}.remap.bam
STAR   --runThreadN 4   --genomeDir /storage/chen/Reference/Ref_cellranger.hg19.premrna/star       --readNameSeparator space      --outStd SAM   --outSAMtype SAM    --readFilesCommand zcat   --outSAMattrRGline ID:10x3_${indi}:0:1:HNGHTDSXX:1 SM:10x3_${indi} LB:0.1 PU:10x3_${indi}:0:1:HNGHTDSXX:1 PL:ILLUMINA --readFilesIn ${DATA_DIR}/${s}.remap.fq.gz --outTmpDir ${DATA_DIR}/${s}_STARtmp  | samtools view -S -b -q 10 - >  ${DATA_DIR}/${s}.remap.bam

# Use filter_remapped_reads.py to create filtered list of reads that correctly
# remap to same position
/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/filter_remapped_reads.py \
       $DATA_DIR/${s}.to.remap.bam \
       $DATA_DIR/${s}.remap.bam \
       $DATA_DIR/${s}.remap.keep.bam

# Create a merged BAM containing [1] reads that did
# not need remapping [2] filtered remapped reads
samtools merge $DATA_DIR/${s}.keep.merged.bam \
	 $DATA_DIR/${s}.keep.bam $DATA_DIR/${s}.remap.keep.bam

# Sort and index the bam file

samtools view $DATA_DIR/${s}.keep.merged.bam -o $DATA_DIR/${s}.keep.merged.noheader.sam
samtools view -H $DATA_DIR/${s}.keep.merged.bam | grep -v "@CO" > $DATA_DIR/${s}.keep.merged.header

cat $DATA_DIR/${s}.keep.merged.header $DATA_DIR/${s}.keep.merged.noheader.sam > $DATA_DIR/${s}.keep.merged.reheader.sam

#samtools sort $DATA_DIR/${s}.keep.merged.bam \
#	 -o $DATA_DIR/${s}.keep.merged.sorted.bam
samtools sort $DATA_DIR/${s}.keep.merged.reheader.sam -@ 4\
	 -o $DATA_DIR/${s}.keep.merged.sorted.bam

	 
samtools index $DATA_DIR/${s}.keep.merged.sorted.bam

# Filter out duplicate reads. Use rmdup_pe.py for paired-end reads,
# rmdup.py for single-end reads.
#/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/rmdup_pe.py $DATA_DIR/${s}.keep.merged.sorted.bam \
#       $DATA_DIR/${s}.keep.rmdup.merged.sorted.bam
#/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/rmdup.py $DATA_DIR/${s}.keep.merged.sorted.bam \
#       $DATA_DIR/${s}.keep.rmdup.merged.sorted.bam

rm $DATA_DIR/${s}.keep.merged.noheader.sam
rm $DATA_DIR/${s}.keep.merged.header
rm $DATA_DIR/${s}.keep.merged.reheader.sam
