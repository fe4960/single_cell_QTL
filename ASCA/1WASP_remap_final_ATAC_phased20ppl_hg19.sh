#!/bin/sh
# Set these environment vars to point to
# your local installation of WASP
#s="out-cluster-snRNA_${1}"
s="out-cluster-ATAC_${1}.filter"
s1="out-cluster-ATAC_${1}"

#s=$1
indi=$2
#@PG     PN:bwa  ID:bwa  VN:0.7.17-r1188 CL:bwa mem -p -t 4 -M -R @RG\tID:10xATAC_D028_13:MissingLibrary:1:HNGHLDSXX:1\tSM:10xATAC_D028_13\tLB:MissingLibrary.1\tPU:10xATAC_D028_13:MissingLibrary:1:HNGHLDSXX:1\tPL:ILLUMINA /storage/chen/Reference/refdata-cellranger-atac-GRCh38-1.1.0/fasta/genome.fa /storage/novaseq/Data/200317_A00431_0172_AHNGHLDSXX/10xATAC_D028_13/SC_ATAC_COUNTER_CS/SC_ATAC_COUNTER/_BASIC_SC_ATAC_COUNTER/_ALIGNER/TRIM_READS/fork0/chnk0-ud584878d65/files/read1.fastq
shellPath=/storage/chen/Pipeline/pipeline/pipeline_restructure
export PATH=${shellPath}/bin/jdk1.8.0_121/bin/:$PATH
bwaPath=${shellPath}/bin/bwa

gt=genotypes_${3}_ref_phased

#gt=genotypes_${3}_ref_phased_hg38
WASP=/storage/chen/home/jw29/software/WASP
DATA_DIR1=/storage/chen/home/jw29/sc_human_retina/data/ASE/new_analysis
####DATA_DIR=/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/
DATA_DIR=/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/

mkdir ${DATA_DIR}
#bam_dir=/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam
#####bam_dir=/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/
bam_dir=/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/

indi1=${indi/_/}
#bam_dir0=/storage/chen/home/jw29/ ###newly add

#samtools view -h ${bam_dir}/${s1}.bam | awk '{if($6 !="*"){print}}' | samtools view -bS - > ${bam_dir}/${s}.bam
samtools view -h ${bam_dir}/${s1}.bam | awk '{if($6 !="*"){print}}' | samtools view -bS - > ${bam_dir}/${s}.bam

/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/find_intersecting_snps.py \
	--is_paired_end \
       --output_dir $DATA_DIR  \
       --snp_index $DATA_DIR1/${gt}/snp_index.h5 \
       --snp_tab $DATA_DIR1/${gt}/snp_tab.h5 \
       --haplotype $DATA_DIR1/${gt}/haps.h5 \
       --samples $indi1 \
        ${bam_dir}/${s}.bam





#/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/find_intersecting_snps.py \
####      --is_paired_end --is_sorted \
####       --output_dir $DATA_DIR  \
###       --snp_index $DATA_DIR1/genotypes_hg38/snp_index.h5 \
###       --snp_tab $DATA_DIR1/genotypes_hg38/snp_tab.h5 \
###       --haplotype $DATA_DIR1/genotypes_hg38/haps.h5 \
###       --samples $indi \
####        ${bam_dir}/${s}.bam


###exit
#       --is_sorted \

#STAR   --runThreadN 4   --genomeDir /storage/chen/Reference/Ref_cellranger.hg19.premrna/star       --readNameSeparator space      --outStd SAM   --outSAMtype SAM      --outSAMunmapped Within   KeepPairs      --outSAMorder PairedKeepInputOrder --outSAMattrRGline ID:10x3_D026_13:0:1:HNGHTDSXX:1 SM:10x3_D026_13 LB:0.1 PU:10x3_D026_13:0:1:HNGHTDSXX:1 PL:ILLUMINA --readFilesIn ${DATA_DIR}/${s}.remap.single.fq.gz   | samtools view -S -b -q 10 - >  ${DATA_DIR}/${s}.remap.bam
#STAR   --runThreadN 4   --genomeDir /storage/chen/Reference/Ref_cellranger.hg19.premrna/star       --readNameSeparator space      --outStd SAM   --outSAMtype SAM    --readFilesCommand zcat   --outSAMattrRGline ID:10x3_${indi}:0:1:HNGHTDSXX:1 SM:10x3_${indi} LB:0.1 PU:10xATAC_${indi}:MissingLibrary:1:HNGHLDSXX:1 PL:ILLUMINA --readFilesIn ${DATA_DIR}/${s}.remap.fq.gz   | samtools view -S -b -q 10 - >  ${DATA_DIR}/${s}.remap.bam

#######${bwaPath} mem -R "@RG\tID:10xATAC_${indi}:MissingLibrary:1:HNGHLDSXX:1\tSM:10xATAC_${indi}\tLB:MissingLibrary.1\tPU:10xATAC_${indi}:MissingLibrary:1:HNGHLDSXX:1\tPL:ILLUMINA" -t 4 -M /storage/chen/Reference/refdata-cellranger-atac-GRCh38-1.1.0/fasta/genome.fa  ${DATA_DIR}/${s}.remap.fq1.gz ${DATA_DIR}/${s}.remap.fq2.gz | samtools  view -bS -o ${DATA_DIR}/${s}.remap.bam -

${bwaPath} mem -R "@RG\tID:10xATAC_${indi}:MissingLibrary:1:HNGHLDSXX:1\tSM:10xATAC_${indi}\tLB:MissingLibrary.1\tPU:10xATAC_${indi}:MissingLibrary:1:HNGHLDSXX:1\tPL:ILLUMINA" -t 4 -M /storage/chen/Reference/refdata-cellranger-atac-hg19-1.2.0/fasta/genome.fa  ${DATA_DIR}/${s}.remap.fq1.gz ${DATA_DIR}/${s}.remap.fq2.gz | samtools  view -bS -o ${DATA_DIR}/${s}.remap.bam -


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

#samtools view $DATA_DIR/${s}.keep.merged.bam -o $DATA_DIR/${s}.keep.merged.noheader.sam
#samtools view -H $DATA_DIR/${s}.keep.merged.bam | grep -v "@CO" > $DATA_DIR/${s}.keep.merged.header

#cat $DATA_DIR/${s}.keep.merged.header $DATA_DIR/${s}.keep.merged.noheader.sam > $DATA_DIR/${s}.keep.merged.reheader.sam

software/samtools-1.11/samtools sort $DATA_DIR/${s}.keep.merged.bam -@ 4 \
	 -o $DATA_DIR/${s}.keep.merged.sorted.bam
#samtools sort $DATA_DIR/${s}.keep.merged.reheader.sam -@ 4\
#	 -o $DATA_DIR/${s}.keep.merged.sorted.bam

	 
samtools index $DATA_DIR/${s}.keep.merged.sorted.bam

# Filter out duplicate reads. Use rmdup_pe.py for paired-end reads,
# rmdup.py for single-end reads.
#/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/rmdup_pe.py $DATA_DIR/${s}.keep.merged.sorted.bam \
#       $DATA_DIR/${s}.keep.rmdup.merged.sorted.bam
#/storage/chen/Software/miniconda3/bin/python3 $WASP/mapping/rmdup.py $DATA_DIR/${s}.keep.merged.sorted.bam \
#       $DATA_DIR/${s}.keep.rmdup.merged.sorted.bam

#rm $DATA_DIR/${s}.keep.merged.noheader.sam
#rm $DATA_DIR/${s}.keep.merged.header
#rm $DATA_DIR/${s}.keep.merged.reheader.sam
#gatk="/storage/chen/Pipeline/pipeline/pipeline_restructure/bin/gatk4/gatk"
#####samtools="/hgsc_software/samtools/samtools-1.9/bin/samtools"
#####${samtools} index ${bam}.bam
####${gatk} --java-options "-Xmx8g" MarkDuplicates -I ${bam}.bam -M ${bam}.metrics_file -O ${bam}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/stornext/snfs5/ruichen/dmhg/fgi/jwang/tmp
#${gatk} --java-options "-Xmx8g" MarkDuplicates -I $DATA_DIR/${s}.keep.merged.sorted.bam -M ${bam}.metrics_file -O $DATA_DIR/${s}.keep.merged.sorted.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/storage/chen/home/jw29/tmp



