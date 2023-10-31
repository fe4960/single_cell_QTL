#!/bin/sh
bam=$1
vcf=$2
#gatk="/stornext/snfs5/ruichen/dmhg/fgi/Tao/Pipeline/pipeline/pipeline_restructure/bin/gatk4/gatk"
#gatk="/hgsc_software/gatk/gatk-4.1.0.0/gatk"
gatk="/storage/chen/Pipeline/pipeline/pipeline_restructure/bin/gatk4/gatk"
#####samtools="/hgsc_software/samtools/samtools-1.9/bin/samtools"
#####${samtools} index ${bam}.bam
####${gatk} --java-options "-Xmx8g" MarkDuplicates -I ${bam}.bam -M ${bam}.metrics_file -O ${bam}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/stornext/snfs5/ruichen/dmhg/fgi/jwang/tmp
${gatk} --java-options "-Xmx8g" MarkDuplicates -I ${bam}.bam -M ${bam}.metrics_file -O ${bam}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/storage/chen/home/jw29/tmp

#${gatkPath} --java-options "-Xmx8g" MarkDuplicates -I ${oPath}/${sName}/BAM/${sName}.matefixed.sorted.bam -M ${oPath}/${sName}/BAM/${sName}.metrics_file -O ${oPath}/${sName}/BAM/${sName}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/storage/novaseq/tmp
##########/stornext/snfs5/ruichen/dmhg/fgi/jwang/miniconda3/bin/umi_tools dedup --extract-umi-method=tag --umi-tag=UB -I ${bam}.bam -S ${bam}.dedup.bam --method unique
#${gatk} AddOrReplaceReadGroups -I ${bam}.dedup.bam  -O ${bam}.dedup.rg.bam -RGID 4  -RGLB lib1  -RGPL ILLUMINA  -RGPU unit1  -RGSM 20
#grep -v "chrM" ${vcf}.vcf | grep -v "chrUn_" | grep -v "random" |grep -v "hap"  > ${vcf}.clean.vcf

#grep "^#" ${vcf}.vcf > ${vcf}.sort.vcf
######(grep "^#" ${vcf}.vcf && grep -v "^#" ${vcf}.vcf  | sort -k 1,1 -k 2n,2n) > ${vcf}.sort.vcf
########(grep "^#" ${vcf}.vcf && grep -v "^#" ${vcf}.vcf  | awk '{if((length($4)==1)&&(length($5)==1)){print }}') > ${vcf}.snp.vcf
######${gatk} --java-options "-Xmx8g" IndexFeatureFile -F ${vcf}.snp.vcf

#${gatk} ASEReadCounter --variant ${vcf}.vcf --input ${bam}.rmdup.bam --output ${bam}.count -R /stornext/snfs5/ruichen/dmhg/fgi/jwang/sc_human_retina/sc_eQTL/hg19_ref/hg19_rmChrM.fasta
######${gatk} SelectVariants --select-type SNP -R /stornext/snfs5/ruichen/dmhg/fgi/jwang/sc_human_retina/sc_eQTL/hg38_ref/hg38.fa

#######${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.bam --output ${bam}.snp.count -R /storage/chen/home/jw29/reference/hg38.fa

#@SQ     SN:GL000192.1   LN:547496
#@SQ     SN:NC_007605    LN:171823
#@SQ     SN:s37d5        LN:35477943

samtools index ${bam}.rmdup.bam
#samtools view -h ${bam}.rmdup.bam | awk  '(($3 != "s37d5") && ($3 !~ "GL")  && ($3 != "NC_007605") && ($3!="chrM") )' | samtools view -Shb - > ${bam}.rmdup.filtered.bam
samtools view -h ${bam}.rmdup.bam | awk  '(($3 != "s37d5") && ($3 !~ "GL")  && ($3 != "NC_007605") && ($3!="chrM") )' | sed -e 's/LN:16569/LN:16571/' | samtools view -Shb - > ${bam}.rmdup.filtered.bam
samtools index ${bam}.rmdup.filtered.bam

#${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.bam --output ${bam}.snp.count -R /storage/chen/home/jw29/Gript/hg19.fasta
${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.chr.vcf --input ${bam}.rmdup.filtered.bam --output ${bam}.snp.count -R /storage/chen/home/jw29/Gript/hg19.fasta
