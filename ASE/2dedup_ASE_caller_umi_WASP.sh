#!/bin/sh
bam=$1
vcf=$2
#gatk="/stornext/snfs5/ruichen/dmhg/fgi/Tao/Pipeline/pipeline/pipeline_restructure/bin/gatk4/gatk"
#gatk="/hgsc_software/gatk/gatk-4.1.0.0/gatk"
gatk="/storage/chen/Pipeline/pipeline/pipeline_restructure/bin/gatk4/gatk"
#samtools="/hgsc_software/samtools/samtools-1.9/bin/samtools"
#####${samtools} index ${bam}.bam
####${gatk} --java-options "-Xmx8g" MarkDuplicates -I ${bam}.bam -M ${bam}.metrics_file -O ${bam}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/stornext/snfs5/ruichen/dmhg/fgi/jwang/tmp
######${gatk} --java-options "-Xmx8g" MarkDuplicates -I ${bam}.bam -M ${bam}.metrics_file -O ${bam}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/storage/chen/home/jw29/tmp
#/stornext/snfs5/ruichen/dmhg/fgi/jwang/miniconda3/bin/umi_tools dedup --extract-umi-method=tag --umi-tag=UB -I ${oPath}/${sName}/${sName}.rmdup.bam -S ${oPath}/${sName}/${sName}.dedup.bam --method unique
####################!!!!!!samtools index ${bam}.bam
/storage/chen/Software/miniconda3/bin/umi_tools dedup --extract-umi-method=tag --umi-tag=UB -I ${bam}.bam -S ${bam}.rmdup.bam 

#/storage/chen/Software/miniconda3/bin/umi_tools dedup --extract-umi-method=tag --umi-tag=UB -I ${bam}.bam -S ${bam}.rmdup.bam --method unique
#${gatkPath} --java-options "-Xmx8g" MarkDuplicates -I ${oPath}/${sName}/BAM/${sName}.matefixed.sorted.bam -M ${oPath}/${sName}/BAM/${sName}.metrics_file -O ${oPath}/${sName}/BAM/${sName}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/storage/novaseq/tmp
##########/stornext/snfs5/ruichen/dmhg/fgi/jwang/miniconda3/bin/umi_tools dedup --extract-umi-method=tag --umi-tag=UB -I ${bam}.bam -S ${bam}.dedup.bam --method unique
#${gatk} AddOrReplaceReadGroups -I ${bam}.dedup.bam  -O ${bam}.dedup.rg.bam -RGID 4  -RGLB lib1  -RGPL ILLUMINA  -RGPU unit1  -RGSM 20
#grep -v "chrM" ${vcf}.vcf | grep -v "chrUn_" | grep -v "random" |grep -v "hap"  > ${vcf}.clean.vcf

#grep "^#" ${vcf}.vcf > ${vcf}.sort.vcf
######(grep "^#" ${vcf}.vcf && grep -v "^#" ${vcf}.vcf  | sort -k 1,1 -k 2n,2n) > ${vcf}.sort.vcf
# (grep "^#" ${vcf}.vcf && grep -v "^#" ${vcf}.vcf  | awk '{if((length($4)==1)&&(length($5)==1)){print }}')  > ${vcf}.snp.vcf
#(grep "^#" ${vcf}.vcf | grep -v "##contig=" && grep -v "^#" ${vcf}.vcf | sed -e 's/chr//g' | awk '{if((length($4)==1)&&(length($5)==1)&&($7=="PASS")){print }}')  > ${vcf}.snp.vcf

#${gatk} --java-options "-Xmx8g" IndexFeatureFile -F ${vcf}.snp.vcf

#${gatk} ASEReadCounter --variant ${vcf}.vcf --input ${bam}.rmdup.bam --output ${bam}.count -R /stornext/snfs5/ruichen/dmhg/fgi/jwang/sc_human_retina/sc_eQTL/hg19_ref/hg19_rmChrM.fasta
######${gatk} SelectVariants --select-type SNP -R /stornext/snfs5/ruichen/dmhg/fgi/jwang/sc_human_retina/sc_eQTL/hg38_ref/hg38.fa
samtools index ${bam}.rmdup.bam
samtools view -h ${bam}.rmdup.bam | awk  '($3 != "MT" && $3 !~ "GL")' | samtools view -Shb - > ${bam}.rmdup.filtered.bam
samtools index ${bam}.rmdup.filtered.bam
############samtools view -h ${bam}.rmdup.filtered.bam | LC_ALL=C grep -v AN:Z | LC_ALL=C grep -v RE:A:I  | samtools view -Shb -o ${bam}.rmdup.filtered.rmAntisenseIntron.bam
############samtools index ${bam}.rmdup.filtered.rmAntisenseIntron.bam
#######samtools view -h ${bam}.rmdup.filtered.bam | LC_ALL=C grep -v AN:Z | LC_ALL=C grep -v RE:A:I | LC_ALL=C grep -v RE:A:N  | samtools view -Shb -o ${bam}.rmdup.filtered.rmAntisenseKeepExon.bam
#######samtools index ${bam}.rmdup.filtered.rmAntisenseKeepExon.bam

#${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.bam --output ${bam}.snp.count -R /storage/chen/home/jw29/reference/hg38.fa
#${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.bam --output ${bam}.snp.count -R /storage/chen/Pipeline/pipeline/pipeline_restructure/reference/hg19/hg19.fasta
${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.filtered.bam --output ${bam}.snp.count -R   /storage/chen/home/jw29/Gript/hg19_rmchr.fasta
# ${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.filtered.rmAntisenseIntron.bam --output ${bam}.snp.rmAntiIntron.count -R   /storage/chen/home/jw29/Gript/hg19_rmchr.fasta
######${gatk} --java-options "-Xmx8g" ASEReadCounter --variant ${vcf}.snp.vcf --input ${bam}.rmdup.filtered.rmAntisenseKeepExon.bam --output ${bam}.snp.rmAntiKeepExon.count -R   /storage/chen/home/jw29/Gript/hg19_rmchr.fasta


#########!!!!! /storage/chen/Software/bcftools/bin/bcftools mpileup -Ou -f Gript/hg19_rmchr.fasta -d 100000000  ${bam}.rmdup.filtered.rmAntisenseKeepExon.bam  | /storage/chen/Software/bcftools/bin/bcftools  call -m -Ob -o ${bam}.mpileup.all.rmAntisenseKeepExon.bcf

#/storage/chen/Software/bcftools/bin/bcftools mpileup -Ou -f Gript/hg19_rmchr.fasta -d 100000000  -R sc_human_retina20/data/WGS20/D013_13.GATK.HaplotypeCaller.mark.snp.vcf_chr1_coverage10_pos ${sample}.rmdup.bam  | /storage/chen/Software/bcftools/bin/bcftools  call -m -Ob -o ${sample}.mpileup.chr1.bcf


#chrM    1612    .       A       G       111     PASS    DP=53;VDB=6.11476e-08;SGB=-0.691153;RPB=2.15592e-06;MQB=1;MQSB=1;BQB=0.952993;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=17,9,5,13;MQ=60        GT:PL   0/1:144,0,255
#/storage/chen/Software/bcftools/bin/bcftools filter  -g3 -G10 -e'%QUAL<30 ||  %MAX(DP)<10 ' ${sample}.mpileup.chr1.bcf  | bcftools view > ${sample}.mpileup.filtered.chr1.vcf

# awk '{if((length($4)==1)&&(length($5)==1)){split($8,a,";"); split(a[13],b,","); if((b[3] + b[4])>=2){print }}}' > ${1}.mpileup.filtered.vcf

##########!!!!!!/storage/chen/Software/bcftools/bin/bcftools filter  -g3 -G10 -e'%QUAL<30 ||  %MAX(DP)<10 ' ${bam}.mpileup.all.rmAntisenseKeepExon.bcf  | bcftools view > ${bam}.mpileup.filtered.all.mcfntisenseKeepExon.vcf


