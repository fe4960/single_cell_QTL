#!/bin/sh
file="sc_human_retina/scripts/single_cell/snATAC/sample_list"
while read line
do
samtools sort -o ${line}_${1}.sorted.bam  -T ${line}_${1}.tmp.bam  ${line}_${1}.bam
done < $file
