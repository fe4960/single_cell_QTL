#!/usr/bin/perl -w
use strict;
#my @sample_list = ("$ARGV[0]");
my $sample = $ARGV[0];
my $sample_region = $ARGV[2];
my $bam_head = "$ARGV[1]";

my $bam_dir = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam_new_merged";

`mkdir $bam_dir`;
#for my $sample(@sample_list){

my $BC_total = "sc_human_retina/data/ASE/ATAC_cell_archr_04062021";
my %hash;
open(INPUT,$BC_total);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
if($info[-2] < 0.5){ ###filter out the cell with predictedScore less than 0.5
next;
}
my @sample = split(/\#/,$info[0]);
my $samplename = "$sample[0]";
my $bc = $sample[1];
#$samplename =~ s/10xATAC_//g;

#if("$sample" eq "$samplename"){
if("$sample_region" eq "$samplename"){

my $label = "ATAC_"."$sample_region"."_$info[-1]";
#my $label = "ATAC_"."$sample"."_$info[-3]";

$hash{$label}->{$bc}=1;
}
}
for my $key (keys %hash){
my $output = "$bam_dir/$key";
open(OUTPUT,">$output");
for my $subkey (keys %{$hash{$key}}){
print OUTPUT "$subkey,$key\n";
}
my $bam = "$bam_head"."$sample"."/outs/possorted_bam.bam";
`cd $bam_dir`;
`/storage/chen/home/jw29/software/split_bam_bc_origin   $output  $bam `
}
#}
