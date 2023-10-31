#!/usr/bin/perl -w
use strict;
my %hash=();
#my $file_list = "sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl_snp_list";
my $file_list = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/phased20ppl_merged_list";
#for my $file_list (@file_list){
open(INPUT,$file_list);
while(my $line = <INPUT>){
chomp $line;
#19_D003_AC.snp.count
$line =~ s/\.snp\.count//g;
my @id = split(/\_/,$line);
my $indi = $id[1]."_".$id[2];
my $cell = $id[0];
#my $input = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-ATAC_$line".".filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP";
my $input = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/$line".".snATAC.merged.filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP";

open(INPUT1,$input); # ||die(print $input);
<INPUT1>;
#print OUTPUT "$line1";
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
shift(@info1);
#contig  position        variantID       refAllele       altAllele       refCount        altCount        totalCount      lowMAPQDepth    lowBaseQDepth   rawDepth       otherBases      improperPairs   pval    q_val
#1       chr1    10583   chr1:10583      G       A       114     1       115     0       0       116     0       1 
#1       chr1    632366  chr1:567746     T       C       94      0       94      0       0       95      0       1       1.00974195868289e-28    3.15083919755244e-24    29      18      1.87173253230728e-10    1
#my $id = "$info1[0]:$info1[1]:$info1[2]:$info1[3]:$info1[4]";
my $id = "$info1[0]:$info1[1]:$info1[3]:$info1[4]";

$hash{$cell}->{$id}=1;
}
}
for my $cell (keys %hash){
#my $output = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/out-cluster-ATAC_celltype_$cell"."_var";
my $output =  "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/$cell".".snATAC.merged_var";
open(OUTPUT,">$output");
for my $id (keys %{$hash{$cell}}){
print OUTPUT "$id\n";
}
}

