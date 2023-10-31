#!/usr/bin/perl -w
my $control_file = "/storage/chenlab/Users/junwang/GWAS_analysis/scripts/MAGMA/file_list_new_control_MAGMA";
my %control;
open(INPUT,$control_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$control{$info[0]} = $info[1];
}
my $output = $control_file."_summary_list_new";

open(OUTPUT,">$output");
print OUTPUT "P-val\tCell_type\tTrait\n";
my @list = ("3MAGMA_list","3MAGMA_list1","3MAGMA_list2");
for my $list (@list){
my $dir = "/storage/chenlab/Users/junwang/GWAS_analysis/scripts/MAGMA/$list";
open(INPUT,$dir);
while(my $line = <INPUT>){
chomp $line;
#######MAGMA output file name
my $file = "$line/MAGMA_result1"; #ad
open(INPUT1,$file)|| die ("$file");
<INPUT1>;
while(my $line1 = <INPUT1>){
#GWAS    Celltype        OBS_GENES       BETA    BETA_STD        SE      P       level   Method  GCOV_FILE       CONTROL CONTROL_label   log10p  genesOutCOND    EnrichmentMode  FDR     Celltype_id
#1       Fritsche-26691988.txt_reform_MungeSumstats.txt.35UP.10DOWN       E n d o        15138   0.0019514       0.024517        0.00065923      0.0015406
#       1       MAGMA   Fritsche-26691988.txt_reform_MungeSumstats.txt.35UP.10DOWN.level1.retina_linear.gsa.out BASELINE        BASELINE        -2.81231010647482       NA      Linear  0.030812        Endo
chomp $line1;
my @info1 = split(/\t/,$line1);
$info1[1] =~ s/_MungeSumstats.txt_reform.35UP.10DOWN//g;
$info1[1] =~ s/_MungeSumstats.txt1_reform.35UP.10DOWN//g;
$info1[1] =~ s/_MungeSumstats_rmFRQ.txt_reform.35UP.10DOWN//g;
$info1[1] =~ s/_MungeSumstats.txt_rmFQS_reform.35UP.10DOWN//g;
if($info1[-3] eq "Linear"){
print OUTPUT "$info1[7]\t$info1[-1]\t$control{$info1[1]}\n";
}
}
}
}
#/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/MAGMA_Files/Fritsche-26691988.txt_reform_MungeSumstats.txt.35UP.10DOWN/MAGMA_result
 

