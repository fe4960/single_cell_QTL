#!/usr/bin/perl -w
use strict;
my @pp = ("19_D003", "19_D005", "19_D006", "19_D007", "19_D008", "19_D009", "19_D010", "19_D011", "19_D019", "D005_13", "D009_13", "D013_13", "D017_13", "D018_13", "D019_13", "D021_13", "D026_13", "D027_13", "D028_13", "D030_13");
my @cell = ("Rod","Cone","RGC","MG","HC","BC","AC");
my $gene ="sc_human_retina/data/eQTL/gencode.v19.genes.v7.patched_contigs.gtf_gene_tss_tes"; 
#1	11869	14362	ENSG00000223972.4	DDX11L1
my %gene;
open(INPUT,$gene);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$gene{$info[-1]}->{tss} = $info[1];
$gene{$info[-1]}->{chr} = $info[0];
}
for my $cell (@cell){
my $exp = "sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_".$cell."_pool_macular_20pp_afterSoupX_new";
#my %gene_info;
#1       11869   14362   ENSG00000223972.4       DDX11L1
open(INPUTg,$exp);
my $header = <INPUTg>;
chomp $header;
my %exp=();
$header =~ s/X//g;
my @header = split(/\s+/,$header);
while(my $lineg = <INPUTg>){
my @infog = split(/\s+/,$lineg);
for(my $i=1;$i<=$#infog; $i++){
$exp{$infog[0]}->{$header[$i-1]} = $infog[$i];
}
}
my $output = $exp."_bed";
open(OUTPUT,">$output");
print OUTPUT "#Chr\tstart\tend\tTargetID";
for(my $i=0; $i<=$#pp;$i++){
print OUTPUT "\t$pp[$i]";
}
print OUTPUT "\n";
for my $gene ( keys %exp){
if(defined $gene{$gene}){
print OUTPUT "$gene{$gene}->{chr}\t",$gene{$gene}->{tss}-1,"\t",$gene{$gene}->{tss},"\t$gene";
for(my $i=0; $i<=$#pp;$i++){
print OUTPUT "\t$exp{$gene}->{$pp[$i]}";
}
print OUTPUT "\n";
}
}
my $output_sort = $output."_sort";
my $output_sort_gz  =$output_sort.".gz";
`(head -n 1 $output && tail -n+2 $output | sort -k 1,1 -k 2n,2n -k 3n,3n) > $output_sort`;
`rm $output_sort_gz*`;
`bgzip $output_sort`;
`tabix -p bed $output_sort_gz`;
} 


