#!/usr/bin/perl -w
use strict;
#my @sample_list=("D028_13", "D027_13", "D026_13", "D021_13", "D019_13", "D018_13", "D017_13", "D013_13", "D009_13","D005_13", "D030_13", "19_D019", "19_D011", "19_D010","19_D009", "19_D008",  "19_D007", "19_D006","19_D005","19_D003","D19D015","D19D016","D19D013","D19D014","19_D020","19_D017");
my @sample_list=("D028_13", "D027_13", "D026_13", "D021_13", "D019_13", "D018_13", "D017_13", "D013_13", "D009_13","D005_13", "D030_13", "19_D019", "19_D011", "19_D010","19_D009", "19_D008",  "19_D007", "19_D006","19_D005","19_D003");



#my @ct_name=("BC","Astro","Rod","MG","Cone","RGC","HC","AC","Mic");
#my @ct_name=("BC","Astro","Rod","MG","Cone","RGC","HC","AC","Mic");
#my @ct_name=("ONBC","OFFBC","Astro","Rod","MG","Cone","RGC","HC","AC");
my @ct_name=("BC","Astro","Rod","MG","Cone","RGC","HC","AC");

my %people;
#for my $sam (@sample_list){
for(my $s=0;$s<=$#sample_list;$s++){
my $file;
my $sam = $sample_list[$s];
if($s<=19){
#$file = "sc_human_retina20/data/snRNA_seq/$sam"."_snRNA_tmp_ave_new";
###########!$file = "/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq/".$sam."_macula_fovea_snRNA_pool_ave_new";
$file = "/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq/".$sam."_macula_fovea_snRNA_pool_ave_new_afterSoupX_new";

}
#elsif(($s>=20)&&($s<=23)){
#$file = "sc_human_retina20/data/snRNA_seq/$sam"."_macula_fovea_snRNA_cpm_ave_per_cell"; 
#}else{
#$file = "/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq/$sam"."_11182019_5sam_ave_cpm_per_cell";
#}
open(INPUT,$file);
my $line = <INPUT>;
chomp $line;
my @ct = split(/\s+/,$line);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
for(my $i=1; $i<=$#info;$i++){
$people{$ct[$i-1]}->{$info[0]}->{$sam} = $info[$i];
}
}
}

for my $key (keys %people){
#my $output = "sc_human_retina20/data/snRNA_seq/$key"."_snRNA_tmp_ave_macular";
#my $output = "sc_human_retina20/data/snRNA_seq/$key"."_snRNA_tmp_ave_macular_24pp";
#########my $output = "sc_human_retina20/data/snRNA_seq/$key"."_snRNA_pool_ave_macular_20pp";
my $output = "sc_human_retina20/data/snRNA_seq/$key"."_snRNA_pool_ave_macular_20pp_afterSoupX_new";

open(OUTPUT,">$output");
my $header = join("\t",@sample_list);
print OUTPUT "$header\n";
for my $gene (sort keys %{$people{$key}}){
print OUTPUT "$gene\t";
    for(my $i=0; $i<=$#sample_list-1; $i++){
        if(!(defined $people{$key}->{$gene}->{$sample_list[$i]})){
                  $people{$key}->{$gene}->{$sample_list[$i]}=0;
        }
     print OUTPUT "$people{$key}->{$gene}->{$sample_list[$i]}\t";
    }

  if(!(defined $people{$key}->{$gene}->{$sample_list[$#sample_list]})){
                  $people{$key}->{$gene}->{$sample_list[$#sample_list]}=0;
   }
     print OUTPUT "$people{$key}->{$gene}->{$sample_list[$#sample_list]}\n";
}
}
