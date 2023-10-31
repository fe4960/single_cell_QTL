#!/usr/bin/perl -w
use strict;
my @pp = ("19_D003", "19_D005", "19_D006", "19_D007", "19_D008", "19_D009", "19_D010", "19_D011", "19_D019", "D005_13", "D009_13", "D013_13", "D017_13", "D018_13", "D019_13", "D021_13", "D026_13", "D027_13", "D028_13", "D030_13");

my @cell = ("Rod","Cone","RGC","MG","HC","BC","AC");

my $file = "sc_human_retina/data/QC/WGS20more_autosome_flt_hapmap.mds";

my %mds;
open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$mds{$info[1]} = $info[4];
}

for my $cell (@cell){
my $exp = "sc_human_retina20/data/snRNA_seq/normalized_avg5cmp_".$cell."_macular_20pp_mean5cpm_3peer_factor_afterSoupX_new";
my %cov=();
open(INPUT,$exp);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$info[0] =~ s/X//g;
for(my $i=1;$i<=3; $i++){
$cov{$info[0]}->{$i} = $info[$i];
}
}

my $output = "sc_human_retina20/data/snRNA_seq/covariate_20pp_".$cell."_3peer_1mds_avg5cpm_afterSoupX_new";
open(OUTPUT,">$output");

print OUTPUT "id";

for(my $i=0; $i<=$#pp; $i++){
print OUTPUT "\t$pp[$i]";
}
print OUTPUT "\n";

for(my $i=1; $i<=3; $i++){
print OUTPUT "PC$i";
#}

for(my $j=0; $j<=$#pp; $j++){
print OUTPUT "\t$cov{$pp[$j]}->{$i}";
}

print OUTPUT "\n";
}
print OUTPUT "PC4";

for(my $i=0; $i<=$#pp; $i++){
print OUTPUT "\t$mds{$pp[$i]}";
}
print OUTPUT "\n";
}
