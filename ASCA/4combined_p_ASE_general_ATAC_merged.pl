#!/usr/bin/perl -w
use strict;
#my $var_list = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/out-cluster-ATAC_celltype_$ARGV[1]"."_var_10000_$ARGV[0]";
my $var_list = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/$ARGV[1]".".snATAC.merged_var_10000_$ARGV[0]";
my %var_list;
open(INPUT,$var_list);
while(my $inputv = <INPUT>){
chomp $inputv;
$var_list{$inputv}=1;
}
my %hash=();
#my $file_list = "/storage/chen/home/jw29/sc_human_retina20/data/snRNA_seq_bam20ppl/WASP_corrected_phased20ppl_snp_list";
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
if($cell eq "$ARGV[1]"){
#my $input = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-ATAC_$line".".filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP";
my $input = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/".$line.".snATAC.merged.filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP";
open(INPUT1,$input); # ||die(print $input);
<INPUT1>;
#print OUTPUT "$line1";
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
shift(@info1);
my $id = "$info1[0]:$info1[1]:$info1[3]:$info1[4]";

if(defined $var_list{$id}){
$hash{$cell}->{$id}->{pval} .="$info1[-2]\t";
$hash{$cell}->{$id}->{pval1} .="$info1[-1]\t";
$hash{$cell}->{$id}->{weight} .=sqrt($info1[5]+$info1[6]+$info1[13]+$info1[14])."\t";

#$hash{$cell}->{$id}->{weight} .=sqrt($info1[5]+$info1[6]+$info1[15]+$info1[16])."\t";
$hash{$cell}->{$id}->{count}++;
}
}
}
}

for my $cell (keys %hash){
#my $output = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-ATAC_cell_$cell".".filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP_10000_$ARGV[0]";
my $output = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/".$cell.".snATAC.merged.filter.keep.merged.sorted.snp.count_PASS_DNAread10_WGS_ratio_binP_10000_$ARGV[0]";
open(OUTPUT,">$output");
for my $id (keys %{$hash{$cell}}){
   my $tmp_p = "sc_human_retina/scripts/ASE/combin_tmp_p_ATAC"."_$ARGV[0]"."_$ARGV[1]"."_merged";
   open(OUTPUT1,">$tmp_p");
  print OUTPUT1 "$hash{$cell}->{$id}->{pval}\n";

   my $tmp_p1 = "sc_human_retina/scripts/ASE/combin_tmp_p1_ATAC"."_$ARGV[0]"."_$ARGV[1]"."_merged";
   open(OUTPUT2,">$tmp_p1");
  print OUTPUT2 "$hash{$cell}->{$id}->{pval1}\n";

   my $tmp_w = "sc_human_retina/scripts/ASE/combin_tmp_weigth_ATAC"."_$ARGV[0]"."_$ARGV[1]"."_merged";
   open(OUTPUT3,">$tmp_w");
  print OUTPUT3 "$hash{$cell}->{$id}->{weight}\n";
#exit;
   my ($com_p,$com_z) = combine_p($tmp_p,$tmp_w);
   my ($com_p1,$com_z1) = combine_p($tmp_p1,$tmp_w);

# print OUTPUT "$id\t$com_p\n";
 my $p=$hash{$cell}->{$id}->{pval};
 my $p1=$hash{$cell}->{$id}->{pval1};

 my $w = $hash{$cell}->{$id}->{weight};
 $p =~ s/\t/\;/g;
 $p1 =~ s/\t/\;/g;
 $w =~ s/\t/\;/g;
 print OUTPUT "$id\t$hash{$cell}->{$id}->{count}\t$p\t$w\t$com_p\t$com_z\tgreater\n";
 print OUTPUT "$id\t$hash{$cell}->{$id}->{count}\t$p1\t$w\t$com_p1\t$com_z1\tless\n";


}
}



sub combine_p{
my ($p_file,$w_file) = (@_);
my $R_script = "/storage/chen/home/jw29/sc_human_retina/scripts/ASE/metap_test_ATAC_$ARGV[0]"."_$ARGV[1]"."_merged.R";
my $R_out = "/storage/chen/home/jw29/sc_human_retina/scripts/ASE/metap_test_ATAC_$ARGV[0]"."_$ARGV[1]"."_merged.out";

my $R_comd = "/storage/chen/home/jw29/software/R-3.5.0/bin/R";
open (OUTFILE, ">$R_script");
print OUTFILE  "library(metap)\n";
print OUTFILE  "p_val = read.table(\"$p_file\")\n";
print OUTFILE  "w_val = read.table(\"$w_file\")\n";
print OUTFILE  "sum=sumz(t(p_val),weights=t(w_val))\n"; 
#print OUTFILE  "comZ=sumz(t(p_val),weights=t(w_val))\$z\n"; 
print OUTFILE  "com=cbind(sum\$p,sum\$z)\n";
print OUTFILE "write.table(com,file=\"$R_out\",sep=\"\\t\",quote=F)\n";
`$R_comd <$R_script --no-save`;
open(INPUTR,"$R_out");
<INPUTR>;
my $lineR = <INPUTR>;
chomp $lineR;
my @info = split(/\t/,$lineR);
#my @infoR = split(/\s+/,$lineR);
`rm $R_out`;

return($info[1],$info[2]);
#return($infoR[3],$infoR[4],$infoR[5]);
}
