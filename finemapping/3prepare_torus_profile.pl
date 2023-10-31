#!/usr/bin/perl 
use strict;
use Math::CDF qw(qnorm);
#my $file = "sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz";
my %G1000_var;
for(my $i=1; $i<=22; $i++){
#my $bim = "sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_chr"."$i"."_1000G.bim"; 
my $bim = $ARGV[0]."_chr$i"."_1000G.bim";
#9       rs4741076       25.284641       10978031        T       G
open(INPUTb,$bim);
while(my $lineb = <INPUTb>){
chomp $lineb;
my @info = split(/\s+/,$lineb);
my $var = "chr$info[0]:$info[3]:$info[4]:$info[5]";
$G1000_var{$var}=$info[1];
}
}
#my $peak_info ="sc_human_retina/data/cell_atlas/allpeak_info"; 

#my @cell=("ONBC", "OFFBC", "RGC", "Cone", "Rod", "MG", "HC","AC","Astro");
my @cell = ("BC", "RGC", "Cone", "Rod", "MG", "HC","AC","Astro");
my %peak_all=();
my %peak_var=();
#for my $cell (@cell){
#my $peak = "sc_human_retina20/data/WGS20/HaplotypeCaller_vcf/cohort_calling/WGS20more_chrAll.GATK.HaplotypeCaller.filtered_rmD017D020_RS_filtered.vcf_wH_flt_lcr_dosage_plink_geno_site_inPeak_archr_new_".$cell."_macular_lobe_20ppl_bed";
my $peak = $ARGV[1]; #####variant peak annotation
open(INPUT,$peak);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $peak = "chr$info[-4]:$info[-3]-$info[-2]";
$peak_all{$peak}=1;
my $pos = "$info[0]:$info[1]";
#print "$pos\n";
#exit;
$peak_var{$pos}=$peak;
}
#}

my %DAR;
for my $cell (@cell){
#my $file = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/$cell"."_DAR_07-23-2021_flt_bin_new_1fpm";
#my $file = "sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/$cell"."_DAR_peak";
my $file = "/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_peak/$cell"."_DAR_06-22-2022_flt_bin_new_2fpkm";
open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\t/,$line);
#my $peak = "$info[1]:$info[3]-$info[4]";
my $peak = "$info[0]:$info[1]-$info[2]";
$DAR{$peak}->{$cell}=1;
}
close(INPUT);
}
#my $peak_all_anno ="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/all_peak_anno_bin_new_1fpm";
my $peak_all_anno = "/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_peak/lobe_macular_macs/all_peak_anno_bin_new_2fpkm";
#my $peak_all_anno ="/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/all_peak_anno";
my %peak_gene;
my %peak_anno;
open(INPUT,$peak_all_anno);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\t/,$line);

my $var_pos = "$info[1]:$info[2]-$info[3]";
$info[1] =~ s/chr//g;
#my $var_pos = "$info[1]:$info[2]";
$peak_gene{$var_pos} = $info[16];
if(($info[6] =~ /Intron/)||($info[6] =~ /Exon/)||($info[6] =~ /Promoter/)||($info[6] =~ /Downstream/)){
my @info1 = split(/\s+/,$info[6]);
$peak_anno{$var_pos} = $info1[0];
}else{
my @info1 = split(/\s+/,$info[6]);
my $tmp = join("_",@info1);
$peak_anno{$var_pos} = $tmp;
}
}


my %CRE;
#my $coAC = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/coAC_cor0.5_resolution1_bin_promoter";

#my $coAC = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/coAC_cor0.5_resolution1_promoter";
#open(INPUT,$coAC);
#while(my $line = <INPUT>){
#chomp $line;
#my @info = split(/\t/,$line);
#my $peak = "$info[29]:$info[30]-$info[31]";
#my $gene = $info[27];
#my $link = "$peak"."_"."$gene";
#$CRE{$peak}->{$gene}=1;
#}

#my $g2p = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/gene2peak_cor0.5_fdr0.01_resolution1_bin_rmPromoter";

#my $g2p = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/gene2peak_cor0.5_fdr0.01_resolution1_rmPromoter";
#open(INPUT,$g2p);
#my %var_peak_link_gene;
#<INPUT>;
#while(my $line = <INPUT>){
#chomp $line;
#my @info = split(/\t/,$line);
#my $gene = $info[22];
#my $gene = $info[-2];
#my $peak = "$info[7]:$info[8]-$info[9]";
#my $link = $peak."_".$gene;
#$CRE{$peak}->{$gene}=1;
#}

my %CRE_gene;
my $CRE_list = "sc_human_retina/data/snATAC_snRNA/LCRE_list_2fpkm";

#my $CRE_list = "sc_human_retina/data/snATAC_snRNA/LCRE_list_1fpm";
open(INPUT,$CRE_list);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $peak = "chr$info[1]:$info[2]-$info[3]";
my $link = $peak."_".$info[0];
$CRE{$peak}->{$info[0]}=1;
$CRE_gene{$info[0]}->{$peak}=1;
}


#my $peak_all = "sc_human_retina/data/data_eQTL_cpm10/fastqtl_250kb_maf01_nsa4_10cpm_cellPeak_20ppl_plot/eQTL_sc_p0.01_bulk_matrix_FDR1_chipseeker";
my $peak_all = $ARGV[2]; ####variant annotation
my %peak;
my %gene;
open(INPUT,$peak_all);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\t/,$line);
$info[1] =~ s/chr//g;
my $var_pos = "$info[1]:$info[2]";
$gene{$var_pos} = $info[16];
#print "$var_pos\n";
#exit;
if(($info[6] =~ /Intron/)||($info[6] =~ /Exon/)||($info[6] =~ /Promoter/)||($info[6] =~ /Downstream/)){
my @info1 = split(/\s+/,$info[6]);
$peak{$var_pos} = $info1[0];
}else{
my @info1 = split(/\s+/,$info[6]);
my $tmp = join("_",@info1);
$peak{$var_pos} = $tmp;
}
}


#my $peak_info = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed.mainChr_rmBlacklist_rmY";
#my $file = "sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8";
my $file = "$ARGV[0]_format.vcf";
my $ldetect = "sc_human_retina/data/GWAS_new/ldetect/EUR/fourier_ls-all.bed.gz"; 
#my $output1=$file.".1000g_flt_20ppl.zscore";
#my $output2=$file.".1000g_flt_20ppl.anno";
my $output1=$file.".1fpm.1000g_flt_20ppl.zscore";
my $output2=$file.".1fpm.1000g_flt_20ppl.anno";

my $output1_gz = "$output1".".gz";
my $output2_gz = "$output2".".gz";
`rm $output1_gz`;
`rm $output2_gz`;
open(OUTPUT1,">$output1");
open(OUTPUT2,">$output2");
print OUTPUT2 "SNP\tannot_d\n";
#chr:6:29027255  a       g       0.4933  0.0068  232207.00       5.977   2.271e-09       ++++    0.0     1.056   3       0.7876
open(INPUT,$file);
<INPUT>;
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my $sign_dir;
if($info[5] > 0){
$sign_dir =1;
}elsif($info[5]<0){
$sign_dir=-1;
}

my $zscore = $sign_dir*abs(qnorm($info[6]/2));

#my $zscore = $info[5];
#my @coor = split(/\:/,$info[0]);
my $chr = "chr$info[0]";
my $var = "$chr:$info[1]:".uc($info[4]).":".uc($info[3]);
my $oppo_var = "$chr:$info[1]:".uc($info[3]).":".uc($info[4]);

#print "$var\n$oppo_var\t$zscore\n";

my $rs;
if(defined $G1000_var{$var}){
$rs=$G1000_var{$var};
}elsif(defined $G1000_var{$oppo_var}){
$zscore = -$zscore;
$rs=$G1000_var{$oppo_var};
}else{
next;
}
#print "$rs\n";
#my $tmp = "sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/20ppl_tmp";
my $tmp = $ARGV[0]."_tmp";
`rm $tmp`;
`less sc_human_retina/data/GWAS_new/ldetect/EUR/fourier_ls-all.bed.gz | awk '{if((\$1==\"$chr\")&&(\$2<=$info[1])&&(\$3>=$info[1])){print \$1\":\"\$2\"-\"\$3}}' > $tmp`;
open(INPUT1,"$tmp");
my $ld = <INPUT1>;
#print "$ld";
#exit;
chomp $ld;
print OUTPUT1 "$chr:$info[1]:$rs:$info[4]:$info[3]\t$ld\t$zscore\n";
#`rm $tmp`;
#######
#assign gene name
#######
my $var_pos = $info[0].":".$info[1];
my $cate=0;
my $assign=0;

#print "$peak{$var_pos}\n";
#exit;
if($peak{$var_pos} =~ /Exon/){
$cate=4; ####label gene Exon
$assign=1;
}

if(($assign==0)&&(($peak{$var_pos} =~ /3'_UTR/ )||( $peak{$var_pos} =~ /5'_UTR/  ))){
$cate=4; #####label gene UTR
$assign=1;
}
if(($assign==0)&&($peak{$var_pos} =~ /Promoter/)){
$cate=3; #####label gene UTR
$assign=1;
}


if(defined $peak_var{$var_pos}){
my $peak = $peak_var{$var_pos};
#print "$peak\n";
#exit;
#####if(($assign==0)&&($peak_anno{$peak} =~ /Promoter/)){

#if(($assign==0)&&($peak_anno{$peak} =~ /Promoter_\(\<\=1kb\)/)){
####$cate=3; #########label promoer
#####$assign=1;
########}

if(($assign==0)&& (defined $CRE{$peak})){
$cate=2; #####label CRE
$assign=1;
}

if($assign==0){
$cate=1; #####label peak
$assign=1;
}
}

if($assign==0){
$cate=0;
$assign=1;
}

print OUTPUT2 "$chr:$info[1]:$rs:$info[4]:$info[3]\t$cate\n";




}

`bgzip $output1`;
`bgzip $output2`;

