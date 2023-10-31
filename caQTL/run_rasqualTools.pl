#!/usr/bin/perl
#`export RASQUALDIR="/storage/chen/home/jw29/software/rasqual/"`;
#`export CFLAGS="-I/storage/chen/home/jw29/software/CLAPACK-3.2.1/INCLUDE -I/storage/chen/home/jw29/software/CLAPACK-3.2.1/F2CLIBS -I/storage/chen/home/jw29/software/gsl-2.6/"`;
#`export LDFLAGS="-L/storage/chen/home/jw29/software/CLAPACK-3.2.1/ -L/storage/chen/home/jw29/software/CLAPACK-3.2.1/F2CLIBS -L/storage/chen/home/jw29/software/gsl-2.6/lib"`;
my $main = "/stornext/snfs130/ruichen/dmhg/fgi/jwang/sc_human_retina/ca_eQTL";
my $cell=$ARGV[0];
my $chr=$ARGV[1];
my $count_table="$main"."/rasqualTool/".$cell."/cellType.expression.bin";
my $sample_offset="$main"."/rasqualTool/".$cell."/cellType.size_factors_gc.bin";
my $cov="$main"."/rasqualTool/".$cell."/cellType.covariate.bin";
my $vcf="$main"."/".$cell."/1000GP_Phase3_20ppl_chr".$chr.".phased.wRef.merge.vcf.correct_ref_new.flt.new.gz";
my $file="$main"."/rasqualTool/cellType.snp_counts.txt";
#1       chr10:100002497-100003133       chr10   *       100002497       100003133       100002297       100003333       0       0
my %hash;
my $cell_peak = "$main/rasqualTool/".$cell."/cellType.expression.txt";
my $output = "$cell_peak"."_$cell"."_$chr";
`rm $output`;

open(INPUT1,$cell_peak);

open(INPUT,$file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);

my $region="$info[2]:$info[6]-$info[7]";
my $geneID=$info[1];

my $num_test=$info[-1];
my $num_feature=$info[-2];
my $start_pos=$info[4];
my $end_pos=$info[5];
my $feature=$info[1];

$hash{$geneID}->{region} = $region;
$hash{$geneID}->{num_test} = $num_test;
$hash{$geneID}->{num_feature} = $num_feature;
$hash{$geneID}->{start_pos} = $start_pos;
$hash{$geneID}->{end_pos} = $end_pos;
#print "$geneID\n";
#exit;
}
my $id=0;
while(my $line1 = <INPUT1>){
$id++;
chomp $line1;
my @info1 = split(/\s+/,$line1);
my @coor = split(/\:|\-/,$info1[0]);
#print "$coor[0]\t$coor[1]\t$coor[2]\n";
my $geneID = $info1[0];
#print "$geneID\n";
#exit;
my $tmp_chr = $coor[0];
$tmp_chr =~ s/chr//g;
if((defined $hash{$geneID})&&($tmp_chr eq $chr)){
#print "enter\n";
my $num_test=$hash{$geneID}->{num_test};
my $num_feature=$hash{$geneID}->{num_feature};
my $start_pos=$hash{$geneID}->{start_pos};
my $end_pos=$hash{$geneID}->{end_pos};
my $feature=$geneID;
my $region = $hash{$geneID}->{region};
#`tabix $vcf $region | /storage/chen/home/jw29/software/rasqual/bin/rasqual -y $count_table -k $sample_offset -n 20 -j $geneID -l $num_test -m $num_feature -s $start_pos -e $end_pos -f $feature -x $cov > $tmp_out`;
`/hgsc_software/tabix/tabix-0.2.6/tabix $vcf $region | /stornext/snfs130/ruichen/dmhg/fgi/jwang/software/rasqual/bin/rasqual -y $count_table -k $sample_offset -n 20 -j $id -l $num_test -m $num_feature -s $start_pos -e $end_pos -f $feature -x $cov >> $output`;

#print "$tmp_out\n";
#exit;
}
}
