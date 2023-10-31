#!/usr/bin/perl -w
use strict;
#my $file_list = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam-snp_count_list_48_0";
my $file_list = "$ARGV[0]";
open(INPUT,$file_list);
while(my $line = <INPUT>){
chomp $line;
#19_D003_AC.snp.count
$line =~ s/\.snp\.count//g;
my @id = split(/\_/,$line);
#my $id = $id[0]."_".$id[1];
my $id = $id[1]."_".$id[2];
my $vcf ="sc_human_retina20/data/WGS20/$id".".GATK.HaplotypeCaller.mark.snp.vcf";
#D027_13.GATK.HaplotypeCaller.mark.snp.vcf
open(INPUT0,$vcf);
my %hash=();
while(my $line0 = <INPUT0>){
if($line0 =~ /^#/){
next;
}
chomp $line0;
my @info0 = split(/\s+/,$line0);
my @gt = split(/\:/,$info0[9]);
my @read = split(/\,/,$gt[1]);
my $read_sum = $read[0] + $read[1];

if(($info0[6] eq "PASS")&&($read_sum>=10)&&($read[0]>0)&&($read[1]>0)){
my $id = "$info0[0]:$info0[1]:$info0[3]:$info0[4]";
$hash{$id}->{filter} = $info0[6];
$hash{$id}->{total} = $read_sum;
$hash{$id}->{ref} = $read[0];
$hash{$id}->{alt} = $read[1];

}
}

#my $input = "/storage/chen/home/jw29/sc_human_retina20/data/snATAC_seq_bam20ppl/WASP_corrected_phased20ppl/out-cluster-ATAC_$line".".filter.keep.merged.sorted.snp.count";
my $input = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular/WASP_corrected_phased20ppl/$id[0]_$id.snATAC.merged.filter.keep.merged.sorted.snp.count";
my $output = $input."_PASS_DNAread10_WGS_ratio";
open(OUTPUT,">$output");
open(INPUT1,$input); #||die(print $input);
my $line1=<INPUT1>;
#print OUTPUT "$line1";
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\s+/,$line1);
#shift(@info1);
#contig  position        variantID       refAllele       altAllele       refCount        altCount        totalCount      lowMAPQDepth    lowBaseQDepth   rawDepth
#        otherBases      improperPairs
#        1       14907   .       A       G       0       2       2       0       0       2       0       0
my $id = "$info1[0]:$info1[1]:$info1[3]:$info1[4]";
$id =~ s/chr//g;
my $RNA_read = $info1[5]+$info1[6];
if(($hash{$id}->{total} >=10) && ($hash{$id}->{filter} eq "PASS")&&($RNA_read>=10)){
print OUTPUT "$line1\t$hash{$id}->{ref}\t$hash{$id}->{alt}\n";
}
}
}
