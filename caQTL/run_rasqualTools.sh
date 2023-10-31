#!/bin/sh
export RASQUALDIR="/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/rasqual/"
#export CFLAGS="-I/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/INCLUDE -I/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/F2CLIBS -I/storage/chen/home/jw29/software/gsl-2.6/"
#export LDFLAGS="-L/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/ -L/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/F2CLIBS -L/storage/chen/home/jw29/software/gsl-2.6/lib"

export CFLAGS="-I/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/INCLUDE -I/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/F2CLIBS "
export LDFLAGS="-L/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/ -L/stornext/snfs130/ruichen/dmhg/fgi/jwang/software/CLAPACK-3.2.1/F2CLIBS "

#export LD_LIBRARY_PATH="/storage/chen/home/jw29/software/gsl-2.6/lib:$LD_LIBRARY_PATH"
perl /stornext/snfs130/ruichen/dmhg/fgi/jwang/sc_human_retina/ca_eQTL/run_rasqualTools.pl $1 $2
#cd $RASQUALDIR
#cell=$1
#chr=$2
#count_table=$3
#sample_offset=$4
#vcf="/storage/chen/home/jw29/sc_human_retina/data/ca_eQTL/${cell}/1000GP_Phase3_20ppl_chr${chr}.phased.wRef.merge.vcf.correct_ref_new.flt.new.gz"
#file=

#while read line

#IFS='\t' read -ra ARRY <<< "$line"
#region=${my_array[1]}
#geneID=${my_array[2]}
#num_test=${my_array[3]}
#num_feature=${my_array[4]}
#start_pos=${my_array[5]}
#end_pos=${my_array[6]}
#feature=${my_array[7]}

#tabix ${vcf} ${region} | bin/rasqual -y ${count_table} -k ${sample_offset} -n 20 -j $geneID -l $num_test -m $num_feature -s $start_pos -e $end_pos -f $feature 

#done < $file
