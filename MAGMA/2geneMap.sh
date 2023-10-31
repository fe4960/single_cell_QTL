#!/bin/sh
export LD_LIBRARY_PATH=/storage/chen/home/jw29/software/nlopt-2.4.2/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/storage/chen/home/jw29/software/nlopt-2.4.2/lib:$LIBRARY_PATH
export PATH=/storage/chen/home/jw29/software/nlopt-2.4.2/bin:$PATH

export LD_LIBRARY_PATH=/storage/chen/Software/gcc-7.3.0_build/lib64:$LD_LIBRARY_PATH
export PATH=/storage/chen/home/jw29/software:/storage/chen/Software/gcc-7.3.0_build/bin:$PATH
#software/R-4.0.0/bin/Rscript --vanilla $R_command $all_gene
#software/R-4.1.0/bin/R < human_meta/scripts/GWAS/2geneMap.R --no-save > human_meta/scripts/GWAS/2geneMap.out 2> human_meta/scripts/GWAS/2geneMap.err & 
#file="2geneMap_list"
#file="2geneMap_list1"
#for file in 2ctd_list_4_0 2ctd_list_4_1 2ctd_list_4_2 2ctd_list_4_3
#for file in 2ctd_list_4_4 2ctd_list_4_5
file=$1

software/R-4.1.0/bin/Rscript --vanilla human_meta/scripts/GWAS/2geneMap.R $file > human_meta/scripts/GWAS/2geneMap_${file}.out 2> human_meta/scripts/GWAS/2geneMap_${file}.err 

#file="2geneMap_list1"

