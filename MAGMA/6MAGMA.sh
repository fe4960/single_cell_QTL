#!/bin/sh

export LD_LIBRARY_PATH=/storage/chen/home/jw29/software/nlopt-2.4.2/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/storage/chen/home/jw29/software/nlopt-2.4.2/lib:$LIBRARY_PATH
export PATH=/storage/chen/home/jw29/software/nlopt-2.4.2/bin:$PATH

export LD_LIBRARY_PATH=/storage/chen/Software/gcc-7.3.0_build/lib64:$LD_LIBRARY_PATH
export PATH=/storage/chen/home/jw29/software:/storage/chen/Software/gcc-7.3.0_build/bin:$PATH

file="all"
#file=$1
software/R-4.1.0/bin/Rscript --no-vanilla human_meta/scripts/GWAS/6MAGMA.R $file > human_meta/scripts/GWAS/6MAGMA_${file}.out 2> human_meta/scripts/GWAS/6MAGMA_${file}.err
