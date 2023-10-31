#!/bin/sh
command=sc_human_retina/scripts/GWAS_new/finemapping/5pip_plot_general
file=$1
software/R-4.0.0/bin/Rscript --vanilla  ${command}.R  $file
