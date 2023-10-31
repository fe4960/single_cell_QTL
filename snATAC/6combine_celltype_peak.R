library(GenomicRanges)
setwd("sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/")
peaks.names=c("Rod_macs3_peaks.narrowPeak","MG_macs3_peaks.narrowPeak","Cone_macs3_peaks.narrowPeak","ONBC_macs3_peaks.narrowPeak","OFFBC_macs3_peaks.narrowPeak","RGC_macs3_peaks.narrowPeak","HC_macs3_peaks.narrowPeak","AC_macs3_peaks.narrowPeak","Astro_macs3_peaks.narrowPeak")
peak.gr.ls = lapply(peaks.names, function(x){
    print(x)
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })
peak.gr = reduce(Reduce(c, peak.gr.ls));
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "/storage/chen/home/jw29/sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/narrowPeaks_macs3.combined.bed",append=FALSE,quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
#write.table(peaks.df,file = "/storage/chen/home/jw29/sc_human_retina/data/ATAC_peak/bulk_peaks.combined.bed",append=FALSE,quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
