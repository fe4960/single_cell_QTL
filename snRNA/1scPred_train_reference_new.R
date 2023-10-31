library(scPred)
library(Seurat)
library(magrittr)

meta.data=read.csv("sc_human_retina/data/single_cell/scPred/Human_retina_combined_all_meta.csv",header=T)
counts=read.csv("sc_human_retina/data/single_cell/scPred/Human_retina_combined_all_expression_matrix.csv.gz",header=T,row.names=1)
cellname=gsub("\\.","-",colnames(counts))
colnames(counts)=cellname
reference <- CreateSeuratObject(counts = counts, meta.data=meta.data)
Idents(reference)=meta.data$Cluster
table(Idents(reference))

#rds=RenameIdents(object=reference,  "BC_0" = "OFFBC","BC_3" = "OFFBC","BC_5" = "OFFBC","BC_7" = "OFFBC","BC_8" = "OFFBC", "BC_12" = "OFFBC", "BC_1" = "ONBC","BC_2" = "ONBC","BC_4" = "ONBC","BC_6" = "ONBC","BC_9" = "ONBC","BC_10" = "ONBC","BC_11" = "ONBC",   "AC_8"="AC",   "AC_7"="AC",  "AC_22"="AC",  "AC_30"="AC",  "AC_13"="AC",   "AC_2"="AC",   "AC_0"="AC",  "AC_1"="AC",  "AC_26"="AC", "AC_12"="AC",  "AC_4"="AC",   "AC_9"="AC",  "AC_16"="AC",  "AC_15"="AC",  "AC_21"="AC",   "AC_5"="AC",  "AC_18"="AC",  "AC_19"="AC",   "AC_3"="AC",  "AC_29"="AC", "AC_37"="AC",  "AC_20"="AC",   "AC_6"="AC",  "AC_28"="AC",  "AC_14"="AC",  "AC_25"="AC",  "AC_17"="AC",  "AC_31"="AC",  "AC_10"="AC",  "AC_11"="AC",  "AC_36"="AC", "AC_24"="AC",  "AC_23"="AC",  "AC_38"="AC",  "AC_32"="AC",  "AC_27"="AC",  "AC_35"="AC",  "AC_34"="AC",  "AC_33"="AC",   "NN_0"="MG",   "NN_1"="Astro",   "NN_2"="NN", "NN_3"="NN",   "HC_2"="HC",   "HC_0"="HC",   "HC_1"="HC",  "RGC_6"="RGC",  "RGC_5"="RGC",  "RGC_4"="RGC",  "RGC_2"="RGC",  "RGC_0"="RGC",  "RGC_1"="RGC",  "RGC_3"="RGC","Rod_0"="Rod", "Cone_0"="Cone", "Cone_1"="Cone")

reference=RenameIdents(object=reference, "Astrocytes"="Astro", "BB+GB*"="BC", "DB1"="BC", "DB2"="BC", "DB3a"="BC", "DB3b"="BC", "DB4"="BC", "DB5*"="BC", "DB6"="BC", "Endothelium"="Endo", "FMB"="BC", "Gaba1"="AC","Gaba10"="AC","Gaba11"="AC","Gaba12"="AC", "Gaba13"="AC", "Gaba14"="AC", "Gaba15"="AC", "Gaba16"="AC", "Gaba17"="AC", "Gaba2"="AC", "Gaba3"="AC", "Gaba4"="AC", "Gaba5"="AC", "Gaba6"="AC", "Gaba7"="AC", "Gaba8"="AC", "Gaba9"="AC", "Gly1"="AC", "Gly2"="AC", "Gly3"="AC", "Gly4"="AC", "Gly5"="AC", "Gly6"="AC", "Gly7"="AC", "Gly8"="AC", "H1"="HC", "H2"="HC", "IMB"="BC", "MG_OFF"="RGC", "MG_ON"="RGC", "MicroGlia"="Mic", "Muller"="MG", "OFFx"="BC", "PG_OFF"="RGC", "PG_ON"="RGC", "RB1"="BC", "RGC10"="RGC", "RGC11"="RGC", "RGC12"="RGC", "RGC5"="RGC", "RGC6"="RGC", "RGC7"="RGC", "RGC8"="RGC", "RGC9"="RGC", "Rods"="Rod", "mlCones"="Cone","sCones"="Cone")

reference$cell_type=Idents(reference)

#normalize reference
reference <- reference %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

pdf("sc_human_retina/data/single_cell/scPred/Human_retina_combined_reference_Sanes_new.pdf")
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)
dev.off()

reference <- getFeatureSpace(reference, "cell_type")

######train the classifiers
reference <- trainModel(reference)
#source("R/x86_64-pc-linux-gnu-library/4.0/scPred/generics.R")
get_probabilities(reference) %>% head()

get_scpred(reference)

pdf("sc_human_retina/data/single_cell/scPred/Human_retina_combined_reference_pred_prob_new.pdf")
plot_probabilities(reference)
dev.off()

saveRDS(reference, "sc_human_retina/data/single_cell/scPred/Human_retina_combined_reference_new.rds")
#reference <- trainModel(reference, model = "mda", reclassify = c("cMono", "ncMono"))

#get_scpred(reference)

#plot_probabilities(reference)

#query <- NormalizeData(query)

#query <- scPredict(query, reference)


