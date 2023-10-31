library(qvalue)
files=c("MAGMA")
for(file in files){
value=read.table(paste0("~/Documents/human_meta/GWAS_file_list_new_control_",file,"_summary_list_new"),header=T)


value$Trait=gsub("POAG_GERA_UKB2018","POAG_GERA",value$Trait)
value$Trait=gsub("AMD_NG2016","AMD",value$Trait)
value$Trait=gsub("Khor_PACG2016","PACG",value$Trait)
value$Trait=gsub("Myopia2018NG","Myopia",value$Trait)
value$Trait=gsub("CA_DA_caucasians_2017","POAG_CA",value$Trait)
value$Trait=gsub("DA_caucasians_2017","POAG_DA",value$Trait)
value$Trait=gsub("IOP_caucasians_2017","POAG_IOP",value$Trait)
value$Trait=gsub("VCDR_caucasians_2017","POAG_VCDR",value$Trait)
value$Trait=gsub("UKBBc_IOP_2018Khawaja","IOP",value$Trait)
value$Trait=gsub("POAG_Gharahkhani2021","POAG",value$Trait)
value$Trait=gsub("Refracive_Error_Hysi2020","Refracive_Error",value$Trait)
value$Trait=gsub("Bone_mineral_density","Bone_mineral_density",value$Trait)
value$Trait=gsub("WBC","WBC",value$Trait)
value$Trait=gsub("UKBB_BMI","BMI",value$Trait)
value$Cell_type=gsub("RBC","BCR",value$Cell_type)
value$Cell_type=gsub("Astrocyte","Astro",value$Cell_type)
value$Cell_type=gsub("Microglia","Micro",value$Cell_type)
value=value[value$Cell_type!="Endo"&value$Cell_type!="Micro"&value$Cell_type!="RPE",]

value=value[value$Trait%in%c("Disorders_of_choroid_and_retina","Retinal_detachments_and_breaks","Retinal_problem","AMD","Refracive_Error","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","Bone_mineral_density","Diabetic_retinopathy","ONL_thickness","IS_thickness","OS_thickness"),]


value$p.adj=p.adjust(value$P.val,method="BH")
value$qval=qvalue(value$P.val)$qvalue

p <- ggplot(value, aes(factor(Cell_type,levels=c("Rod","Cone","BC","HC","AC","RGC","MG","Astro")),factor(Trait,levels=c("Disorders_of_choroid_and_retina","Retinal_detachments_and_breaks","Retinal_problem","AMD","Refracive_Error","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","Diabetic_retinopathy","ONL_thickness","IS_thickness","OS_thickness","Bone_mineral_density"))))

p_dar=p + geom_point(aes(size=-log(value$P.val,base=10)),color="red")+theme(
    axis.line=element_line(colour="black")
    ,panel.background = element_rect(fill="transparent") # bg of the panel
    , plot.background = element_blank() # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    #, legend.position = "none"
    , text=element_text(size=20)
   # axis.text.x = element_text(size=10)  # get rid of legend bg
  )+xlab("")+ylab("")+ggtitle("snRNA-seq Gene expression")+
  guides(color = guide_colorbar(order = 1),fill = guide_legend(order = 1))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(1,15))
 pdf(paste0("~/Documents/human_meta/file_list_new_control_",file,"_summary_list_new.pdf"),width=10.5)
 print(p_dar)
 dev.off()
 #write.table(value,file=paste0("~/Documents/sc_human_retina/analysis_2fpkm/file_list_new_control_",file,"_summary_list_new_padj"),quote=F,sep="\t")
 write.table(value,file=paste0("~/Documents/human_meta/file_list_new_control_",file,"_summary_list_new_padj1"),quote=F,sep="\t")
 }

