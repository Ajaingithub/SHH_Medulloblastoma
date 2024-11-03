library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
shh_snRNA <- readRDS(paste(savedir,"saveRDS_obj/SHH.merged_all_samples.hg38.qc.filtered.rds",sep = ""))

dir.create(paste(savedir,"dimplot",sep = ""), showWarnings = FALSE)

## QC and creating the Seurat Object made a minute change in scRNA_QC so it will read table
split_obj <- SplitObject(shh_snRNA, split.by = "sample")
obj_name <- names(split_obj)
obj_name_2 <- gsub("-","_",obj_name)


### Adding the Run information to it
library(stringr)
sample_batch <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_batch.txt", header = TRUE)

shh_snRNA@meta.data$SampleID <- gsub("_.*.","",shh_snRNA@meta.data$sample_id)
patterns = sample_batch$SampleID
replacements = sample_batch$Run
shh_snRNA@meta.data$Run <- shh_snRNA@meta.data$SampleID
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  shh_snRNA@meta.data$Run <- str_replace_all(shh_snRNA@meta.data$Run, pattern, replacements[i])
}

#### Need to make some changes in the samplename in the orig.ident
all_metadata <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/snRNA_metadata.txt", header = TRUE, sep = "\t")
all_metadata$Sample.ID2 <- gsub("7316-","",all_metadata$Sample.ID)
all_metadata$Tumor.type <- gsub("Primary","P",all_metadata$Tumor.type) %>% gsub("Recurrent","R",.)
all_metadata$Batch..snRNA.seq. <- gsub("Batch","Run",all_metadata$Batch..snRNA.seq.)
all_metadata$orig.ident <- paste(all_metadata$Project, all_metadata$Sample.ID2,"_", 
      all_metadata$Age,"_", all_metadata$Sex,"_", 
      all_metadata$Tumor.type,"_", all_metadata$Batch..snRNA.seq.,sep="")

patterns = all_metadata$Sample.ID
replacements = all_metadata$orig.ident
shh_snRNA@meta.data$orig.ident <- shh_snRNA@meta.data$SampleID
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  shh_snRNA@meta.data$orig.ident <- str_replace_all(shh_snRNA@meta.data$orig.ident, pattern, replacements[i])
}

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
objname = "SHH"
Assay = "RNA"
process = "sctransform"
shh_sctransformed <- sctransform_V2_integration(obj = shh_snRNA,
                                                saveDir = savedir,
                                                ngenes = 4000,
                                                regress = c("nCount_RNA"),
                                                dims = 30,
                                                Assay = Assay,
                                                process = process,
                                                objname = objname,
                                                split_by = "orig.ident",
                                                reference = c(9),
                                                sample_tree = NULL)

saveRDS(shh_sctransformed, paste(savedir,"saveRDS_obj/SHH_sample_integrated_obj.RDS",sep = ""))

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_RNA_intergation.R")
objname = "SHH"
process = "integration"
Assay = "RNA"
shh_integrated_RNA <- RNA_integration(shh_sctransformed, savedir,
                                      dims = 30, RNA_features = c("CD4","CD8A"),
                                      Assay=Assay, process=process, objname=objname,
                                      ncol = 4, ndims = 50)

shh_integrated_RNA_NN <- FindNeighbors(shh_integrated_RNA, dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
source("/diazlab/data3/.abhinav/resources/all_scripts/R/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.6
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  shh_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = shh_integrated_RNA_NN, dims = 30, res = res[i], 
                                                       saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA"), 
                                                       objname = "shh_sample_integrated_RNA_NN_cluster",
                                                       process = process, col_sel = c("orig.ident"))
}

saveRDS(shh_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/shh_integrated_RNA_NN_cluster_0.6.RDS",sep = ""))

