#### 2nd Time
### Some of the sample were not good quality
library(Seurat)
library(dplyr)
library(glmGamPoi)

obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.6.RDS")


DefaultAssay(obj) <- "RNA"
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
dir.create(savedir, showWarnings=FALSE)

cell_names <- rownames(obj@meta.data[grep("SF10961|SF12930|SF9232",obj@meta.data$orig.ident,invert=TRUE,value=FALSE),])
obj_subset <- subset(obj, cells = cell_names)

DefaultAssay(obj_subset) <- "RNA"
obj_subset <- NormalizeData(obj_subset, normalization.method = "LogNormalize", scale.factor = 10000)
obj_subset <- ScaleData(obj_subset, features = rownames(obj_subset))

s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
obj_subset <- CellCycleScoring(obj_subset, g2m.features = g2m_genes, s.features = s_genes)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
objname = "obj"
Assay = "RNA"
process = "sctransform"
CD4_RNA_sctransformed <- sctransform_V2_integration(obj = obj_subset, saveDir = savedir, ngenes = 4000,
                                                    regress = c("nCount_RNA"),
                                                    dims = 30,
                                                    Assay = Assay, process = process, objname = objname,
                                                    split_by = "orig.ident",
                                                    reference = c(2),
                                                    sample_tree = NULL)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_RNA_intergation.R")
objname = "obj"
process = "integration"
Assay = "RNA"
CD4_integrated_RNA <- RNA_integration(CD4_RNA_sctransformed, savedir, dims = 30, RNA_features = c("ATOH1","PTCH1"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

CD4_integrated_RNA_NN <- FindNeighbors(CD4_integrated_RNA, dims = 1:30)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/cluster_UMAP_QC.R")
res = c(1.2,1.4,1.6,1.8,2)
# res = 0.8
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD4_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = obj_impute, dims = 30, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA","S.Score","G2M.Score"), objname = "obj_subset",
                                                       process = process, col_sel = c("orig.ident","tumor_type","sex","age"))
}

obj_impute <- CD4_integrated_RNA_NN_cluster

#### Splitting the clusters
clusters <- levels(obj_impute@meta.data$seurat_clusters)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
  cellnames <- obj_impute@meta.data[grep(paste("^",clusters[i],"$",sep=""), obj_impute@meta.data$seurat_clusters),] %>% rownames()
  p <- DimPlot(obj_impute,
               cells.highlight = cellnames,
               reduction = "umap",
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(clusters[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/cluster_splitted_res_0.8.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","KIF21B","DISP3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")

nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP,GNearly, GNlate, Neuronal, Stem_cells, nonmalignant_cells)

DefaultAssay(CD4_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(CD4_integrated_RNA_NN_cluster, feature = genes, group.by = "seurat_clusters", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_marker_genes_reS_.pdf", width = 24, height =11)
p
dev.off()

dir.create(paste(savedir,"saveRDS_obj",sep=""),showWarning=FALSE)
saveRDS(CD4_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/obj_integrated_RNA_NN_cluster_0.6.RDS",sep = ""))

#### BO reference mapping
BO_reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/BO_human_cerebellar/human.dev_cerebellum.glutamatergic.lineage_cells_SCT.rds")
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "SCT"
DefaultAssay(BO_reference) <- "SCT"
anchors <- FindTransferAnchors(reference = BO_reference,
                                query = shh_integrated_RNA_NN_cluster,
                                normalization.method = "SCT",
                                dims = 1:50)
    
predictions.assay <- TransferData(anchorset = anchors,
                                  refdata = BO_reference$celltype,
                                  prediction.assay = TRUE,
                                  weight.reduction = shh_integrated_RNA_NN_cluster[["pca"]],
                                  dims = 1:30, k.weight = 50)
    
shh_integrated_RNA_NN_cluster[["BO_predictions"]] <- predictions.assay

celltypes <- as.vector(unique(BO_reference@meta.data$celltype))
    
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "BO_predictions"
pdf(paste(savedir,"/featureplot/featureplot_reference_BO_ref_mapped.pdf",sep = ""), width = 10, height = 9)
FeaturePlot(shh_integrated_RNA_NN_cluster, features = celltypes)
dev.off()
    
prediction <- TransferData(anchorset = anchors,
                           refdata = BO_reference$celltype,
                           prediction.assay = FALSE,
                           weight.reduction = shh_integrated_RNA_NN_cluster[["pca"]],
                           dims = 1:30)
                              
shh_integrated_RNA_NN_cluster <- AddMetaData(shh_integrated_RNA_NN_cluster, metadata = prediction)    
p <- DimPlot(shh_integrated_RNA_NN_cluster,group.by = "predicted.id", reduction = "umap")
    
dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_celltypes_umap_BO_ref.pdf",sep = ""))
print(p)
dev.off()

#### Human Cerebellar dev
reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/human_cerebellar_SCT.RDS")

DefaultAssay(shh_integrated_RNA_NN_cluster) <- "SCT"
DefaultAssay(reference) <- "SCT"
anchors <- FindTransferAnchors(reference = reference,
                                query = shh_integrated_RNA_NN_cluster,
                                normalization.method = "SCT",
                                dims = 1:50)
    
predictions.assay <- TransferData(anchorset = anchors,
                                  refdata = reference$fig_cell_type,
                                  prediction.assay = TRUE,
                                  weight.reduction = shh_integrated_RNA_NN_cluster[["pca"]],
                                  dims = 1:30, k.weight = 50)
    
shh_integrated_RNA_NN_cluster[["cerebellar_dev_predictions"]] <- predictions.assay

celltypes <- as.vector(unique(reference@meta.data$fig_cell_type))
    
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "cerebellar_dev_predictions"
pdf(paste(savedir,"/featureplot/featureplot_reference_human_cerebellum_ref_mapped.pdf",sep = ""), width = 10, height = 9)
FeaturePlot(shh_integrated_RNA_NN_cluster, features = celltypes)
dev.off()
    
prediction <- TransferData(anchorset = anchors,
                           refdata = reference$fig_cell_type,
                           prediction.assay = FALSE,
                           weight.reduction = shh_integrated_RNA_NN_cluster[["pca"]],
                           dims = 1:30)

human_cerebellum_pred <- prediction
                              
shh_integrated_RNA_NN_cluster <- AddMetaData(shh_integrated_RNA_NN_cluster, metadata = human_cerebellum_pred)    
p <- DimPlot(shh_integrated_RNA_NN_cluster,group.by = "predicted.id", reduction = "umap", 
             cols = c("H-Glia"='#1f77b4',"H-RL"='#ff7f0e',"H-GCP"='#279e68',"H-PC"='#d62728',"H-Meninges"='#aa40fc',
                          "H-eCN/UBC"='#8c564b',"H-BG"='#e377c2',"H-Ast"='#b5bd61',"H-Ast/Ependymal"='#b5bd61',"H-OPC"='#17becf',
                          "H-GN"='#aec7e8',"H-PIP"='#ffbb78',"H-Microglia"='cadetblue3',"H-iCN"='cornsilk4',"H-Choroid"='plum2',
                          "H-MLI"='yellow',"H-Endothelial"='hotpink',"H-Pericytes" = "black","H-Committed OPC" = "cyan"))
    
dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_celltypes_umap_cerebellum_reference_ref.pdf",sep = ""))
print(p)
dev.off()

DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"

genes <- c("EOMES","LMX1A","SMAD3","OTX2","UNC5D")
pdf(paste(savedir,"/featureplot/UBC_genes.pdf",sep = ""))
FeaturePlot(shh_integrated_RNA_NN_cluster, features = genes)
dev.off()

DefaultAssay(CD4_integrated_RNA_NN_cluster) <- "RNA"
markers <- FindAllMarkers(CD4_integrated_RNA_NN_cluster)

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","GRIN2C")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","KIF21B","DISP3","BRINP1")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]

genes <- c(GCP,GNearly, GNlate, Stem_cells, nonmalignant_cells)

DefaultAssay(obj) <- "RNA"
p <- DotPlot(obj, feature = genes) + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  + scale_size(range = c(3, 6)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/marker_genes.pdf", height = 20, width =9)
p
dev.off()


pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/cellcycling.pdf")
DimPlot(obj ,group.by = "Phase")
dev.off()

obj$CT_phase <- paste(obj$seurat_clusters, obj$Phase, sep = "_")
cycling_table <- table(obj$CT_phase) %>% as.data.frame()

cycling_table$CT <- paste("C",gsub("_.*.","",cycling_table$Var1),sep="")

data2 <- cycling_table[,-1]
# Calculate the percentages
data2 <- data2 %>%
  group_by(CT) %>%
  mutate(total = sum(Freq),
         percentage = (Freq / total) * 100) %>%
  select(-total)

data2$phase <- gsub(".*._","",cycling_table$Var1)
# Recast the table

write.table(recast_data, "/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/cycling_distribution.txt", sep = "\t", quote = F, col.names = T, row.names = T)

library(stringr)
sample_batch <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_batch.txt", header = TRUE)

GCP <- c(0,1,18,8)
GCP_cycling <- c(4,5,7,16)
GNearly <- c(6,14,2)
GNlate <- c(3,9,10,21)
Pericytes <- 12
Astrocytes <- 13
Immune <- 15
Oligo <- 17
OPC <- 20
Endo <- 19
Unknown <- 11

patterns <- c(GCP,GCP_cycling,GNearly,GNlate,Pericytes,Astrocytes,Immune,Oligo,OPC,Endo,Unknown)
replacements <- c(rep("GCP",length(GCP)),rep("GCP_cycling",length(GCP_cycling)),
  rep("GNearly",length(GNearly)),rep("GNlate",length(GNlate)), rep("Pericytes",length(Pericytes)),rep("Astrocytes",length(Astrocytes)),
  rep("Immune",length(Immune)),rep("Oligo",length(Oligo)),rep("OPC",length(OPC)),rep("Endo",length(Endo)),rep("Unknown",length(Unknown)))


obj@meta.data$celltypes <- obj@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  obj@meta.data$celltypes <- str_replace_all(obj@meta.data$celltypes, pattern, replacements[i])
}

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/celltype.pdf")
DimPlot(obj,group.by = "celltypes",
cols = c("GCP"='#1f77b4',"GCP_cycling"='#ff7f0e',"GNearly"='#279e68',"GNlate"='#d62728',"Pericytes"='#aa40fc',
                          "Astrocytes"='#8c564b',"Immune"='black',"Oligo"='#b5bd61',"OPC"='hotpink',"Endo"='#17becf',"Unknown"="gray33"))
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","GRIN2C","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]

genes <- c(GCP, Stem_cells, GNearly, GNlate, nonmalignant_cells)


DefaultAssay(obj) <- "RNA"
p <- DotPlot(obj, feature = genes, group.by = "celltypes") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  + scale_size(range = c(3, 6)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/marker_genes_celltype.pdf", height = 20, width =9)
p
dev.off()

#### Longitudinal analysis
celltypes_ind <- table(obj@meta.data$sample,obj@meta.data$celltypes)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_individual3.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_individual3.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
  df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample","celltype","percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c("7316_4529","7316_5881","7316_278","7316_333","7316_737","7316_2118",
  "7316_931","7316_3023","7316_2978","7316_311","7316_1666",
  "7316_1676","DOD4182","SF10961","SF12930","SF7994",
  "SF8368","SF8539","SF9232") -> sample_levels

df_melted$sample <- gsub("-","_",df_melted$sample)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#279e68','#d62728','#ff7f0e','#1f77b4','#aa40fc',
                               '#8c564b','#e377c2','#b5bd61','#17becf',
                               '#aec7e8','#ffbb78','darkseagreen3','cornsilk4','plum2',
                               'yellow', "black", "blanchedalmond","blue","white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype.pdf",width =10, height = 7)
p
dev.off()

p <- DimPlot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#279e68','#d62728','#ff7f0e','#1f77b4','#aa40fc',
                               '#8c564b','#e377c2','#b5bd61','#17becf',
                               '#aec7e8','#ffbb78','darkseagreen3','cornsilk4','plum2',
                               'yellow', "black", "blanchedalmond","blue","white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype.pdf",width =10, height = 7)
p
dev.off()

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",df_melted$sample),]
p <- ggplot(df_melted_paired, aes(fill=celltype, y=percentage, x=timepts)) + 
  geom_bar(position="fill", stat="identity") + 
  # geom_text(aes(label = paste(percentage,"%")), position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values = c('#279e68','#d62728','#ff7f0e','#1f77b4','#aa40fc',
                               '#8c564b','#e377c2','#b5bd61','#17becf',
                               '#aec7e8','#ffbb78','darkseagreen3','cornsilk4','plum2',
                               'yellow', "black", "blanchedalmond","blue","white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype_timepoint_paired.pdf",width =10, height = 7)
p
dev.off()


pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/celltype.pdf")
DimPlot(obj,group.by = "celltypes",
cols = c("GCP"='#ff7f0e',"GCP_cycling"='#1f77b4',"GNearly"='#aa40fc',"GNlate"='#8c564b',"Pericytes"='#aec7e8',
                          "Astrocytes"='#279e68',"Immune"='#e377c2',"Oligo"='#b5bd61',"OPC"='#17becf',"Endo"='#d62728',"Unknown"="#ffbb78"))
dev.off()

# p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
#   geom_bar(position="fill", stat="identity") + 
#   scale_fill_manual(values = c('#279e68','#d62728','#ff7f0e','#1f77b4','#aa40fc',
#                                '#8c564b','#e377c2','#b5bd61','#17becf',
#                                '#aec7e8','#ffbb78','darkseagreen3','cornsilk4','plum2',
#                                'yellow', "black", "blanchedalmond","blue","white")) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype_GCP.pdf",width =10, height = 7)
# p
# dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header= TRUE,sep="\t")

patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

### Patient ID
patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Patient.ID
# names(replacement) <- patterns
df_melted$Patient.ID <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$Patient.ID <- str_replace_all(df_melted$Patient.ID, pattern, replacements[i])
}

df_melted$celltype2 <- gsub("GCP_cycling","GCP",df_melted$celltype) %>% gsub("GNearly|GNlate","GN",.)

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",df_melted$sample),]
df_melted_paired$Patient_Timepts <- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts, sep = "_")
df_wide <- dcast(df_melted_paired, celltype ~ Patient_Timepts, value.var = "percentage", fun.aggregate=mean)

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test$p <- (stat.test$p)/2

bxp <- ggboxplot(
  df_melted_paired, x = "celltype", y = "percentage", 
  color = "timepts", palette = c("#00AFBB", "#E7B800")
) 

bxp2 <- bxp +   geom_dotplot(
  aes(fill = timepts, color = timepts), trim = FALSE,
  binaxis='y', stackdir='center', dotsize = 0.15,
  position = position_dodge(0.8)
)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +   theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white"))

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"Table/t_test_paired_primary_recurrent_samples_bxplot_two_sided_celltype2.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

#### Genotype #####
var_anno <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/total_variant_summary.txt", header=TRUE, sep = "\t")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Filter data to keep only Primary1 and Recurrent
filtered_data <- var_anno %>% filter(!is.na(Primary1) & !is.na(Recurrent))

# Perform paired Wilcox test and t-test for each var_type
test_results <- filtered_data %>%
  group_by(var_type) %>%
  summarize(
    wilcox_p = wilcox.test(Primary1, Recurrent, paired = TRUE)$p.value,
    ttest_p = t.test(Primary1, Recurrent, paired = TRUE, alternative = "less")$p.value
  )

# Reshape the data for plotting
plot_data <- filtered_data %>%
  gather(key = "Condition", value = "Value", -var_type, -Patients) %>%
  filter(Condition %in% c("Primary1", "Recurrent"))
  
p <- ggplot(plot_data, aes(x = Condition, y = Value, color = Condition)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Condition), position = position_dodge(width = 0.7), fill = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.25) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = Patients), size = 0.5, color = "black") +
  scale_color_manual(values = c("cornflowerblue", "lightcoral")) +
  scale_fill_manual(values = c("cornflowerblue", "lightcoral")) +
  facet_wrap(~var_type) +
  labs(title = "Mean Number of Variants by Type and Condition",
       x = "Variant Type",
       y = "Mean Number of Variants") +
  theme_minimal()  +  
  geom_text(data = test_results, aes(x = 1.5,label = sprintf("T-test p=%.3f", ttest_p), y = max(plot_data$Value) + 10),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 5, color = "black", inherit.aes = FALSE)

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/total_variant_information.pdf")
p
dev.off()


#### Since some of the variant can be different built i.e. hg19 and hg38
variant <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/checking_genome_version.txt", sep = '\t')
variant$gene <- gsub(":.*.","",variant$V3)

unmatch_variants <- variant[!(variant$V1 == variant$gene),]

write.table(unmatch_variants, "/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/unmatched_variants.txt", sep = "\t" ,quote=F, row.names = F, col.names=F)

#### Filtered variants
var_anno <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/Filtered_variant_summary.txt", header=TRUE, sep = "\t")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Filter data to keep only Primary1 and Recurrent
filtered_data <- var_anno %>% dplyr::filter(!is.na(Primary1) & !is.na(Recurrent))

# Perform paired Wilcox test and t-test for each var_type
test_results <- filtered_data %>%
  group_by(var_type) %>%
  summarize(
    wilcox_p = wilcox.test(Primary1, Recurrent, paired = TRUE)$p.value,
    ttest_p = t.test(Primary1, Recurrent, paired = TRUE, alternative = "less")$p.value,
    ttest_u = t.test(Primary1, Recurrent, paired = FALSE, alternative = "less")$p.value
  )

# Reshape the data for plotting
plot_data <- filtered_data %>%
  gather(key = "Condition", value = "Value", -var_type, -Patients) %>%
  filter(Condition %in% c("Primary1", "Recurrent"))
  
p <- ggplot(plot_data, aes(x = Condition, y = Value, color = Condition)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Condition), position = position_dodge(width = 0.7), fill = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.25) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = Patients), size = 0.5, color = "black") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  facet_wrap(~var_type) +
  labs(title = "Mean Number of Variants by Type and Condition",
       x = "Variant Type",
       y = "Mean Number of Variants") +
  theme_minimal()  +  
  geom_text(data = test_results, aes(x = 1.5,label = sprintf("T-test p=%.3f", ttest_p), y = max(plot_data$Value) + 2.5),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 5, color = "black", inherit.aes = FALSE)

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_information.pdf")
p
dev.off()

gene_table <- list.files("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/", pattern = "_genes_2.txt", recursive = TRUE, full.names = TRUE)
gene_table <- grep("PT_H3WWDMW9",gene_table,invert=TRUE,value=TRUE)
samplename <- gsub("_genes_2.txt","",basename(gene_table))


for(i in 1:length(gene_table)){
  try({
  sample_table <- read.table(gene_table[i])
  sample_table$samplename <- samplename[i]
  colnames(sample_table) <- c("gene_number","gene","samplename")
  assign(paste(samplename[i],"_gene",sep=""), sample_table)
  })
}

sample_tablename <- ls(pattern = "_gene")

combined_table <- rbind(get(sample_tablename[1]), get(sample_tablename[2]), get(sample_tablename[3]), get(sample_tablename[4]), get(sample_tablename[5]),
      get(sample_tablename[6]),get(sample_tablename[7]), get(sample_tablename[8]), get(sample_tablename[9]), get(sample_tablename[10]))

library(reshape2)
result <- combined_table %>%
  pivot_wider(names_from = gene, values_from = gene_number, values_fill = list(gene_number = 0))
  
write.table(result, "/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/final_result.txt", quote=F, row.names = F, col.names = T, sep = "\t")

heatmap_matrix <- as.data.frame(result[,-1])
row.names(heatmap_matrix) <- result$samplename

# Create the heatmap
library(ComplexHeatmap)
heatmap_matrix2 <- heatmap_matrix[c(1,2,5,6,8,7,10,9,3,4),]
p <- Heatmap(as.matrix(heatmap_matrix2),  
        name = "Gene Number", 
        col = colorRampPalette(c("white", "red"))(50),
        cluster_rows = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal"))

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/variant_heatmap.pdf", width =11, height =6)
draw(p, heatmap_legend_side = "bottom")
dev.off()

### RNA GATK variant
#### Filtered variants
var_anno <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_summary.txt", header=TRUE, sep = "\t")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Filter data to keep only Primary1 and Recurrent
filtered_data <- var_anno %>% dplyr::filter(!is.na(Primary1) & !is.na(Recurrent))

# Perform paired Wilcox test and t-test for each var_type
test_results <- filtered_data %>%
  group_by(var_type) %>%
  summarize(
    wilcox_p = wilcox.test(Primary1, Recurrent, paired = TRUE)$p.value,
    ttest_p = t.test(Primary1, Recurrent, paired = TRUE, alternative = "less")$p.value,
    ttest_u = t.test(Primary1, Recurrent, paired = FALSE, alternative = "less")$p.value
  )

# Reshape the data for plotting
plot_data <- filtered_data %>%
  gather(key = "Condition", value = "Value", -var_type, -Samples) %>%
  filter(Condition %in% c("Primary1", "Recurrent"))
  
p <- ggplot(plot_data, aes(x = Condition, y = Value, color = Condition)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Condition), position = position_dodge(width = 0.7), fill = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.25) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = Samples), size = 0.5, color = "black") +
  scale_color_manual(values = c("cornflowerblue", "lightcoral")) +
  scale_fill_manual(values = c("cornflowerblue", "lightcoral")) +
  facet_wrap(~var_type) +
  labs(title = "Mean Number of Variants by Type and Condition",
       x = "Variant Type",
       y = "Mean Number of Variants") +
  theme_minimal()  +  
  geom_text(data = test_results, aes(x = 1.5,label = sprintf("T-test p=%.3f", ttest_p), y = max(plot_data$Value) + 2.5),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 5, color = "black", inherit.aes = FALSE)

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_information.pdf", width = 4)
p
dev.off()

### All the samples
var_anno <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/nonsynon_synon_ratio.txt", header=TRUE, sep = "\t")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Filter data to keep only Primary1 and Recurrent
# filtered_data <- var_anno %>% dplyr::filter(!is.na(Primary1) & !is.na(Recurrent))
filtered_data <- var_anno 

# Perform paired Wilcox test and t-test for each var_type
test_results <- filtered_data %>%
  group_by(var_type) %>%
  summarize(
    wilcox_p = wilcox.test(Primary1, Recurrent, paired = FALSE)$p.value,
    ttest_p = t.test(Primary1, Recurrent, paired = FALSE, alternative = "two.sided")$p.value,
    ttest_u = t.test(Primary1, Recurrent, paired = FALSE, alternative = "two.sided")$p.value
  )

# Reshape the data for plotting
plot_data <- filtered_data %>%
  gather(key = "Condition", value = "Value", -var_type, -Samples) %>%
  filter(Condition %in% c("Primary1", "Recurrent"))
  
p <- ggplot(plot_data, aes(x = Condition, y = Value, color = Condition)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Condition), position = position_dodge(width = 0.7), fill = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.25) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = Samples), size = 0.5, color = "black") +
  scale_color_manual(values = c("cornflowerblue", "lightcoral")) +
  scale_fill_manual(values = c("cornflowerblue", "lightcoral")) +
  facet_wrap(~var_type) +
  labs(title = "Synonymous to Non Synonymous Variants Ratio",
       x = "Tumor Type",
       y = "Non Synonymous to Synonymous Ratio") +
  theme_minimal()  +  
  geom_text(data = test_results, aes(x = 1.5,label = sprintf("T-test p=%.3f", ttest_p), y = max(plot_data$Value) + 0.5),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 5, color = "black", inherit.aes = FALSE)

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/nonsynonymous_to_synonymous_ratio_paired.pdf", width = 4)
p
dev.off()

### Synonymous and Nonsynonymous variant
### All the samples
var_anno <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_summary_2.txt", header=TRUE, sep = "\t")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Filter data to keep only Primary1 and Recurrent
filtered_data <- var_anno 

# Perform paired Wilcox test and t-test for each var_type
test_results <- filtered_data %>%
  group_by(var_type) %>%
  summarize(
    wilcox_p = wilcox.test(Primary, Recurrent, paired = TRUE)$p.value,
    ttest_p = t.test(Primary, Recurrent, paired = FALSE, alternative = "less")$p.value,
    ttest_u = t.test(Primary, Recurrent, paired = FALSE, alternative = "less")$p.value
  )

# Reshape the data for plotting
plot_data <- filtered_data %>%
  gather(key = "Condition", value = "Value", -var_type, -Samples) %>%
  filter(Condition %in% c("Primary", "Recurrent"))
  
p <- ggplot(plot_data, aes(x = Condition, y = Value, color = Condition)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Condition), position = position_dodge(width = 0.7), fill = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.25) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  # geom_line(aes(group = Samples), size = 0.5, color = "black") +
  scale_color_manual(values = c("cornflowerblue", "lightcoral")) +
  scale_fill_manual(values = c("cornflowerblue", "lightcoral")) +
  facet_wrap(~var_type) +
  labs(title = "Mean Number of Variants by Type and Condition",
       x = "Variant Type",
       y = "Mean Number of Variants") +
  theme_minimal()  +  
  geom_text(data = test_results, aes(x = 1.5,label = sprintf("T-test p=%.3f", ttest_p), y = max(plot_data$Value) + 2.5),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 5, color = "black", inherit.aes = FALSE)

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_information_all_samples.pdf", width = 4)
p
dev.off()


### No Frameshift
var_anno <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_summary_noframeshift.txt", header=TRUE, sep = "\t")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Filter data to keep only Primary1 and Recurrent
filtered_data <- var_anno %>% dplyr::filter(!is.na(Primary1) & !is.na(Recurrent))

# Perform paired Wilcox test and t-test for each var_type
test_results <- filtered_data %>%
  group_by(var_type) %>%
  summarize(
    wilcox_p = wilcox.test(Primary1, Recurrent, paired = TRUE)$p.value,
    ttest_p = t.test(Primary1, Recurrent, paired = TRUE, alternative = "less")$p.value,
    ttest_u = t.test(Primary1, Recurrent, paired = FALSE, alternative = "less")$p.value
  )

# Reshape the data for plotting
plot_data <- filtered_data %>%
  gather(key = "Condition", value = "Value", -var_type, -Samples) %>%
  filter(Condition %in% c("Primary1", "Recurrent"))
  
p <- ggplot(plot_data, aes(x = Condition, y = Value, color = Condition)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Condition), position = position_dodge(width = 0.7), fill = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.25) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = Samples), size = 0.5, color = "black") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  facet_wrap(~var_type) +
  labs(title = "Mean Number of Variants by Type and Condition",
       x = "Variant Type",
       y = "Mean Number of Variants") +
  theme_minimal()  +  
  geom_text(data = test_results, aes(x = 1.5,label = sprintf("T-test p=%.3f", ttest_p), y = max(plot_data$Value) + 2.5),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 5, color = "black", inherit.aes = FALSE)

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/Filtered_variant_information_noframeshift.pdf")
p
dev.off()

gene_table <- list.files("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/", pattern = "_gene_only.txt", recursive = TRUE, full.names = TRUE)
#gene_table <- grep("PT_H3WWDMW9",gene_table,invert=TRUE,value=TRUE)
samplename <- gsub("_annovar_result.hg38_multianno_selected_col_combined_selected_gene_only.txt","",basename(gene_table))


for(i in 1:length(gene_table)){
  try({
  sample_table <- read.table(gene_table[i])
  sample_table$samplename <- samplename[i]
  colnames(sample_table) <- c("gene_number","gene","samplename")
  assign(paste(samplename[i],"_gene",sep=""), sample_table)
  })
}

sample_tablename <- ls(pattern = "_gene")

combined_table <- rbind(get(sample_tablename[1]), get(sample_tablename[2]), get(sample_tablename[3]), get(sample_tablename[4]), get(sample_tablename[5]),
      get(sample_tablename[6]),get(sample_tablename[7]), get(sample_tablename[8]), get(sample_tablename[9]), get(sample_tablename[10]),get(sample_tablename[11]),
      get(sample_tablename[12]),get(sample_tablename[13]), get(sample_tablename[14]), get(sample_tablename[15]), get(sample_tablename[16]),get(sample_tablename[17]))

library(reshape2)
library(tidyr)
library(dplyr)
result <- combined_table %>%
  pivot_wider(names_from = samplename, values_from = gene_number, values_fill = list(gene_number = 0))

write.table(result, "/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/final_result.txt", quote=F, row.names = F, col.names = T, sep = "\t")

# heatmap_matrix <- as.data.frame(result[,-1])
# row.names(heatmap_matrix) <- result$gene

heatmap_matrix <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/final_result.txt", header = TRUE, row.names = 1)
# Create the heatmap
library(ComplexHeatmap)
grep("SF12930",colnames(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[,-14] ## Please check if this SF12930
samplename <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/samplename.txt")
samplename$V1 <- gsub("7316","X7316",samplename$V1)
all(colnames(heatmap_matrix) == samplename$V1)

colnames(heatmap_matrix) <- samplename$V2

genes <- c("ATOH1","PCNA","JAG1","NES","SOX2","SMO","SUFU","PTCH1","PTCH2","GLI1","GLI2","MYCN","RELN","TERT","ELP1","DDX3X","EYA1",
           "KMT2D","BCOR","GSE1","PTEN","CREBBP","KMT2C","SHH","PRKACB","MYC","GNAS","PPM1D","CCDN1","CCND2","RBFOX3",
           "NEUROD1","ZIC1","ZIC2","TUBB3","MEG3","GRIN2B","GRIN2C","CHD7","MGP","TCERG1L","TP53","MKI67","TOP2A")

prim_index <- grep("Primary",colnames(heatmap_matrix))
recur_index <- grep("Recurrent",colnames(heatmap_matrix))
heatmap_matrix2 <- heatmap_matrix[genes,c(prim_index,recur_index)]
heatmap_matrix3 <- na.omit(heatmap_matrix2)
# heatmap_matrix4 <- heatmap_matrix3[,c(1,2,3,4,8,7,9,10,5,6,11,12,13,14)]

p <- Heatmap(as.matrix(heatmap_matrix3), 
        name = "Gene Number", 
        col = colorRampPalette(c("white", "red"))(50),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal"))

pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/variant_heatmap_unclus.pdf", width =6, height =8)
p
dev.off()

#### RAS RAF pathway
c("MAP2K1","MAP2K2","MAPK3","MAPK1","KRAS","NRAS","HRAS",)

### Preparing an Oncoprint
# grep "[0-9]" */*_multianno_selected_col_combined_selected_gene.txt | sed 's/:/\t/1' | cut -f33-37 | sort | uniq > all_variants
# grep -wFf all_variants */*.vcf  | sed 's/:/\t/1' | cut -f2-6 | sort | uniq -c | awk '{if($1>1) print $0}' > variant_sample_repeating.txt
# grep -wFf remaining_variants */*.vcf | awk '{print $NF}' | sed 's/:/\t/g' | cut -f2 | sed 's/,/\t/g' | awk '{$3=$1+$2; print $3}' | awk '{print NR" "$0}' | awk '{if($2 > 19) print $0}' | awk '{print $1"p;"}' | tr -s "\n" "\t" | sed 's/\t//g'
# sed -n '2p;7p;19p;22p;23p;28p;29p;34p;35p;36p;37p;41p;42p;44p;45p;46p;47p;50p;51p;53p;54p;55p;57p;58p;62p;64p;65p;66p;70p;71p;72p;74p;75p;76p;77p;78p;79p;84p;87p;88p;89p;93p;94p;95p;96p;97p;98p;99p;103p;104p;105p;107p;109p;110p;111p;112p;113p;117p;121p;122p;126p;130p;131p;134p;136p;137p;139p;141p;143p;144p;145p;147p;148p;149p;152p;153p;155p;156p;157p;159p;160p;161p;167p;168p;169p;170p;172p;173p;176p;177p;179p;182p;183p;186p;188p;189p;190p;194p;195p;199p;200p;204p;205p;206p;207p;215p;216p;217p;218p;219p;220p;223p;224p;226p;228p;229p;230p;231p;232p;233p;234p;235p;241p;249p;250p;251p;255p;260p;269p;270p;277p;278p;282p;283p;289p;290p;291p;294p;297p;298p;299p;302p;303p;304p;305p;306p;307p;308p;310p;311p;312p;313p;314p' remaining_variants.vcf
# Generating the dataframe for the genes and variants
gene_variant_file <- list.files("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/", 
           recursive = TRUE, 
           pattern = "gt_eq_20_final_genelist_gene_type.txt", full.names = TRUE)

samplename <- gsub("_annovar_result.hg38_multianno_selected_col_combined_allele_count_gt_eq_20_final_genelist_gene_type.txt","",basename(gene_variant_file))

for(i in 1:length(samplename)){
  try({
  variant_file <- read.table(gene_variant_file[i], header=FALSE, sep ="\t")
  colnames(variant_file) <- c("gene","variant")
  variant_file$samplename <- samplename[i]
  assign(paste("variant_file_",samplename[i],sep=""), variant_file)
  })
}

variant_files <- ls(pattern="variant_file_")

variants <- rbind(get(variant_files[1]),get(variant_files[2]),get(variant_files[3]),get(variant_files[4]),get(variant_files[5]),get(variant_files[6]),
      get(variant_files[7]),get(variant_files[8]),get(variant_files[9]),get(variant_files[10]),get(variant_files[11]),get(variant_files[12]),
      get(variant_files[13]))

variants$variant <- gsub(" ","_",variants$variant)

df <- variants 
df <- df[grep("^synonymous_SNV",df$variant,invert=TRUE),]
df <- df[grep("^NF1",df$gene,invert=TRUE),]

library(tidyr)
library(dplyr)
df_grouped <- df %>%
  group_by(gene, samplename) %>%
  summarise(variant = paste(variant, collapse = "; ")) %>%
  ungroup()

df_wide <- df_grouped %>%
  pivot_wider(names_from = samplename, values_from = variant)

df_wide <- as.data.frame(df_wide)
rownames(df_wide) <- (df_wide$gene)
df_wide_2 <- df_wide[,-1]
df_wide_2[is.na(df_wide_2)] <- ""

colnames(df_wide_2)[grep("DOD4182",colnames(df_wide_2))] <- "7316_4162"

write.table(df_wide_2, "/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/variant_calling.txt", quote =F, col.names= T, row.names = T, sep = "\t")

df_wide_2 <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/variant_calling_2.txt", header= TRUE, row.names = 1, sep = "\t")
# Load necessary libraries
library(ComplexHeatmap)

# Data preparation
# metadata <- read.csv("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/oncoprint.csv", stringsAsFactors = FALSE)
metadata <- read.table("/diazlab/data3/.abhinav/projects/SHH/metadata/Project_overview_oncoprint.txt", sep = "\t", stringsAsFactors = FALSE, header= TRUE)

metadata$Sample.ID<- gsub("-","_",metadata$Sample.ID)
metadata$Sample.ID <- gsub("SF7994R","SF7994",metadata$Sample.ID)
rem_sample <- grep(paste(colnames(df_wide_2),sep="",collapse="|"), metadata$Sample.ID,value=TRUE,invert=TRUE)
empty_df <- data.frame(matrix(ncol = length(rem_sample), nrow = nrow(df_wide_2)))
colnames(empty_df) <- rem_sample
rownames(empty_df) <- rownames(df_wide_2)

merge_df <- cbind(empty_df, df_wide_2)
merge_df[is.na(merge_df)] <- ""

metadata$Sample.ID <- gsub("7316","X7316",metadata$Sample.ID)
metadata1 <- metadata[match(colnames(df_wide_2),metadata$Sample.ID,nomatch=0),]

# Load necessary libraries
library(ComplexHeatmap)

# Ensure relevant columns are character vectors
metadata1$snRNA.seq <- as.character(metadata1$snRNA.seq)
metadata1$snATAC.seq <- as.character(metadata1$snATAC.seq)
metadata1$Visium <- as.character(metadata1$Visium)

# metadata1$snRNA.seq <- gsub("TRUE","Available",metadata1$snRNA.seq)
# metadata1$snATAC.seq <- gsub("TRUE","Available",metadata1$snATAC.seq)
# metadata1$Visium <- gsub("TRUE","Available",metadata1$Visium)

# metadata1[is.na(metadata1)] <- ""

stopifnot(all(metadata1$Sample.ID == colnames(df_wide_2)))
metadata1$Tumor.Site <- gsub("Cerebellum/Posterior Fossa","Cerebellum",metadata1$Tumor.Site)
column_anno <- HeatmapAnnotation(Type = metadata1$Tumor.type,
                                 Tumor.Site = metadata1$Tumor.Site,
                                 Status = metadata1$Status,
                                 snRNA = metadata1$snRNA.seq,
                                 snATAC = metadata1$snATAC.seq,
                                 Visium = metadata1$Visium,
                                 Sex = metadata1$Sex,
                                 Age = anno_points(metadata1$Age, ylim = c(0, max(metadata1$Age, na.rm = TRUE)), axis = TRUE),
                                 col = list(Type = c("Recurrent" = "#E7B800","Primary"="#00AFBB"," N.A."="grey"),
                                            Tumor.Site = c("Cerebellum" = "indianred1", "Ventricles"="lightskyblue3"," N.A."="grey"),
                                            Status = c("Deceased" = "lightpink1", "Alive" = "lightgreen"," N.A."="grey"),
                                            snRNA = c("Available" = "darkgreen", " N.A."="grey"),
                                            snATAC = c("Available" = "darkgreen", " N.A."="grey"),
                                            Visium = c("Available" = "darkgreen"," N.A."="grey"),
                                            Sex = c("M" = "Blue","F"="Pink")),
                                annotation_height = unit(c(3,3,3,3,3,3,3,8), "mm")  
)

col = c(frameshift_deletion = "mediumorchid2", nonsynonymous_SNV = "steelblue1", frameshift_insertion = "lightpink2", "stopgain"="red4")

ht1 <- oncoPrint(df_wide_2,
    alter_fun = list(
      background = function(x, y, w, h) grid.rect(x, y, w, h,
            gp = gpar(fill = "white")),
      frameshift_deletion = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.5, 
            gp = gpar(fill = col["frameshift_deletion"], col = NA), just = c("right", "top")),
        frameshift_insertion = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.5, 
            gp = gpar(fill = col["frameshift_insertion"], col = NA),just = c("left", "top") ),
        nonsynonymous_SNV = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.5, 
            gp = gpar(fill = col["nonsynonymous_SNV"], col = NA) , just = c("right", "bottom")),
        stopgain = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.5, 
            gp = gpar(fill = col["stopgain"], col = NA), just = c("left", "bottom"))), 
    col = col,
    show_column_names = TRUE,
    bottom_annotation = column_anno)

# pdf("/diazlab/data3/.abhinav/projects/SHH/Genotype/RNA_variants/SHH/GATK_RNAseq_Variant_Calling/oncoprint_snRNA_RNA_varaints.pdf", width = 8, height =9)
# ht1
# dev.off()

# pdf("/diazlab/data3/.abhinav/projects/SHH/metadata/oncoprint_snRNA_variants.pdf")
# ht1
# dev.off()

tiff("/diazlab/data3/.abhinav/projects/SHH/metadata/oncoprint_snRNA_variants.tiff", units="in", width=10, height=9, res=600)
ht1
dev.off()


tiff("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/sub_celltype_0.8_res_finalcelltypes.tiff", units="in", width=5.8, height=5, res=600)
DimPlot(RNA_integrated,group.by = "final_celltypes", label = TRUE,
cols = c("GCP_cyc"='darkslategray2',"GCP_Rib"='darkturquoise',"GN_Rib"='#436957',"GN_premig"='#5fe3a7',"GCP_HSP" = "deepskyblue4",
         "GN_mig"='#30b377',"GN_postmig" = "#107044","GCP" = "#6fb9ed","Pericyte"='#aa40fc',
         "Astro"='#8c564b',"Immune"='black',"Oligo"='#b5bd61',"OPC"='hotpink',
         "Endothelial"='#17becf',"GCP-GN_int" = "cyan4","Other" = "brown3"))
dev.off()


#### CNVs
### Generating script for each sample
# List of sample names
dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/")
setwd("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/")
samples <- c("obj_7316_1666", "obj_7316_1676", "obj_7316_2118", "obj_7316_278",
             "obj_7316_2978", "obj_7316_3023", "obj_7316_311", "obj_7316_333",
             "obj_7316_4529", "obj_7316_5881", "obj_7316_737", "obj_7316_931")

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sgene_position_add_cytoband.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation_marker_2.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- paste("X",gsub("-", ".", row.names(annot_file)),sep="")
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/markerannot_P_Q_arm",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/markerannot_P_Q_arm_2",  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             analysis_mode = "subclusters",
                             denoise=T,
                             tumor_subcluster_partition_method = "leiden",
                             HMM=T,
                             HMM_type = "i6",
                             plot_steps = T,
                             write_expr_matrix = T,
                             write_phylo = T,
                             no_prelim_plot=TRUE,
                             output_format = "pdf",
                             num_threads = 32)

saveRDS(infercnv_obj, "%s%s/CNV/markerannot_P_Q_arm/%s_infercnv.RDS")

', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)
  
  # Write the script content to a file
  script_file <- paste0("marker_script_", sample, ".R")
  writeLines(script_content, script_file)
}

### for other samples
samples <- c("obj_DOD4182","obj_SF7994", "obj_SF8368", "obj_SF8539")

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
num_cores <- detectCores()
cat("Number of cores detected:", num_cores, "\n")

# Register the parallel backend
# cl <- makeCluster(10)
# registerDoParallel(cl)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sgene_position_add_cytoband.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation_marker_2.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", row.names(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/markerannot_P_Q_arm",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/markerannot_P_Q_arm_2",  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             analysis_mode = "subclusters",
                             denoise=T,
                             tumor_subcluster_partition_method = "leiden",
                             HMM=T,
                             HMM_type = "i6",
                             plot_steps = T,
                             write_expr_matrix = T,
                             write_phylo = T,
                             no_prelim_plot=TRUE,
                             output_format = "pdf",
                             num_threads = 32)

saveRDS(infercnv_obj, "%s%s/CNV/markerannot_P_Q_arm/%s_infercnv.RDS")

# Stop the cluster
# stopCluster(cl)
', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)
  
  # Write the script content to a file
  script_file <- paste0("marker_script_", sample, ".R")
  writeLines(script_content, script_file)
}

### Running it on Slurm
# at this location /diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV
# marker_submit_template.sh
# #!/bin/bash
# #SBATCH --job-name=%s
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=50G
# #SBATCH --time=1-00:00:00
# #SBATCH --output=%s_%j.out
# #SBATCH --error=%s_%j.err
# #SBATCH --mail-user=abhinav.jain@ucsf.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50

# source ~/.bashrc
# conda activate
# R CMD BATCH marker_script_%s.R

#### cat submit_job.sh
#!/bin/bash

# List of sample names
samples=("obj_7316_2978" "obj_7316_3023" "obj_7316_5881" "obj_DOD4182" "obj_SF7994" "obj_SF8368" "obj_SF8539" "obj_7316_1666" "obj_7316_1676" "obj_7316_2118" "obj_7316_278" "obj_7316_2978" "obj_7316_3023" "obj_7316_311" "obj_7316_333" "obj_7316_4529" "obj_7316_5881" "obj_7316_737" "obj_7316_931")

# Base directory path for R scripts
r_script_dir="."

# Loop through each sample and create and submit the SLURM job script
for sample in "${samples[@]}"; do
  r_script="${r_script_dir}/script_${sample}.R"

  # Create the SLURM job script from the template
  slurm_script="submit_${sample}.sh"
  cp marker_submit_template.sh $slurm_script
  sed -i "s/%s/${sample}/g" $slurm_script
  sed -i "s|%s|${r_script}|g" $slurm_script

  # Submit the SLURM job
  sbatch $slurm_script
done

#### Calculating the score
library(dplyr)
source("/diazlab/data3/.abhinav/resources/functions/R/score_calc.R")
full_path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", pattern = ".hmm_mode-subclusters.observations.txt", recursive = TRUE, full.names = TRUE)
infercnv_dir <- gsub("infercnv.17_HMM_predHMMi6.*.","",full_path)
samplename <- gsub(".*.subset//","",full_path) %>% gsub("/CNV/.*.","",.)
xlim <- c(2000, 2000, 2500, 2000, 
            2000, 2000, 2500, 1000, 
            1000, 2000, 2200, 3000,
            5000, 2000, 2000, 1500)

for(i in 1:length(samplename)){
  final_score <- calc_scores(samplename[i],infercnv_dir[i], xlim[i])
  assign(paste("cnvscore",samplename[i],sep=""),final_score)
}

cnv_score <- ls(pattern = "cnvscore")

rbind(get(cnv_score[1]),get(cnv_score[2]), get(cnv_score[3]), get(cnv_score[4]), get(cnv_score[5]), get(cnv_score[6]),
      get(cnv_score[7]), get(cnv_score[8]),get(cnv_score[9]), get(cnv_score[10]), get(cnv_score[11]), get(cnv_score[12]), get(cnv_score[13]),
      get(cnv_score[14]), get(cnv_score[15])) -> combined_cnvscore

write.table(combined_cnvscore, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/combined_cnv_score_2.txt", quote=F, row.names=T, col.names=T, sep = "\t")

# Guassian distribution for K means clustering to separate out the tumnor and normal
score  <- get(cnv_score[i])

# set.seed(42) 
# kmeans_result <- kmeans(score$cnv_score, centers = 2)

# score$cluster <- kmeans_result$cluster

# normal_cluster <- which.min(kmeans_result$centers)
# tumor_cluster <- which.max(kmeans_result$centers)

# # Get the cutoff as the maximum value in the normal cluster
# cutoff <- max(score$cnv_score[score$cluster == normal_cluster])
# cat("Cutoff value for separating normal and tumor cells:", cutoff, "\n")

### K means giving very far fetched using the logistic regression

library(glmnet)
library(pROC)
library(ggplot2)
library(ggpubr)
library(ggrepel)

objname <- gsub("cnvscore","",cnv_score)
for(i in 1:length(cnv_score)){
  score  <- get(cnv_score[i])

  score$type_binary <- ifelse(score$type == "Tumor", 1, 0)
  score$type_binary <- as.factor(score$type_binary)

  logistic_model <- glm(type_binary ~ cnv_score, data = score, family = "binomial")
  predicted_probs <- predict(logistic_model, type = "response")
  
  roc_obj <- roc(score$type, predicted_probs)
  optimal_cutoff <- coords(roc_obj, "best", ret = "threshold",transpose = TRUE)
  
  pdf(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"ROC_curve.pdf",sep =""))
  plot.roc(roc_obj, main = "ROC Curve for CNV score")
  dev.off()
  
  score$predicted_prob <- predicted_probs
  pdf(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"ROC_line.pdf",sep =""))
  ggplot(score, aes(x=predicted_prob, y=cnv_score)) + geom_line()
  dev.off()
  
  # Perform interpolation
  y_value <- approx(score$predicted_prob, score$cnv_score, xout = optimal_cutoff[["threshold"]])$y  
  # fit_values_cnvscore <- data.frame(logistic_model$data$cnv_score,logistic_model$fitted.values)
  # cutoff_values <- fit_values_cnvscore[between(fit_values_cnvscore$logistic_model.fitted.values, round(optimal_cutoff[["threshold"]],2)-0.1, round(optimal_cutoff[["threshold"]],2)+0.1),"logistic_model.data.cnv_score"] %>% unique()
  
  # y_positions <- rep(0, length(cutoff_values))  # Adjust the y positions if needed
  
  # labels_df <- data.frame(cutoff_values = cutoff_values, y_positions = y_positions)
  
  p <- gghistogram(score,
        x = "cnv_score", y = "..density..",
        add = "mean", rug = FALSE,
        fill = "type", palette = c("#848484", "#F9BF1A"),
        add_density = TRUE, bins = 150, xlim = c(0, 2000)
    ) +
        geom_vline(xintercept = y_value, linetype = "dashed", color = "#62A7D5") +
        geom_text(aes(x = y_value, y = 0.002, label = y_value)) +
        ggtitle(objname[i])
        
  pdf(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"hist_with_cutoff_2.pdf",sep=""))
  print(p)
  dev.off()

  write.table(y_value, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"_cutoff_values.txt",sep = ""))

  ## and then extract just the information that we want from that variable.
  roc.df <- data.frame(
  tpp=roc_obj$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc_obj$specificities)*100, ## fpp = false positive precentage
  thresholds=roc_obj$thresholds)
    
  ## Calculate the area under the curve...
  pdf(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"ROC_curve_wid_AUC.pdf",sep =""))
  roc(score$type_binary, logistic_model$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
  dev.off()
}

cutoff_val <- as.numeric(c("396","932","501","270","469.5","442.5","996","23.5","240","487","569","1700","3236.5","399","726.5","532"))

for(i in 1:length(cnv_score)){
  scores <- get(cnv_score[i])
  scores$predicted_CNV <- "normal"
  scores[scores$cnv_score > cutoff_val[i],"predicted_CNV"] <- "Tumor"
  conf_mat <- table(scores$predicted_CNV, scores$type)
  assign(paste(objname[i],"_CNV",sep=""),scores)
  write.table(conf_mat, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"_confusion_mat_2.txt",sep =""))
  write.table(scores, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",objname[i],"/CNV/markerannot_P_Q_arm_2/",objname[i],"_cnv_score_df.txt",sep =""))
}

combined_CNVscore <- rbind(get(CNV_score[1]),get(CNV_score[2]),get(CNV_score[3]),get(CNV_score[4]),get(CNV_score[5]),get(CNV_score[6]),
      get(CNV_score[7]),get(CNV_score[8]),get(CNV_score[9]),get(CNV_score[10]),get(CNV_score[11]),get(CNV_score[12]),
      get(CNV_score[13]),get(CNV_score[14]),get(CNV_score[15]),get(CNV_score[16]))

write.table(combined_CNVscore, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/combined_CNVscore_931_cutoff.txt",quote=T, row.names = T, col.names = T, sep ="\t")

combined_CNVscore<- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/combined_CNVscore.txt", header=TRUE, row.names = 1)
combined_CNVscore_ordered <- combined_CNVscore[match( rownames(obj_impute@meta.data), rownames(combined_CNVscore),nomatch=0),]
stopifnot(rownames(combined_CNVscore_ordered) == rownames(obj_impute@meta.data)) ## Sanity Check
obj_impute@meta.data$cnv_score <- combined_CNVscore_ordered$cnv_score
obj_impute@meta.data$predicted_CNVtype <- combined_CNVscore_ordered$predicted_CNV
obj_impute@meta.data$cnv_type <- combined_CNVscore_ordered$type

p <- DimPlot(obj_impute,group.by = "predicted_CNVtype", reduction = "umap", split.by = "predicted_CNVtype",
             cols = c("normal"='gray60',"Tumor"='Red'))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/"    
dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_predicted_CNV_2.pdf",sep = ""))
print(p)
dev.off()

p <- DimPlot(obj_impute,group.by = "predicted_CNVtype", reduction = "umap",
             cols = c("normal"='gray60',"Tumor"='Red'))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/"    
dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_predicted_CNV_combined.pdf",sep = ""))
print(p)
dev.off()

obj_impute@meta.data$predicted_CNVtype2 <- obj_impute@meta.data$predicted_CNVtype
obj_impute@meta.data[grep("7316-1666|7316-1676|7316-2118|7316-333",rownames(obj_impute@meta.data)),"predicted_CNVtype2"] <- NA

p <- DimPlot(obj_impute,group.by = "predicted_CNVtype2", reduction = "umap", split.by = "predicted_CNVtype2",
             cols = c("normal"='gray60',"Tumor"='Red', "Confused" = "darkgoldenrod4"))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/"    
dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_predicted_CNV_confused.pdf",sep = ""))
print(p)
dev.off()

obj_impute@meta.data$predicted_CNVtype2 <- obj_impute@meta.data$predicted_CNVtype
obj_impute@meta.data[grep("7316-1666|7316-1676|7316-2118|7316-333",rownames(obj_impute@meta.data)),"predicted_CNVtype2"] <- NA

p <- DimPlot(obj_impute,group.by = "predicted_CNVtype2", reduction = "umap",
             cols = c("normal"='gray60',"Tumor"='Red'))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/"    
dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_predicted_CNV_confused_combined.pdf",sep = ""))
print(p)
dev.off()

### Identifying the purity of cancer
cm_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", pattern = "_confusion_mat_2.txt", recursive = TRUE, full.names = TRUE)
samplename <- gsub("_confusion_mat.txt","",basename(cm_files)) %>% gsub("obj_","",.)
df_tumor <- as.data.frame((matrix(nrow=length(cm_files), ncol=3)))
colnames(df_tumor) <- c("samplename","Tumor_purity","Normal_purity")

df_tumor$samplename <- samplename
for(i in 1:length(cm_files)){
  conf_file <- read.table(cm_files[i], header = T, row.names = 1)
  conf_file$total_row <- rowSums(conf_file)
  conf_file["total_col",] <- colSums(conf_file)
  df_tumor[i,2] <- conf_file["Tumor","Tumor"]/conf_file["total_col","Tumor"]
  df_tumor[i,3] <- conf_file["normal","Normal"]/conf_file["total_col","Normal"]
}

write.table(df_tumor,"/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF8539/CNV/markerannot_P_Q_arm_2/Tumor_purity.txt", 
            sep ="\t", row.names = F, col.names = T, quote=F)

# combined_cnvscore$cellnames <- gsub("X","",rownames(combined_cnvscore))
obj@meta.data$cnv_score <- combined_cnvscore[match(rownames(obj@meta.data), combined_cnvscore$cellnames),"cnv_score"]
obj@meta.data$infercnv_tag <- combined_cnvscore[match(rownames(obj@meta.data), combined_cnvscore$cellnames),"type"]

p <- FeaturePlot(obj, features = "cnv_score", raster = TRUE) + scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,1500))
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"featureplot/infercnv_score_1500.pdf",sep=""))
p
dev.off()

# currently we donot have SF7994 we will subset it
cellnames <- rownames(obj@meta.data[grep("UCSFSF7994R_11_F_R_Run2", obj@meta.data$orig.ident, invert=TRUE),])
objsubset <- subset(obj, cells = cellnames)

p <- FeaturePlot(objsubset, features = "cnv_score", raster = TRUE) + scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,1500))
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"featureplot/infercnv_score_removed_SF7994_1500.pdf",sep=""))
p
dev.off()

p <- VlnPlot(objsubset, features = "cnv_score", group.by = "seurat_clusters", pt.size = 0, y.max = 2500) + geom_boxplot()
pdf(paste(savedir,"vlnplot/infercnv_score_removed_SF7994_2500.pdf",sep=""))
p
dev.off()

p <- VlnPlot(objsubset, features = "cnv_score", group.by = "seurat_clusters", pt.size = 0, y.max = 2500)
pdf(paste(savedir,"vlnplot/infercnv_score_removed_SF7994_all_noboxplot_2500.pdf",sep=""), width = 10, height = 8.5)
p
dev.off()

### Subsetting cluster 8 and 13
cellnames <- rownames(obj@meta.data[grep("^8$|^13$",obj@meta.data$seurat_clusters),])
req_cellnames <- grep("SF7994",cellnames,invert=TRUE,value=TRUE)
obj_subset <- subset(obj, cells = req_cellnames)

obj_subset@meta.data$cnv_score <- combined_cnvscore[match(rownames(obj_subset@meta.data), combined_cnvscore$cellnames),"cnv_score"]

p <- FeaturePlot(obj_subset, features = "cnv_score", raster = TRUE) + scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,2000))
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"featureplot/infercnv_score_removed_SF7994_1500_clus_8_13.pdf",sep=""))
p
dev.off()

tumor_cellname <- obj_subset@reductions$umap@cell.embeddings[obj_subset@reductions$umap@cell.embeddings[,"umap_1"] < 2.5,] %>% rownames()
normal_cellname <- obj_subset@reductions$umap@cell.embeddings[obj_subset@reductions$umap@cell.embeddings[,"umap_1"] > 2.5,] %>% rownames()

tumor_df <- data.frame(tumor_cellname)
tumor_df$type <- "tumor"

obj_subset@meta.data$tumor_normal <- tumor_df[match(rownames(obj_subset@meta.data),tumor_cellname),"type"]
obj_subset@meta.data$tumor_normal[is.na(obj_subset@meta.data$tumor_normal)] <- "normal"
obj_subset@meta.data$tumor_normal_clus <- paste(obj_subset@meta.data$tumor_normal,obj_subset@meta.data$seurat_clusters,sep="_")

p <- DimPlot(obj_subset, group.by = "tumor_normal_clus")
pdf(paste(savedir,"UMAP/tumor_normal_clus_8_13.pdf",sep=""))
p
dev.off()

# p <- DimPlot(objsubset, group.by = "seurat_clusters2")
# pdf(paste(savedir,"UMAP/obj_Seuratclusters.pdf",sep=""))
# p
# dev.off()

p <- VlnPlot(obj_subset, features = "cnv_score", group.by = "tumor_normal_clus", pt.size = 0)
pdf(paste(savedir,"vlnplot/objsubset_clus_8_13_all_noboxplot.pdf",sep=""))
p
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","GRIN2C")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","KIF21B","DISP3","BRINP1")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]

genes <- c(GCP,GNearly, GNlate, Stem_cells, nonmalignant_cells)

DefaultAssay(obj_subset) <- "RNA"
p <- DotPlot(obj_subset, feature = genes, group.by = "tumor_normal_clus", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  + scale_size(range = c(3, 6)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/clus_8_13_marker_genes.pdf", height = 20, width =9)
p
dev.off()


#### Now we are changing the celltype to different seurat clusters
obj_subset@meta.data$cellnames <- rownames(obj_subset@meta.data)
clus_8_13_cellname <- select(obj_subset@meta.data, c(tumor_normal_clus,cellnames))

write.table(clus_8_13_cellname, 
            "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/clus_8_13_cellnames.txt", 
            quote=F, 
            row.names = T, col.names = T, sep ="\t")

clus_8_13_cellname <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/clus_8_13_cellnames.txt")

objsubset@meta.data$seurat_clusters2 <- paste("clus",objsubset@meta.data$seurat_clusters,sep="_")
obj_index <- match(clus_8_13_cellname$cellnames,rownames(objsubset@meta.data))
objsubset@meta.data$seurat_clusters2[obj_index] <- clus_8_13_cellname$tumor_normal_clus

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","GRIN2C")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","KIF21B","DISP3","BRINP1")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]

genes <- c(GCP,GNearly, GNlate, Neuronal, Stem_cells, nonmalignant_cells)

DefaultAssay(objsubset) <- "MAGIC_RNA"
p <- DotPlot(objsubset, feature = genes, group.by = "seurat_clusters2", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_remove_SF7994_8_13_split_marker_genes.pdf", height = 20, width =9)
p
dev.off()

p <- VlnPlot(objsubset, features = "cnv_score", group.by = "seurat_clusters2", pt.size = 0, y.max = 2500) + geom_boxplot()
pdf(paste(savedir,"vlnplot/infercnv_score_removed_SF7994_clus2_2500.pdf",sep=""))
p
dev.off()

p <- VlnPlot(objsubset, features = "cnv_score", group.by = "seurat_clusters2", pt.size = 0, y.max = 2500)
pdf(paste(savedir,"vlnplot/infercnv_score_removed_SF7994_clus2_all_noboxplot_2500.pdf",sep=""), width = 10, height = 8.5)
p
dev.off()

p <- VlnPlot(objsubset, features = "cnv_score", group.by = "seurat_clusters2", pt.size = 0)
pdf(paste(savedir,"vlnplot/infercnv_score_removed_SF7994_clus2_all_noboxplot_all.pdf",sep=""), width = 10, height = 8)
p
dev.off()

### Based on the CNV score we can mark cells
# grep("clus_12|clus_15|clus_16|clus_17|clus_18|clus_19|clus_20|normal_13|tumor_8|tumor_8", objsubset@meta.data$seurat_clusters2)
objsubset@meta.data$CNV_tags <- "tumor"
objsubset@meta.data[grep("clus_12|clus_15|clus_16|clus_17|clus_18|clus_19|clus_20|normal_13|tumor_8|tumor_8", objsubset@meta.data$seurat_clusters2),"CNV_tags"] <- "normal"

p <- DimPlot(objsubset, group.by = "CNV_tags", cols = c("grey","red"))
pdf(paste(savedir,"UMAP/CNV_cancer.pdf",sep=""))
p
dev.off()

CT_table <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/seurat_cluster2_CT.txt", header=TRUE)

patterns = CT_table$clus
replacements = CT_table$celltypes

objsubset@meta.data$celltypes <- objsubset@meta.data$seurat_clusters2
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  objsubset@meta.data$celltypes <- str_replace_all(objsubset@meta.data$celltypes, pattern, replacements[i])
}

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/celltype_2.pdf")
DimPlot(obj,group.by = "celltypes",
cols = c("GCP"='#ff7f0e',"GCP_cycling"='#1f77b4',"GNearly"='#aa40fc',"GN_late"='#8c564b',"Pericytes"='#aec7e8',
                          "Astrocytes"='#279e68',"Immune"='#e377c2',"Oligo"='#b5bd61',"OPC"='#17becf',"Endo"='#d62728',"Neuronal"="#ffbb78"))
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","GRIN2C")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","KIF21B","DISP3","BRINP1")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]

genes <- c(GCP,GNearly, GNlate, Neuronal, Stem_cells, nonmalignant_cells)

DefaultAssay(objsubset) <- "MAGIC_RNA"
p <- DotPlot(objsubset, feature = genes, group.by = "celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_remove_SF7994_split_marker_genes.pdf", height = 24, width =9)
p
dev.off()

for(i in 1:22){
  print(paste("chr",i,"p chr",i,"q",sep = ""))
}

### generating CNV heatmap
infercnv <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/result/infercnv_table", header = TRUE, row.names = 1, sep = "\t")
infercnv2 <- infercnv %>%
  select(where(~ !all(is.na(.))))

infercnv2[is.na(infercnv2)] <- 0

library(ComplexHeatmap)

p <- Heatmap(t(infercnv2),  
        name = "Gene Number", 
        col = colorRampPalette(c("cornflowerblue","white", "lightcoral"))(50),
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal"))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/result/heatmap.pdf")
p
dev.off()

# Add Module Score UCell
# Since the object is huge so decreasing it to 50k
malignant_cells <- rownames(obj_impute@meta.data[grep(c("12|13|15|19|17|20"),obj_impute@meta.data$seurat_clusters,invert=TRUE),]) ## removing the nonmalignant cells
set.seed(123)
cells_selected <- sample(rownames(obj_impute@meta.data),60000)
malignant_cells_60k <- malignant_cells[match(cells_selected,malignant_cells, nomatch=0)]
obj_subset <- subset(obj_impute, cells = malignant_cells_60k)
obj <- obj_subset
obj <- AddModuleScore_UCell(obj, features = markers, ncores = cores, slot="data", assay = assay)

filename_UCell <- colnames(obj@meta.data)[29:62]

# p1 <- VlnPlot(obj, features = filename, group.by = "seurat_clusters")
p3 <- VlnPlot(obj, features = filename_UCell, group.by = "seurat_clusters", pt.size = 0)
p4 <- VlnPlot(obj, features = filename_UCell, group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()

rm(plot_list)
plot_list <- list()
for(i in 1:length(filename_UCell)){
  plot_list[[i]] <- FeaturePlot(obj, features = filename_UCell[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_Ucell.pdf",sep = ""), width = 14, height = 12)
print(p1)
dev.off()

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_Ucell.pdf",sep = ""), width = 14, height = 20)
print(p3)
dev.off()

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_bxp_Ucell.pdf",sep = ""), width = 14, height = 20)
print(p4)
dev.off()


dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/SHH_MB_het_gene_score_combined_Ucell.pdf",sep = ""))
plot_list
dev.off()

### Add Module score Seurat
  library(ArchR)
  library(UCell)
  library(ggplot2)
  # source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
  DefaultAssay(obj) <- "MAGIC_RNA"
  rm(markers)
  markers <- list()
  files <- list.files("/diazlab/data3/.abhinav/resources/gene_list/SHH_Gold_Taylor_2024_Nat_comm/",pattern=".txt", full.names = TRUE)
  filename <- gsub(".txt","",basename(files))
  gsub(".txt","",basename(files))
  for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[,1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset,rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i],Tcellsubset)
  }

obj <- AddModuleScore(obj, features = markers, slot="data")
colnames(obj@meta.data)[28:62] <- filename


p1 <- VlnPlot(obj, features = filename, group.by = "seurat_clusters")
p3 <- VlnPlot(obj, features = filename, group.by = "seurat_clusters", pt.size = 0)
p4 <- VlnPlot(obj, features = filename, group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()

rm(plot_list)
plot_list <- list()
for(i in 1:length(filename)){
  plot_list[[i]] <- FeaturePlot(obj, features = filename[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score.pdf",sep = ""), width = 14, height = 12)
print(p1)
dev.off()

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2.pdf",sep = ""), width = 14, height = 20)
print(p3)
dev.off()

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_bxp.pdf",sep = ""), width = 14, height = 20)
print(p4)
dev.off()


dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/SHH_MB_het_gene_score_combined.pdf",sep = ""))
plot_list
dev.off()

### Cerebellar development celltypes
  library(ArchR)
  library(UCell)
  library(ggplot2)
  # source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
  DefaultAssay(obj) <- "MAGIC_RNA"
  rm(markers)
  markers <- list()
  files <- list.files("/diazlab/data3/.abhinav/resources/gene_list/cerebellar_development/",pattern="genes.txt", full.names = TRUE)
  filename <- gsub(".txt","",basename(files))
  gsub(".txt","",basename(files))
  for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[,1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset,rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i],Tcellsubset)
  }

obj <- AddModuleScore(obj, features = markers, slot="data")
celltypes <- gsub("_genes","",filename)
colnames(obj@meta.data)[63:82] <- celltypes

# p1 <- VlnPlot(obj, features = celltypes, group.by = "seurat_clusters")
p3 <- VlnPlot(obj, features = celltypes, group.by = "seurat_clusters", pt.size = 0)
p4 <- VlnPlot(obj, features = celltypes, group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()

rm(plot_list)
plot_list <- list()
for(i in 1:length(celltypes)){
  plot_list[[i]] <- FeaturePlot(obj, features = celltypes[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_cerebellar_dev.pdf",sep = ""), width = 14, height = 20)
print(p3)
dev.off()

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_bxp_cerebellar_dev.pdf",sep = ""), width = 14, height = 20)
print(p4)
dev.off()


dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/SHH_MB_het_gene_score__cerebellar_dev_combined.pdf",sep = ""))
plot_list
dev.off()

#### Add Module Score Ucell
  library(ArchR)
  library(UCell)
  library(ggplot2)

  DefaultAssay(obj) <- "MAGIC_RNA"
  rm(markers)
  markers <- list()
  files <- list.files("/diazlab/data3/.abhinav/resources/gene_list/cerebellar_development/",pattern="genes.txt", full.names = TRUE)
  filename <- gsub(".txt","",basename(files))
  gsub(".txt","",basename(files))
  for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[,1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset,rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i],Tcellsubset)
  }

obj <- AddModuleScore_UCell(obj, features = markers, ncores = cores, slot="data", assay = assay)
# colnames(obj@meta.data)[28:62] <- filename

featurename <- paste(filename,"UCell",sep = "_")

# p1 <- VlnPlot(obj, features = featurename, group.by = "seurat_clusters")
p3 <- VlnPlot(obj, features = featurename, group.by = "seurat_clusters", pt.size = 0)
p4 <- VlnPlot(obj, features = featurename, group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()

rm(plot_list)
plot_list <- list()
for(i in 1:length(featurename)){
  plot_list[[i]] <- FeaturePlot(obj, features = featurename[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_cerebellar_dev_Ucell.pdf",sep = ""), width = 14, height = 20)
print(p3)
dev.off()

dir.create(paste(savedir,"vlnplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/SHH_MB_het_gene_score_2_bxp_cerebellar_dev_Ucell.pdf",sep = ""), width = 14, height = 20)
print(p4)
dev.off()


dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/SHH_MB_het_gene_score_combined_cerebellar_dev_Ucell.pdf",sep = ""))
plot_list
dev.off()

genes <- c("OTX2","OTX2-AS1","PRDM6", "KDM6A","CBFA2T2","CBFA2T3","SLIT2","WLS","TRPM3", "DPYD",
          "LMX1A", "EOMES", "PAX6", "MKI67")

genes <- c("TCF7L1","RFX2","SOX6","TCF7L2","OTX2","BRCA1","SUZ12","EZH2",
           "BARHL1","CTCF","EOMES","UNCX","LMX1A","SOX4","PPARA","NRF1","ZBTB37")

DefaultAssay(obj_impute3) <- "MAGIC_RNA"
genes <- c("HOXA3","HOXA10","HOXA9","OTX2","NTF3","UNC5D","GRIN2D","PRDM6","LINGO2")
genes <- c("RPGRIP1","SLC6A5","OTX2-AS1","SLC17A8","KHDRBS2","LEMD1")

genes <- c("OTX2","OTX2-AS1","PRDM6","UNC5D","TRPM3")

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- FeaturePlot(obj_remove_SF7994, features = genes[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"scdifferential/featureplot/RL_genes_3_removed_SF7994.pdf",sep = ""))
plot_list
dev.off()

pdf(paste(savedir,"dotplot/RL_genes.pdf",sep = ""),width = 9, height = 8)
DotPlot(obj, features = genes) + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP,GNearly, GNlate, Neuronal, Stem_cells, nonmalignant_cells)

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "seurat_clusters", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_marker_genes_0.8.pdf", height = 24, width =11)
p
dev.off()

genes <- c("SERPINF1","PATZ1","MIR503HG","AIMP1", "SPARC", "RGCC", "STAT1")
genes <- c("IGFBP7","MT-ND4", "EGLN3", "MT-ND5", "ANGPT2", "KCNJ8", "ITGA2", "MT-CO2", "HK2", "HSP90B1")
genes <- unique(c("EGLN3", "KCNJ8", "ITGA2","AIMP1", "MED1", "CLIC4", "FLT1", "ANGPT2", 
          "PDCD10", "APLNR", "CCL2", "TGFBI", "FGF1", "APLN",
          "ITGB1", "FLT1", "ANGPT2", "AKT3", "APLNR", "FGF1", "CD34", "HK2","KDR"))
genes <- c("PDIA3", "DNAJC19", "HSPA4", "PTGES3", "CCDC47", "CLPX", "HSP90B1", "CCT6A", 
           "DNAJB11", "CALR", "DNAJC10", "EIF2AK2", "CANX", "TBCE", "ANP32E", "FKBP9")

# genes <- c("NTRK2", "FLT1", "RYK", "IGFBP3", "ITGA1", "IGF2", "EIF2AK2", "NENF", "GPR37L1")
# genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/recurrent_high_genes.txt", header = FALSE)[,1]
DefaultAssay(snRNA) <- "MAGIC_RNA"
p <- DotPlot(snRNA, feature = genes, group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/Public_bulk_RNA_test_recurrent_ER_stress_unimputed.pdf")
p
dev.off()

genes <- c("IGFBP7","MT-ND4", "EGLN3", "MT-ND5", "ANGPT2", "KCNJ8", "ITGA2", "MT-CO2", "HK2", "HSP90B1")
# genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/recurrent_high_genes.txt", header = FALSE)[,1]
DefaultAssay(snRNA) <- "MAGIC_RNA"
p <- DotPlot(snRNA, feature = genes, group.by = "celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/Public_bulk_RNA_test_recurrent_lfc_gt_0.5_2_imputed.pdf")
p
dev.off()


#### BO reference mapping
BO_reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/BO_human_cerebellar/human.dev_cerebellum.glutamatergic.lineage_cells_SCT.rds")
DefaultAssay(obj_impute) <- "SCT"
DefaultAssay(BO_reference) <- "SCT"
anchors <- FindTransferAnchors(reference = BO_reference,
                                query = obj_impute,
                                normalization.method = "SCT",
                                dims = 1:50)
    
predictions.assay <- TransferData(anchorset = anchors,
                                  refdata = BO_reference$celltype,
                                  prediction.assay = TRUE,
                                  weight.reduction = obj_impute[["pca"]],
                                  dims = 1:30, k.weight = 50)
    
obj_impute[["BO_predictions"]] <- predictions.assay

celltypes <- as.vector(unique(BO_reference@meta.data$celltype))

DefaultAssay(obj_impute) <- "BO_predictions"
pdf(paste(savedir,"/featureplot/featureplot_reference_BO_ref_mapped.pdf",sep = ""), width = 12, height = 12)
FeaturePlot(obj_impute, features = celltypes)
dev.off()

prediction <- TransferData(anchorset = anchors,
                           refdata = BO_reference$celltype,
                           prediction.assay = FALSE,
                           weight.reduction = obj_impute[["pca"]],
                           dims = 1:30)

obj_impute <- AddMetaData(obj_impute, metadata = prediction)
p <- DimPlot(obj_impute,group.by = "predicted.id", reduction = "umap", label = TRUE)

dir.create(paste(savedir,"/UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"/UMAP/dimplot_celltypes_umap_BO_ref.pdf",sep = ""))
print(p)
dev.off()


#### Based on the Add Module Score we have subdivided the clustering

# Divide the cluster 5 and 14 into 2 below 0 and above 0 and then rename it
cells_5_14 <- rownames(obj_impute@meta.data[grep("^5$|^14$",obj_impute@meta.data$seurat_clusters),])
umap_coord <- as.matrix(obj_impute@reductions$umap@cell.embeddings)
umap_coord_df <- as.data.frame(unlist(umap_coord))
cells_less_than0 <- rownames(umap_coord_df[umap_coord_df$umap_1 < 0,])
cells_more_than0 <- rownames(umap_coord_df[umap_coord_df$umap_1 > 0,])
cells_5_14_lessthan0 <- cells_5_14[match(cells_less_than0,cells_5_14,nomatch=0)]
cells_5_14_morethan0 <- cells_5_14[match(cells_more_than0,cells_5_14,nomatch=0)]

obj_impute@meta.data$umap_dir_0 <- ""
obj_impute@meta.data[match(cells_5_14_lessthan0,rownames(obj_impute@meta.data)),"umap_dir_0"] <- "less_0"
obj_impute@meta.data[match(cells_5_14_morethan0,rownames(obj_impute@meta.data)),"umap_dir_0"] <- "more_0"
obj_impute@meta.data$seurat_clusters2 <- paste(obj_impute@meta.data$seurat_clusters, obj_impute@meta.data$umap_dir_0, sep = "")

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/cluster_res_0.8_divide.pdf")
DimPlot(obj_impute, group.by = "seurat_clusters2", label = TRUE)
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6")
HSP <- c("HSPH1","HSPA1B","HSPA1A","DNAJB1","HSPB1")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP, Stem_cells, HSP, GNearly, GNlate, Neuronal, nonmalignant_cells)

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "seurat_clusters2", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_marker_2_genes_0.8.pdf", height = 24, width =11)
p
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
RB <- c("RPL35","RPS5","RPL4","RPL26","RPL29","RPS17")
HSP <- c("HSPH1","HSPA1B","HSPA1A","DNAJB1","HSPB1")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP, Stem_cells,HSP,RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_celltype_marker_2_genes_0.8.pdf", height = 24, width =9)
p
dev.off()

GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
HSP <- c("HSPH1","HSPA1B","HSPA1A","DNAJB1","HSPB1")
RB <- c("RPL35","RPS5","RPL4","RPL26","RPL29","RPS17")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP, Stem_cells,HSP,RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_sub_celltype_marker_2_genes_0.8.pdf", height = 26, width =9)
p
dev.off()

### Making for publication
GCP <- c("PCNT","BOC","NTRK3","SHROOM2","SFRP1","PTCH2","TCERG1L","ATOH1","GLI2","GLI1","PTCH1")
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
HSP <- c("HSPH1","HSPA1B","HSPA1A","DNAJB1","HSPB1")
RB <- c("RPL35","RPS5","RPL4","RPL26","RPL29","RPS17")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP, Stem_cells,HSP,RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "sub_celltypes", assay = "MAGIC_RNA") + 
scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_sub_celltype_marker_2_genes_0.8.pdf", height = 26, width =9)
p
dev.off()

GCP <- c("ATOH1","GLI1")
Stem_cells <- c("MKI67","TOP2A")
HSP <- c("HSPH1","HSPA1B")
RB <- c("RPL35","RPS5")
GNearly <- c("NHLH1","STMN2")
GNlate <- c("NEUROD1","NEUROD2","RBFOX3","GRIN2C")
Astro <- c("GFAP","AQP4")
Endo <- c("CD34","VWF")
Immune <- c("CD14", "CD4")
OPC <- c("MYT1","VCAN")
Oligo <- c("MOBP" ,"PLP1")
Pericytes <- c("PDGFRB","RGS5")

genes <- c(GCP, Stem_cells,HSP,RB, GNearly, GNlate, Astro, Endo, Immune, OPC, Oligo, Pericytes)

DefaultAssay(snRNA) <- "MAGIC_RNA"
p <- DotPlot(snRNA, feature = genes, group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  
        coord_flip()  + scale_size(range = c(1, 10)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_sub_celltype_marker_2_paper_imputed.pdf", height = 12, width =9)
p
dev.off()



markers2 <- FindAllMarkers(obj_impute, group.by = "seurat_clusters2"
write.table(markers2, "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/cluster_marker_res_0.8_merged.txt", 
            sep = "\t", quote = F, col.names = T, row.names = T)

### based on the celltype and markers
Neuron_other <- c(21,"5more_0")
GCP_cycling_1 <- c(2,10)
GCP_HSP <- c(19)
GCP <- c(0,1,23)
GN_Premigratory <- c(3,4,15,"5less_0")
GN_migratory <- c(9,12, 26, "14less_0")
GCP_cycling_2 <- c(6,20)
GN_postmigratory <- c(7,8)
RL_like <- c(11)
GN_cycling <- c(13)
Immune <- c(16)
Pericyte <- c(17,18)
Oligo <- c(22)
OPC <- c(25)
Endo <- c(24)
Astro <- c("14more_0")

patterns <- c(Neuron_other, GCP_cycling_1, GCP_HSP, GCP, GN_Premigratory, 
              GN_migratory, GCP_cycling_2, GN_postmigratory, 
              RL_like,GN_cycling, Immune, Pericyte, Oligo, OPC,Endo, Astro)
replacements <- c(rep("Neuron_other",length(Neuron_other)),rep("GCP_cycling_1",length(GCP_cycling_1)),
                  rep("GCP_HSP",length(GCP_HSP)),rep("GCP",length(GCP)),rep("GN_Premigratory",length(GN_Premigratory)), rep("GN_migratory",length(GN_migratory)),
                  rep("GCP_cycling_2",length(GCP_cycling_2)),rep("GN_postmigratory",length(GN_postmigratory)),
                  rep("RL_like",length(RL_like)),rep("GN_cycling",length(GN_cycling)),rep("Immune",length(Immune)),rep("Pericyte",length(Pericyte)),
                  rep("Oligo",length(Oligo)),rep("OPC",length(OPC)),rep("Endo",length(Endo)),rep("Astro",length(Astro)))


obj_impute@meta.data$celltypes <- obj_impute@meta.data$seurat_clusters2
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  obj_impute@meta.data$sub_celltypes <- str_replace_all(obj_impute@meta.data$celltypes, pattern, replacements[i])
}

GCP_cycling <- c(2,10,6,20)
GN_cycling <- c(13)
GCP <- c(0,1,19,23) 
GN_early <- c(3,4,15,"5less_0")
GN_late <- c(7,8,9,12,26,"14less_0")
RL_like <- c(11)
Immune <- c(16)
Pericyte <- c(17,18)
Oligo <- c(22)
OPC <- c(25)
Endo <- c(24)
Astro <- c("14more_0")
Neuron_other <- c(21,"5more_0")

patterns <- c(GCP_cycling, GN_cycling, GCP, GN_early, GN_late, 
              RL_like, Immune, Pericyte, Oligo, OPC,Endo, Astro,Neuron_other)

replacements <- c(rep("GCP_cycling",length(GCP_cycling)),rep("GN_cycling",length(GN_cycling)),
                  rep("GCP",length(GCP)),rep("GN_early",length(GN_early)),rep("GN_late",length(GN_late)), 
                  rep("RL_like",length(RL_like)),rep("Immune",length(Immune)),rep("Pericyte",length(Pericyte)),
                  rep("Oligo",length(Oligo)),rep("OPC",length(OPC)),rep("Endo",length(Endo)),rep("Astro",length(Astro)),
                  rep("Neuron_other",length(Neuron_other)))


obj_impute@meta.data$celltypes <- obj_impute@meta.data$seurat_clusters2
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  obj_impute@meta.data$celltypes <- str_replace_all(obj_impute@meta.data$celltypes, pattern, replacements[i])
}

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/celltype_0.8_res.pdf")
DimPlot(obj_impute,group.by = "celltypes", label = TRUE,
cols = c("GCP_cycling"='#1a69a1',"GN_cycling"='#436957',"GN_early"='#5fe3a7',
         "GN_late"='#107044',"GCP" = "#6fb9ed","Pericyte"='#aa40fc',
         "Astro"='#8c564b',"Immune"='black',"Oligo"='#b5bd61',"OPC"='hotpink',
         "Endo"='#17becf',"Neuron_Other"="brown3","RL_like" = "cornsilk3"))
dev.off()


pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/sub_celltype_0.8_res.pdf")
DimPlot(obj_impute,group.by = "sub_celltypes", label = TRUE,
cols = c("GCP_cycling_1"='#1a69a1',"GCP_cycling_2"='#0a3d61',"GN_cycling"='#436957',"GN_Premigratory"='#5fe3a7',"GCP_HSP" = "#3b95d4",
         "GN_migratory"='#30b377',"GN_postmigratory" = "#107044","GCP" = "#6fb9ed","Pericyte"='#aa40fc',
         "Astro"='#8c564b',"Immune"='black',"Oligo"='#b5bd61',"OPC"='hotpink',
         "Endo"='#17becf',"Neuron_Other"="gray33","RL_like" = "cornsilk3","Neuron_other" = "brown3"))
dev.off()

RNA_integrated@meta.data$final_celltypes <- RNA_integrated@meta.data$sub_celltypes
RNA_integrated@meta.data$final_celltypes <- gsub("RL_like","GCP-GN_int",RNA_integrated@meta.data$final_celltypes) %>% 
                                            gsub("GCP_cycling_2","GCP_Rib",.) %>% 
                                            gsub("GN_cycling","GN_Rib",.) %>% 
                                            gsub("GCP_cycling_1","GCP_cyc",.) %>% 
                                            gsub("Neuron_other","Other",.) %>%
                                            gsub("Endo","Endothelial",.) %>%
                                            gsub("GN_postmigratory","GN_postmig",.) %>%
                                            gsub("GN_Premigratory","GN_premig",.) %>%
                                            gsub("GN_migratory","GN_mig",.)

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/sub_celltype_0.8_res_finalcelltypes.pdf")
DimPlot(RNA_integrated,group.by = "final_celltypes", label = TRUE,
cols = c("GCP_cyc"='darkslategray2',"GCP_Rib"='darkturquoise',"GN_Rib"='#436957',"GN_premig"='#5fe3a7',"GCP_HSP" = "deepskyblue4",
         "GN_mig"='#30b377',"GN_postmig" = "#107044","GCP" = "#6fb9ed","Pericyte"='#aa40fc',
         "Astro"='#8c564b',"Immune"='black',"Oligo"='#b5bd61',"OPC"='hotpink',
         "Endothelial"='#17becf',"GCP-GN_int" = "cyan4","Other" = "brown3"))
dev.off()

tiff("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/UMAP/sub_celltype_0.8_res_finalcelltypes.tiff", units="in", width=5.8, height=5, res=600)
DimPlot(RNA_integrated,group.by = "final_celltypes", label = TRUE,
cols = c("GCP_cyc"='darkslategray2',"GCP_Rib"='darkturquoise',"GN_Rib"='#436957',"GN_premig"='#5fe3a7',"GCP_HSP" = "deepskyblue4",
         "GN_mig"='#30b377',"GN_postmig" = "#107044","GCP" = "#6fb9ed","Pericyte"='#aa40fc',
         "Astro"='#8c564b',"Immune"='black',"Oligo"='#b5bd61',"OPC"='hotpink',
         "Endothelial"='#17becf',"GCP-GN_int" = "cyan4","Other" = "brown3"))
dev.off()


#### Longitudinal analysis
celltypes_ind <- table(obj_impute@meta.data$sample,obj_impute@meta.data$celltypes)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_individual_res_0.8.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_individual_res_0.8.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
  df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample","celltype","percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c("7316_4529","7316_5881","7316_278","7316_333","7316_737","7316_2118",
  "7316_931","7316_3023","7316_2978","7316_311","7316_1666",
  "7316_1676","DOD4182","SF10961","SF12930","SF7994",
  "SF8368","SF8539","SF9232") -> sample_levels

df_melted$sample <- gsub("-","_",df_melted$sample)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

df_melted$celltype <- factor(df_melted$celltype, levels = c("GCP_cycling","GN_cycling","GN_early","GN_late","GCP",
                                                            "Pericyte","Astro","Immune","Oligo","OPC","Endo","Neuron_other","RL_like"))

library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header= TRUE,sep="\t")

patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

### Patient ID
patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Patient.ID
# names(replacement) <- patterns
df_melted$Patient.ID <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$Patient.ID <- str_replace_all(df_melted$Patient.ID, pattern, replacements[i])
}

df_melted$Patient_timepts <- paste(df_melted$Patient.ID, df_melted$timepts,sep = "_")

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=Patient_timepts)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype_res_0.8_pts_timepoints.pdf",width =10, height = 7)
p
dev.off()

#### Performing the celltype distribution 
#####
df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",df_melted$sample),]
# df_melted_paired <- df_melted
# df_melted_paired$celltype <- factor(df_melted_paired$celltype, levels = c("GCP_cycling","GN_cycling","GN_early","GN_late","GCP",
#                                                             "Pericyte","Astro","Immune","Oligo","OPC","Endo","Neuron_other","RL_like"))

# df_melted_paired$Patient_timepts<- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts,sep ="_")

# p <- ggplot(df_melted_paired, aes(fill=celltype, y=percentage, x=Patient_timepts)) + 
#   geom_bar(position="fill", stat="identity") + 
#   # geom_text(aes(label = paste(percentage,"%")), position = position_stack(vjust = 0.5)) + 
# scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
#                                '#8c564b','black','#b5bd61','hotpink','#17becf',
#                                "brown3","cornsilk3")) + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype_timepoint_individuals_paired_res_0.8.pdf",width =7, height = 6)
# p
# dev.off()

df_melted_paired$Patient_Timepts <- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts, sep = "_")
df_wide <- dcast(df_melted_paired, celltype ~ Patient_Timepts, value.var = "percentage", fun.aggregate=mean)

write.table(df_wide, "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/df_melted_all_res_0.8_celltypes.txt",sep = "\t",
            row.names = F, col.names = T, quote=F)

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test$p <- (stat.test$p)/2

bxp <- ggboxplot(
  df_melted_paired, x = "celltype", y = "percentage", 
  color = "timepts", palette = c("#00AFBB", "#E7B800")
) 

bxp2 <- bxp +   geom_dotplot(
  aes(fill = timepts, color = timepts), trim = FALSE,
  binaxis='y', stackdir='center', dotsize = 0.15,
  position = position_dodge(0.8)
)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +   theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white"))

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"Table/t_test_unpaired_all_samples_bxplot_two_sided_celltype_res_0.8.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

#### subcelltypes 
#### Longitudinal analysis
celltypes_ind <- table(obj_impute@meta.data$sample,obj_impute@meta.data$sub_celltypes)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/sub_celltype_individual_res_0.8.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/sub_celltype_individual_res_0.8.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
  df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample","celltype","percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c("7316_4529","7316_5881","7316_278","7316_333","7316_737","7316_2118",
  "7316_931","7316_3023","7316_2978","7316_311","7316_1666",
  "7316_1676","DOD4182","SF10961","SF12930","SF7994",
  "SF8368","SF8539","SF9232") -> sample_levels

df_melted$sample <- gsub("-","_",df_melted$sample)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

df_melted$celltype <- factor(df_melted$celltype, levels = c("GCP_cycling_1","GCP_cycling_2","GN_cycling","GN_Premigratory", "GCP_HSP",
                                                            "GN_migratory","GN_postmigratory","GCP","Pericyte","Astro","Immune","Oligo","OPC",
                                                            "Endo","Neuron_other","RL_like"))

library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header= TRUE,sep="\t")

patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

### Patient ID
patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Patient.ID
# names(replacement) <- patterns
df_melted$Patient.ID <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$Patient.ID <- str_replace_all(df_melted$Patient.ID, pattern, replacements[i])
}

df_melted$Patient_timepts <- paste(df_melted$Patient.ID, df_melted$timepts,sep = "_")

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=Patient_timepts)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#0a3d61','#436957','#5fe3a7', "#3b95d4",
         '#30b377', "#107044", "#6fb9ed",'#aa40fc',
         '#8c564b','black','#b5bd61','hotpink',
         '#17becf', "brown3", "cornsilk3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_sub_celltype_res_0.8_pts_timepoints.pdf",width =10, height = 7)
p
dev.off()

#### Performing the celltype distribution 
#####
df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",df_melted$sample),]
df_melted_paired <- df_melted
# df_melted_paired$celltype <- factor(df_melted_paired$celltype, levels = c("GCP_cycling","GN_cycling","GN_early","GN_late","GCP",
#                                                             "Pericyte","Astro","Immune","Oligo","OPC","Endo","Neuron_other","RL_like"))

# df_melted_paired$Patient_timepts<- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts,sep ="_")

# p <- ggplot(df_melted_paired, aes(fill=celltype, y=percentage, x=Patient_timepts)) + 
#   geom_bar(position="fill", stat="identity") + 
#   # geom_text(aes(label = paste(percentage,"%")), position = position_stack(vjust = 0.5)) + 
# scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
#                                '#8c564b','black','#b5bd61','hotpink','#17becf',
#                                "brown3","cornsilk3")) + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//marker_celltype_timepoint_individuals_paired_res_0.8.pdf",width =7, height = 6)
# p
# dev.off()

df_melted_paired$Patient_Timepts <- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts, sep = "_")
df_wide <- dcast(df_melted_paired, celltype ~ Patient_Timepts, value.var = "percentage", fun.aggregate=mean)

write.table(df_wide, "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/df_melted_subcelltypes_all_samples_res_0.8.txt",sep = "\t",
            row.names = F, col.names = T, quote=F)

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test$p <- (stat.test$p)/2 ## for one sided unpaired

bxp <- ggboxplot(
  df_melted_paired, x = "celltype", y = "percentage", 
  color = "timepts", palette = c("#00AFBB", "#E7B800")
) 

bxp2 <- bxp +   geom_dotplot(
  aes(fill = timepts, color = timepts), trim = FALSE,
  binaxis='y', stackdir='center', dotsize = 0.15,
  position = position_dodge(0.8)
)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +   theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white"))

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"Table/t_test_unpaired_all_samples_bxplot_two_sided_celltype_res_0.8_subcelltypes.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

### Performing the differential for Primary and Recurrent
obj_impute@meta.data$sample2 <- gsub("-","_",obj_impute@meta.data$sample)
paired_cells <- row.names(obj_impute@meta.data[grep(c("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931"),obj_impute@meta.data$sample2),])
obj_paired <- subset(obj_impute, cells = paired_cells)
obj_paired <- NormalizeData(obj_paired) 

R_vs_P_all <- FindMarkers(obj_paired, ident.1 = "Recurrent", ident.2 = "Primary", group.by="tumor_type")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
dir.create(paste(savedir,"scdifferential"), showWarning=FALSE)
write.table(R_vs_P_all, paste(savedir,"scdifferential/R_vs_P_all_lognorm.txt",sep = ""),sep = "\t", col.names = T, quote=F)

obj_paired@meta.data$tumor_subcelltype <- paste(obj_paired@meta.data$tumor_type, obj_paired@meta.data$sub_celltypes, sep = "_")
sub_celltypes <- unique(obj_paired@meta.data$sub_celltypes)
dir.create(paste(savedir,"scdifferential/subcelltypes/",sep =""), showWarnings = FALSE)
for(i in 1:length(sub_celltypes)){
  R_vs_P_sub <- FindMarkers(obj_paired, ident.1 = paste("Recurrent_",sub_celltypes[i], sep=""), 
                            ident.2 = paste("Primary_",sub_celltypes[i],sep = ""), group.by="tumor_subcelltype")
  write.table(R_vs_P_sub, paste(savedir,"scdifferential/subcelltypes/",sub_celltypes[i],"_Recurrent_vs_Primary.txt",sep =""),sep = "\t", quote=F, col.names=T)
}

obj_paired@meta.data$tumor_celltype <- paste(obj_paired@meta.data$tumor_type, obj_paired@meta.data$celltypes, sep = "_")
celltypes <- unique(obj_paired@meta.data$celltypes)
dir.create(paste(savedir,"scdifferential/celltypes/",sep =""), showWarnings = FALSE)
for(i in 1:length(celltypes)){
  R_vs_P_cell <- FindMarkers(obj_paired, ident.1 = paste("Recurrent_",celltypes[i], sep=""), 
                            ident.2 = paste("Primary_",celltypes[i],sep = ""), group.by="tumor_celltype")
  write.table(R_vs_P_cell, paste(savedir,"scdifferential/celltypes/",celltypes[i],"_Recurrent_vs_Primary.txt",sep =""), sep = "\t", quote=F, col.names=T)
}

### Checking does gene expressed in all the samples

obj_paired@meta.data$patient_tumortype <- paste(obj_paired@meta.data$patient, obj_paired@meta.data$tumor_type,sep ="_")

rm(plot_list)
plot_list <- list()
genes <- c("NID1","COL4A1","COL11A1","COL1A1","COL1A2","COL3A1","DCN","LAMA1","LAMA2","FBN2", "SPARC","TIMP3","COL21A1","FREM3")
genes <- c("FBN1","JUN","FOS","EGFR","FGFR1","TGFBR3","PDGFD","ID1","SOX5","CUX2","BCL11B","UNC5D","GRIN2B")
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_paired,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"scdifferential/vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/vlnplot/recurrent_boxplot_public_bulk.pdf",sep = ""))
plot_list
dev.off()

DefaultAssay(obj_paired) <- "MAGIC_RNA"
genes <- c("NID1","COL4A1","COL11A1","COL1A1","COL1A2","COL3A1","DCN","LAMA1","LAMA2","FBN2", "SPARC","TIMP3","COL21A1","FREM3",
           "FBN1","JUN","FOS","EGFR","FGFR1","TGFBR3","PDGFD","ID1","SOX5","CUX2","BCL11B","UNC5D","GRIN2B")

# genes <- c("GRIN2B")
rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- FeaturePlot(obj_paired, features = genes[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"scdifferential/featureplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/featureplot/primary_genes.pdf",sep = ""))
plot_list
dev.off()

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/dotplot/primary_genes.pdf",sep = ""), width = 9, height = 8)
p
dev.off()

### While Performing differential getting genes mostly from the non malignant celltypes 
### Removing the non malignant genes 
obj_paired@meta.data[grep("GN_early|GN_cycling|GCP|GN_late|GCP_cycling|RL_like",obj_paired@meta.data$celltypes),] %>% rownames() -> paired_malignant_cells
obj_paired_malignant <- subset(obj_paired, cells = paired_malignant_cells)
DefaultAssay(obj_paired_malignant) <- "RNA"
obj_paired_malignant <- NormalizeData(obj_paired_malignant)
obj_paired_malignant <- ScaleData(obj_paired_malignant)

R_vs_P_mal <- FindMarkers(obj_paired_malignant, ident.1 = "Recurrent", ident.2 = "Primary", group.by="tumor_type")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
dir.create(paste(savedir,"scdifferential"), showWarning=FALSE)
write.table(R_vs_P_mal, paste(savedir,"scdifferential/R_vs_P_only_malignant_lognorm.txt",sep = ""),sep = "\t", col.names = T, quote=F)

DefaultAssay(obj_paired_malignant) <- "MAGIC_RNA"
genes <- c("NID1","COL4A1","COL11A1","COL1A1","COL1A2","COL3A1","DCN","LAMA1","LAMA2","FBN2", "SPARC","TIMP3","COL21A1","FREM3",
           "FBN1","JUN","FOS","EGFR","FGFR1","TGFBR3","PDGFD","ID1","SOX5","CUX2","BCL11B","UNC5D","GRIN2B")

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- FeaturePlot(obj_paired_malignant, features = genes[i], reduction = "umap", raster = TRUE) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"scdifferential/featureplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/featureplot/primary_mal_genes_2.pdf",sep = ""))
plot_list
dev.off()

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_paired_malignant, feature = genes, group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/dotplot/primary_genes_malignant_collagen_subcelltypes_RNA.pdf",sep = ""), width = 9, height = 8)
p
dev.off()

DefaultAssay(obj_paired_malignant) <- "MAGIC_RNA"

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(snRNA,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"scdifferential/vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/vlnplot/recurrent_genes_bulk_protein_folding.pdf",sep = ""))
plot_list
dev.off()

### Calcium Mediated Signaling
genes <- unique(c("RYR3","DMD","ITPR2","RYR2","GRIN2B","SLC8B1","PLCG2",
           "KCNC2","PRKG1","PDE3A","CHRM3","KCNC2","PRKD1","SCN3B","OPRM1"))

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_paired_malignant,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"scdifferential/vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/vlnplot/primary_genes_calcium_malignant_only.pdf",sep = ""))
plot_list
dev.off()

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_paired_malignant, feature = genes, group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/dotplot/primary_genes_calcium_malignant_subcelltypes.pdf",sep = ""), width = 9, height = 8)
p
dev.off()

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- FeaturePlot(obj_paired_malignant, features = genes[i], reduction = "umap", raster = TRUE) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"scdifferential/featureplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/featureplot/primary_mal_calcium_genes_2.pdf",sep = ""))
plot_list
dev.off()

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_paired_malignant, feature = genes, group.by = "patient_tumortype", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/dotplot/primary_genes_calcium_malignant_individuals.pdf",sep = ""), width = 9, height = 8)
p
dev.off()


# Recurrent Genes 
genes <- c("TUBD1","PTRH2","HEATR6","BCAS3","USP32","APPBP2","PPM1D","RPS5KB1","VMP1","CLTC") ### specific to one person

genes <- c("HMGCR","HMGCS1","INSIG1","INSIG2","MSMO1","ABCG1","ACSL3","ACSS2",
           "HSPA1A","HSPB1","HSPH1","HOXA10","HOXA9","DNAJB1","DNAJC12","SOX13",
           "CCND2","FOXO1","LEF1","VEGFA")

genes <- grep("HOXA",rownames(obj_paired_malignant),value=TRUE)

genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/cBioportal/PT_7WYPEC3Q_7316_5881_Recurrent/CNVs_genes.txt")[,1]
genes <- grep(paste0("^",genes,"$",collapse = "|"),row.names(obj_paired_malignant),value=TRUE)
genes <- c("KRIT1","ANKIB1","GATAD1","AKAP9","MTERF1","CYP51A1","PEX1")

genes <- c("AIF1")

genes <- unique(c("FLT1","KDR","FLT4"))
genes <- c("EDN1","SPARC","SLC16A1","ANGPT2","FLT1","KCNJ8","ITGA1","ELTD1",
          "SLC7A11","IFTI1","APLN","SLCO1A2","GPR116","APLNR","PLS3","IGFBP1",
          "FKBP9","CD34","GBP3")

genes <- c("SPP1","POSTN","IL16","CD4","CD163")

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"dotplot/SPP1_celltypes_MAGIC_RNA.pdf",sep = ""))
p
dev.off()

p <- DotPlot(obj_impute, feature = genes, group.by = "patient_tumortype", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/dotplot/AIF1_celltypes_individuals_MAGIC_RNA.pdf",sep = ""), width = 15, height = 8)
p
dev.off()

# DefaultAssay(obj_impute) <- "MAGIC_RNA"
genes <- c("PROM1","FUT4","SOX2","NGFR","NES")
p <- DotPlot(obj_impute, feature = genes, group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"dotplot/SPP1_subcelltypes_RNA.pdf",sep = ""), width = 9, height = 8)
p
dev.off()

p <- DotPlot(obj_malignant, feature = genes, group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create(paste(savedir,"scdifferential/dotplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/dotplot/RL_like_malignant_celltypes_RNA.pdf",sep = ""), width = 9, height = 8)
p
dev.off()


# rm(plot_list)
# plot_list <- list()
# for(i in 1:length(genes)){
#   plot_list[[i]] <- FeaturePlot(obj_paired_malignant, features = genes[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
# }

# dir.create(paste(savedir,"scdifferential/featureplot",sep =""), showWarnings = FALSE)
# pdf(paste(savedir,"scdifferential/featureplot/PT_7WYPEC3Q_7316_5881_Recurrent_genes.pdf",sep = ""))
# plot_list
# dev.off()

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_impute,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/VEGFA_genes.pdf",sep = ""))
plot_list
dev.off()

### Since for one sample most of them are coming in the RL like cells
RL_cellnames <- rownames(obj_impute@meta.data[grep("RL_like",obj_impute@meta.data$celltypes),])
obj_RL <- subset(obj_impute, cells = RL_cellnames)
RL_mat_transpose <- t(obj_RL@assays$MAGIC_RNA@data)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scRNA_correlation.R")
gene <- c("OTX2","OTX2-AS1","UNC5D","LINGO2")
for(i in 1:length(gene)){
  corre <- scpearson_correlation(RL_mat_transpose,gene[i],"/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/RL_like_clus_11/",gene[i])
  assign(paste(gene[i],"_corr",sep =""),corre)
}

write.table(RL_mat_transpose, "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/RL_like_clus_11/RL_mat_transpose.txt",
            sep = "\t", row.names = T, col.names = T, quote=F)

malignant_transpose <- t(obj_subset@assays$MAGIC_RNA@data)
gene <- c("OTX2","OTX2-AS1","UNC5D","LINGO2")
for(i in 1:length(gene)){
  corre <- scpearson_correlation(malignant_transpose,gene[i],"/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/RL_like_clus_11/",paste(gene[i],"malignant",sep =""))
  assign(paste(gene[i],"_corr_mal",sep =""),corre)
}

write.table(RL_mat_transpose, "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/RL_like_clus_11/RL_mat_transpose.txt",
            sep = "\t", row.names = T, col.names = T, quote=F)

### Since we are noticing new CNVs in the recurrent individuals when comparing with the primary
### Since RL like cell are an artifact so we are removing those cells
obj_impute@meta.data$celltype_samples <- paste(obj_impute@meta.data$celltypes, obj_impute@meta.data$orig.ident, sep = "_")
obj_impute_remove <- subset(obj_impute, cells = cells)

DefaultAssay(obj_impute_remove) <- "RNA"
obj_impute_remove <- NormalizeData(obj_impute_remove)
RL_markers <- FindMarkers(obj_impute_remove, ident.1 = "RL_like", group.by = "sub_celltypes")

genes <- rownames(head(RL_markers[order(-RL_markers$avg_log2FC),],50))

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_impute_remove,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"scdifferential/vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/vlnplot/RL_removed_like_marker_genes.pdf",sep = ""))
plot_list
dev.off()

pdf(paste(savedir,"scdifferential/dotplot/OTX2_MAGIC_RNA.pdf",sep = ""))
DotPlot(obj_impute, c("OTX2-AS1","OTX2","PRDM6","ITGA3","TRPM3"),group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste(savedir,"scdifferential/dotplot/OTX2_remove_SF7994M_AGIC_RNA.pdf",sep = ""))
DotPlot(obj_impute_remove, c("OTX2-AS1","OTX2","PRDM6","ITGA3","TRPM3"),group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

saveRDS(obj_impute_remove, "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed_remove_RL.RDS")


### Astrocytes turns out to be very different between the primary and recurrent
astrocells <- rownames(obj_impute@meta.data[grep("Astro",obj_impute@meta.data$celltypes),])
obj_astro <- subset(obj_impute, cells = astrocells)

genes <- c("SOX2","NES","FOXA1","GATA2","GATA3","MITF","SOX9","LAMA2")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"scdifferential/dotplot/stem_cells_astro_samples_3_orig.pdf",sep = ""))
DotPlot(obj_astro, genes,group.by = "orig.ident", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# genes <- c("LCN2","EOMES","MEIS2","SOX6","NR1H3","GATA2","GATA3","FOXA1","E2F1","E2F4","FOXA2","MYC","NR1H2","TCF12")

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- FeaturePlot(obj_astro, features = genes[i], reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir,"scdifferential/featureplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/featureplot/stem_cells_astro_samples_2.pdf",sep = ""))
plot_list
dev.off()

rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_astro,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"scdifferential/vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/vlnplot/Astro_like_markers_genes.pdf",sep = ""))
plot_list
dev.off()

### Making a plot for Primary and recurrent cells
library(dplyr)
library(ggplot2)
subcellfiles <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/subcelltypes/", pattern = "_Recurrent_vs_Primary.txt", full.names = TRUE)
cellname <- gsub("_Recurrent_vs_Primary.txt","",basename(subcellfiles))
rm(df_sub)
df_sub <- as.data.frame((matrix(nrow = length(cellname), ncol = 3)))
colnames(df_sub) <- c("cellname","Recurrent","Primary")
df_sub$cellname <- cellname

for(i in 1:length(subcellfiles)){
  subcell <- read.table(subcellfiles[i], header=TRUE, sep = "\t")
  df_sub[i,"Recurrent"] <- subcell[subcell$avg_log2FC > 1.5 & subcell$pct.1 > 0.05,] %>% nrow()
  df_sub[i,"Primary"] <- subcell[subcell$avg_log2FC < -1.5 & subcell$pct.2 > 0.05,] %>% nrow()
}


df_sub_melted <- melt(df_sub)
colnames(df_sub_melted) <- c("cellname","tumortype","Gene_number")
df_sub_melted$tumortype <- factor(df_sub_melted$tumortype,levels = c("Primary","Recurrent"))

df_sub_melted$cellname <- factor(df_sub_melted$cellname, levels = c("all","malignant","GCP","GCP_cycling_1","GCP_cycling_2","GCP_HSP","GN_cycling","GN_Premigratory",
                                          "GN_migratory","GN_postmigratory","Neuron_other","RL_like","Astro","Immune","Pericyte","Endo","OPC"))

p <- ggplot(df_sub_melted, aes(fill=tumortype, y=Gene_number, x=cellname)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(size =16)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/subcelltypes/cell_number.pdf",width =10, height = 8)
p
dev.off()

### Celltypes
library(dplyr)
library(ggplot2)
subcellfiles <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/celltypes/", pattern = "_Recurrent_vs_Primary.txt", full.names = TRUE)
cellname <- gsub("_Recurrent_vs_Primary.txt","",basename(subcellfiles))
rm(df_sub)
df_sub <- as.data.frame((matrix(nrow = length(cellname), ncol = 3)))
colnames(df_sub) <- c("cellname","Recurrent","Primary")
df_sub$cellname <- cellname

for(i in 1:length(subcellfiles)){
  subcell <- read.table(subcellfiles[i], header=TRUE, sep = "\t")
  df_sub[i,"Recurrent"] <- subcell[subcell$avg_log2FC > 1.5 & subcell$pct.1 > 0.05,] %>% nrow()
  df_sub[i,"Primary"] <- subcell[subcell$avg_log2FC < -1.5 & subcell$pct.2 > 0.05,] %>% nrow()
}


df_sub_melted <- melt(df_sub)
colnames(df_sub_melted) <- c("cellname","tumortype","Gene_number")
df_sub_melted$tumortype <- factor(df_sub_melted$tumortype,levels = c("Primary","Recurrent"))

df_sub_melted$cellname <- factor(df_sub_melted$cellname, levels = c("all","malignant","GCP","GCP_cycling","GN_cycling","GN_early",
                                          "GN_late","Neuron_other","RL_like","Astro","Immune","Pericyte","Endo","OPC"))

p <- ggplot(df_sub_melted, aes(fill=tumortype, y=Gene_number, x=cellname)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16), axis.text.y = element_text(size =16)) 

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/celltypes/cell_number.pdf",width =10, height = 8)
p
dev.off()

paired_cells <- row.names(obj_impute@meta.data[grep(c("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931"),obj_impute@meta.data$sample2),])
obj_paired <- subset(obj_impute, cells = paired_cells)
obj_paired <- NormalizeData(obj_paired) 

genes <- c("B4GALNT3","NINJ2","RAD52","WNT5B","ERC1","FBX14","CACNA1C","TSPAN9","PRMT8","CCND2-AS1","FGF23")
rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_paired,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir,"scdifferential/vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"scdifferential/vlnplot/recurrent_MYCN_high_markers_genes.pdf",sep = ""))
plot_list
dev.off()

obj_paired@meta.data[grep("GN_early|GN_cycling|GCP|GN_late|GCP_cycling|RL_like",obj_paired@meta.data$celltypes),] %>% rownames() -> paired_malignant_cells
obj_paired_malignant <- subset(obj_paired, cells = paired_malignant_cells)
DefaultAssay(obj_paired_malignant) <- "RNA"
obj_paired_malignant <- NormalizeData(obj_paired_malignant)
obj_paired_malignant <- ScaleData(obj_paired_malignant)

genes <- c("MYCN","MYCNUT","GACAT3")
pdf(paste(savedir,"scdifferential/dotplot/recurrent_MYC_RNA.pdf",sep = ""))
DotPlot(obj_paired_malignant, genes,group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Recurrent High genes
DefaultAssay(obj_paired_malignant) <- "RNA"

# primary and recurrent samples
patientid <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

for(i in 1:length(patientid)){
  group1 <- paste(patientid[i],"_Recurrent",sep = "")
  group2 <- paste(patientid[i],"_Primary",sep = "")
  markers <- FindMarkers(obj_paired_malignant, ident.1 = group1, ident.2 = group2, group.by = "patient_tumortype", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
  assign(paste(patientid[i],"markers_MAST",sep = "_"), markers)
  write.table(markers, paste(savedir,"Table/",patientid[i],"_recurrent_vs_primary_markers_MAST_nocutoff.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

### Combining the dataframe checking it has same
all(rownames(PT_7WYPEC3Q_markers_MAST) %in% rownames(PT_9S6WMQ92_markers_MAST)) # should be TRUE
all(rownames(PT_7WYPEC3Q_markers_MAST) %in% rownames(PT_H3WWDMW9_markers_MAST)) # should be TRUE
all(rownames(PT_7WYPEC3Q_markers_MAST) %in% rownames(PT_XA98HG1C_markers_MAST)) # should be TRUE

# Reorder data frames to have the same row order as PT_7WYPEC3Q_markers_MAST
PT_9S6WMQ92_markers_MAST_ordered <- PT_9S6WMQ92_markers_MAST[rownames(PT_7WYPEC3Q_markers_MAST), ]
PT_H3WWDMW9_markers_MAST_ordered <- PT_H3WWDMW9_markers_MAST[rownames(PT_7WYPEC3Q_markers_MAST), ]
PT_XA98HG1C_markers_MAST_ordered <- PT_XA98HG1C_markers_MAST[rownames(PT_7WYPEC3Q_markers_MAST), ]

## Sanity Check
all(rownames(PT_9S6WMQ92_markers_MAST_ordered) == rownames(PT_H3WWDMW9_markers_MAST_ordered))
all(rownames(PT_H3WWDMW9_markers_MAST_ordered) == rownames(PT_XA98HG1C_markers_MAST_ordered))
all(rownames(PT_7WYPEC3Q_markers_MAST) == rownames(PT_XA98HG1C_markers_MAST_ordered))

## Combined the dataframes
combined_df <- cbind("PT_7WYPEC3Q"=PT_7WYPEC3Q_markers_MAST,
                    "PT_H3WWDMW9"=PT_H3WWDMW9_markers_MAST_ordered, 
                    "PT_9S6WMQ92"=PT_9S6WMQ92_markers_MAST_ordered, 
                    "PT_XA98HG1C"=PT_XA98HG1C_markers_MAST_ordered)

# Calculating the Fisher Combined Probability using method H (pi) = −2ln( pi),
# Combining p-values using Fisher's method
combine_pvalues_fisher <- function(pvalues) {
  chisq_stat <- -2 * sum(log(pvalues))
  df <- 2 * length(pvalues)
  combined_pvalue <- pchisq(chisq_stat, df, lower.tail = FALSE)
  return(combined_pvalue)
}

combined_df$combined_adj_pvalue <- apply(combined_df[, selected_cols], 1, function(x) combine_pvalues_fisher(x))

# Define the columns to select for adjusted pvalue
selected_cols <- grep("p_val_adj",colnames(combined_df),value=TRUE)

# Use apply function to perform the calculation on each row
library(qvalue) # https://github.com/StoreyLab/qvalue
combined_df$combined_adj_pvalue <- apply(combined_df[, selected_cols], 1, function(x) combine_pvalues_fisher(x))

selected_cols <- grep(".avg_log2FC$",colnames(combined_df),value=TRUE)
combined_df$mean_logFC <- apply(combined_df[, selected_cols], 1, function(x) mean(x))
combined_df$median_logFC <- apply(combined_df[, selected_cols], 1, function(x) median(x))
combined_df$combined_adj_pvalue_BH_adjusted <- p.adjust(combined_df$combined_adj_pvalue, method = "BH")
combined_df$combined_adj_pvalue_Bonferroni_adjusted <- p.adjust(combined_df$combined_adj_pvalue, method = "bonferroni")
combined_df$combined_adj_pvalue_Storey_adjusted <- qvalue(combined_df$combined_adj_pvalue)$qvalues

markersname <- ls(pattern="_markers_MAST")

# extracting the primary and recurrent markers 
for(i in 1:length(markersname)){
  PTmarkers <- get(markersname[i])
  recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0.5,]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
  assign(paste(markersname[i],"_recurrent",sep = ""),recurrentmarkers)
  primarymarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC < -0.5,])
  assign(paste(markersname[i],"_primary",sep = ""),primarymarkers)
}

recurrentname <- ls(pattern = "_MAST_recurrent")

H <-list("PT_7WYPEC3Q_MAST_recurrent"=PT_7WYPEC3Q_markers_MAST_recurrent,
         "PT_9S6WMQ92_MAST_recurrent"=PT_9S6WMQ92_markers_MAST_recurrent,
         "PT_H3WWDMW9_MAST_recurrent"=PT_H3WWDMW9_markers_MAST_recurrent,
         "PT_XA98HG1C_MAST_recurrent"=PT_XA98HG1C_markers_MAST_recurrent)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
dir.create(paste(savedir,"Table/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"Table/MAST_recurrent_high_venny.pdf",sep = ""))
a
dev.off()

### Making an upset plot

n <- max(length(PT_7WYPEC3Q_markers_MAST_recurrent),length(PT_9S6WMQ92_markers_MAST_recurrent),
         length(PT_H3WWDMW9_markers_MAST_recurrent), length(PT_XA98HG1C_markers_MAST_recurrent))

length(PT_7WYPEC3Q_markers_MAST_recurrent) = n                    
length(PT_9S6WMQ92_markers_MAST_recurrent) = n
length(PT_H3WWDMW9_markers_MAST_recurrent) = n
length(PT_XA98HG1C_markers_MAST_recurrent) = n

df = as.data.frame(cbind((PT_7WYPEC3Q_markers_MAST_recurrent), (PT_9S6WMQ92_markers_MAST_recurrent), 
                        (PT_H3WWDMW9_markers_MAST_recurrent), (PT_XA98HG1C_markers_MAST_recurrent)))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

pdf(paste(savedir,"Table/MAST_recurrent_markers_upset_plot.pdf",sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list("PT_7WYPEC3Q"=PT_7WYPEC3Q_markers_MAST_recurrent,
                    "PT_9S6WMQ92"=PT_9S6WMQ92_markers_MAST_recurrent,
                    "PT_H3WWDMW9"=PT_H3WWDMW9_markers_MAST_recurrent,
                    "PT_XA98HG1C"=PT_XA98HG1C_markers_MAST_recurrent)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

df_int$int2 <- gsub("\\|","__",df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3),1]
morethan_3 <- df_int[(df_int$sample_num >= 3),1]

write.table(df_int, paste(savedir,"Table/MAST_recurrent_gene_intersect.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"Table/MAST_recurrent_intersected_gene_number.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir,"Table/MAST_recurrent_all_samples.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir,"Table/MAST_recurrent_morethan_3.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

#### Performing with the linent cutoff for logFC == 0
library(UpSetR)
library(dplyr)
filename <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table",
                        pattern = "_recurrent_vs_primary_markers_MAST_nocutoff.txt", full.names = TRUE) %>% 
                        grep("_Endo_|_HSP_",.,value=TRUE, invert = TRUE)
markersname <- gsub("_recurrent_vs_primary_markers_MAST_nocutoff.txt","",basename(filename))

# extracting the primary and recurrent markers 
for(i in 1:length(markersname)){
  PTmarkers <- read.table(filename[i], sep = "\t", header = TRUE)
  recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0,]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
  assign(paste(markersname[i],"_markers_MAST_recurrent",sep = ""),recurrentmarkers)
  primarymarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC < 0,])
  assign(paste(markersname[i],"_markers_MAST_primary",sep = ""),primarymarkers)
}

recurrentname <- ls(pattern = "_MAST_recurrent")

H <-list("PT_7WYPEC3Q_MAST_recurrent"=PT_7WYPEC3Q_markers_MAST_recurrent,
         "PT_9S6WMQ92_MAST_recurrent"=PT_9S6WMQ92_markers_MAST_recurrent,
         "PT_H3WWDMW9_MAST_recurrent"=PT_H3WWDMW9_markers_MAST_recurrent,
         "PT_XA98HG1C_MAST_recurrent"=PT_XA98HG1C_markers_MAST_recurrent)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
dir.create(paste(savedir,"Table/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"Table/MAST_recurrent_high_venny_lfc_0.pdf",sep = ""))
a
dev.off()

### Making an upset plot
n <- max(length(PT_7WYPEC3Q_markers_MAST_recurrent),length(PT_9S6WMQ92_markers_MAST_recurrent),
         length(PT_H3WWDMW9_markers_MAST_recurrent), length(PT_XA98HG1C_markers_MAST_recurrent))

length(PT_7WYPEC3Q_markers_MAST_recurrent) = n                    
length(PT_9S6WMQ92_markers_MAST_recurrent) = n
length(PT_H3WWDMW9_markers_MAST_recurrent) = n
length(PT_XA98HG1C_markers_MAST_recurrent) = n

df = as.data.frame(cbind((PT_7WYPEC3Q_markers_MAST_recurrent), (PT_9S6WMQ92_markers_MAST_recurrent), 
                        (PT_H3WWDMW9_markers_MAST_recurrent), (PT_XA98HG1C_markers_MAST_recurrent)))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

pdf(paste(savedir,"Table/MAST_recurrent_markers_upset_plot_lfc_0.pdf",sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list("PT_7WYPEC3Q"=PT_7WYPEC3Q_markers_MAST_recurrent,
                    "PT_9S6WMQ92"=PT_9S6WMQ92_markers_MAST_recurrent,
                    "PT_H3WWDMW9"=PT_H3WWDMW9_markers_MAST_recurrent,
                    "PT_XA98HG1C"=PT_XA98HG1C_markers_MAST_recurrent)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
library(stringr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

df_int$int2 <- gsub("\\|","__",df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3),1]
morethan_3 <- df_int[(df_int$sample_num >= 3),1]

write.table(df_int, paste(savedir,"Table/MAST_recurrent_gene_intersect_lfc_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"Table/MAST_recurrent_intersected_gene_number_lfc_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir,"Table/MAST_recurrent_all_samples_lfc_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir,"Table/MAST_recurrent_morethan_3_lfc_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")


### Performing the similar process with the primary genes
mast_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/",pattern = "_MAST_nocutoff.txt", full.names = TRUE)
markersname <- paste(gsub("_recurrent_vs_primary_markers_MAST_nocutoff.txt","",basename(mast_files)),"_markers_MAST",sep="")

for(i in 1:length(mast_files)){
  PTmarkers <- read.table(mast_files[i], header = TRUE, sep = "\t")
  recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0.5,]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
  assign(paste(markersname[i],"_recurrent",sep = ""),recurrentmarkers)
  primarymarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC < -0.5,])
  assign(paste(markersname[i],"_primary",sep = ""),primarymarkers)
}


primaryname <- ls(pattern = "_MAST_primary")

H <-list("PT_7WYPEC3Q_MAST_primary"=PT_7WYPEC3Q_markers_MAST_primary,
         "PT_9S6WMQ92_MAST_primary"=PT_9S6WMQ92_markers_MAST_primary,
         "PT_H3WWDMW9_MAST_primary"=PT_H3WWDMW9_markers_MAST_primary,
         "PT_XA98HG1C_MAST_primary"=PT_XA98HG1C_markers_MAST_primary)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
dir.create(paste(savedir,"Table/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"Table/MAST_primary_high_venny_lFC_gt_1.pdf",sep = ""))
a
dev.off()

### Making an upset plot
n <- max(length(PT_7WYPEC3Q_markers_MAST_primary),length(PT_9S6WMQ92_markers_MAST_primary),
         length(PT_H3WWDMW9_markers_MAST_primary), length(PT_XA98HG1C_markers_MAST_primary))

length(PT_7WYPEC3Q_markers_MAST_primary) = n                    
length(PT_9S6WMQ92_markers_MAST_primary) = n
length(PT_H3WWDMW9_markers_MAST_primary) = n
length(PT_XA98HG1C_markers_MAST_primary) = n

df = as.data.frame(cbind((PT_7WYPEC3Q_markers_MAST_primary), (PT_9S6WMQ92_markers_MAST_primary), 
                        (PT_H3WWDMW9_markers_MAST_primary), (PT_XA98HG1C_markers_MAST_primary)))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

pdf(paste(savedir,"Table/MAST_primary_markers_upset_plot_lfc_gt_1.pdf",sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list("PT_7WYPEC3Q"=PT_7WYPEC3Q_markers_MAST_primary,
                    "PT_9S6WMQ92"=PT_9S6WMQ92_markers_MAST_primary,
                    "PT_H3WWDMW9"=PT_H3WWDMW9_markers_MAST_primary,
                    "PT_XA98HG1C"=PT_XA98HG1C_markers_MAST_primary)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

df_int$int2 <- gsub("\\|","__",df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3),1]
morethan_3 <- df_int[(df_int$sample_num >= 3),1]

write.table(df_int, paste(savedir,"Table/MAST_primary_gene_intersect_lfc_gt_1.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"Table/MAST_primary_intersected_gene_number_lfc_gt_1.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir,"Table/MAST_primary_all_samples_lfc_gt_1.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir,"Table/MAST_primary_morethan_3_lfc_gt_1.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Performing with the primary gene less than 0
primaryname <- ls(pattern = "_MAST_primary")

H <-list("PT_7WYPEC3Q_MAST_primary"=PT_7WYPEC3Q_markers_MAST_primary,
         "PT_9S6WMQ92_MAST_primary"=PT_9S6WMQ92_markers_MAST_primary,
         "PT_H3WWDMW9_MAST_primary"=PT_H3WWDMW9_markers_MAST_primary,
         "PT_XA98HG1C_MAST_primary"=PT_XA98HG1C_markers_MAST_primary)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
dir.create(paste(savedir,"Table/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"Table/MAST_primary_high_venny_lFC_ls_0.pdf",sep = ""))
a
dev.off()

### Making an upset plot
n <- max(length(PT_7WYPEC3Q_markers_MAST_primary),length(PT_9S6WMQ92_markers_MAST_primary),
         length(PT_H3WWDMW9_markers_MAST_primary), length(PT_XA98HG1C_markers_MAST_primary))

length(PT_7WYPEC3Q_markers_MAST_primary) = n                    
length(PT_9S6WMQ92_markers_MAST_primary) = n
length(PT_H3WWDMW9_markers_MAST_primary) = n
length(PT_XA98HG1C_markers_MAST_primary) = n

df = as.data.frame(cbind((PT_7WYPEC3Q_markers_MAST_primary), (PT_9S6WMQ92_markers_MAST_primary), 
                        (PT_H3WWDMW9_markers_MAST_primary), (PT_XA98HG1C_markers_MAST_primary)))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

pdf(paste(savedir,"Table/MAST_primary_markers_upset_plot_lfc_ls_0.pdf",sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list("PT_7WYPEC3Q"=PT_7WYPEC3Q_markers_MAST_primary,
                    "PT_9S6WMQ92"=PT_9S6WMQ92_markers_MAST_primary,
                    "PT_H3WWDMW9"=PT_H3WWDMW9_markers_MAST_primary,
                    "PT_XA98HG1C"=PT_XA98HG1C_markers_MAST_primary)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

df_int$int2 <- gsub("\\|","__",df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3),1]
morethan_3 <- df_int[(df_int$sample_num >= 3),1]

write.table(df_int, paste(savedir,"Table/MAST_primary_gene_intersect_lfc_ls_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"Table/MAST_primary_intersected_gene_number_lfc_ls_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir,"Table/MAST_primary_all_samples_lfc_ls_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir,"Table/MAST_primary_morethan_3_lfc_ls_0.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")


### Identify the genes involved in the Ribosomal genes
genes <- c("EIF3E","EIF4EBP1","EIF2A","EIF2AK3","ATF4","ATF6","NFE2L2","EIF2S1")
genes <- unique(genes)
genes <- c("AKT1","FOXO1","RICTOR","CDKN1B","AKT3","VEGFA","HSPB1","CCND2","TP53","RPTOR","PTEN","FLT1")
# genes <- c("ATF4","EIF2A")
genes <- c("MAPK13","ELK1","DUSP8")
genes <- c("WDR49","NES","LAMA4","CCDC177")

genes <- c("DNAJA2","DNAJB1","DNAJC12","HSPA13","HSPA4L","HSPB1","HSPH1", "WDR49","NES","LAMA4","CCDC177")
# "EIF2AK3","EIF2AK4","EIF3E","EIF3A","BCL2","FOXO3","ATF4",
genes <- unique(c("PIK3CA","AKT3","FOXO1","RICTOR","VEGFA", "HSPB1","HIF1A","VHL"))

genes <- unique(c("FLT1","PTPRB","EGFL7","CLDN5","VWF","CDH5","CD93","TIE1","ERG","ABCB1","ENG","EGR1","CD34","PLVAP",
                "PDGFRB","RGS5","ACTA2","CDH6","COL1A1","COL1A2","COL6A2","COL3A1","EDNRA","DCN","BGN","NR2F2","PRKG1"))

genes <- c("CD133","NES","CD15","CD271")

genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/recurrent_important_genes", header=FALSE)[,1]


obj@meta.data$patient_tumortype <- paste(obj@meta.data$patient, obj@meta.data$tumor_type, sep = "_")
rm(plot_list)
plot_list <- list()
for(i in 1:length(genes)){
  plot_list[[i]] <- VlnPlot(obj_paired_malignant,  genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

genes <- c("PDK1","PIK3CA","VEGFA","FOXO1","AKT3","RICTOR")
genes <- c("FGFR1","EGFR","IGF1R","FLT1","KDR","ITGA5","ITGB1","ITGA6","ITGB3","ITGA2")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(obj_paired_malignant, genes[i],
        pt.size = 0,
        group.by = "patient_tumortype",
        assay = "MAGIC_RNA",
        cols = c("lightblue3", "lightcoral", "lightblue3", "lightcoral", "lightblue3", "lightcoral", "lightblue3", "lightcoral")
    ) + geom_boxplot()
}

dir.create(paste(savedir,"vlnplot",sep =""), showWarnings = FALSE)
pdf(paste(savedir,"vlnplot/Receptor_malignant.pdf",sep = ""))
plot_list
dev.off()

genes <- c("PPP1R14A","EDN1","SPARC","SLC16A1","ANGPT2","FLT1","KCNJ8","ITGA1","IGF2","ELTD1","ABCC9","HBA1",
          "SLC7A11","IFTI1","APLN","RGCC","SLCO1A2","GPR116","APLNR","PLS3","IGFBP1","FKBP9","CD34","GBP3")

pdf(paste(savedir,"dotplot/ER_stress_pathway_recurrent_MAGIC_RNA.pdf",sep = ""))
DotPlot(obj_paired_malignant, genes,group.by = "sub_celltypes", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + 
        coord_flip()  +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste(savedir,"dotplot/ER_stress_pathway_recurrent_RNA.pdf",sep = ""))
DotPlot(obj_paired_malignant, genes,group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + 
        coord_flip()  + scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf(paste(savedir,"dotplot/recurrent_RNA.pdf",sep = ""))
DotPlot(obj_paired_malignant, genes,group.by = "sub_celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + 
        coord_flip()  +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
genes <- c("FLT4", "KDR", "NRP1", "FLT1", "ESM1", "APLN", "APLNR", "CLDN5", "RAMP3", "NOTCH1", "JAG1", "DLL4")
pdf(paste(savedir,"dotplot/Endothelial_Angiogenesis_tip_cells_high_recurrent_RNA.pdf",sep = ""))
DotPlot(snRNA, genes, group.by = "celltypes", assay = "RNA") + 
        scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
        # scale_color_gradientn(colors = ArchRPalettes$solarExtra) + 
        coord_flip()  + scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

c("SP9","SP8","SP4","SP3","SP2","SP1",
  "KLF6","KLF5","KLF3","KLF2","KLF16","KLF15","KLF11",
  "EGR4","EGR3","EGR2","EGR1","FOXN3","FOXO3","FOXG1","FOXD1",
  "FOXP2","FOXA2","FOXA3","FOXK1","FOXP1","FOXC2","FOXK2","FOXP3",
  "FOXO6","FOXO4","FOXL1","FOXI1","FOXF2","HIF1A") -> genes


library(ArchR)
genes <- c("ANGTPL2","CD34","XBP1","GRN","SPARC","JUP","PRKCB","VASH2","PTPRM","NRG1",
"VASH1","THBS2","RHOB","AMOT","CDH5","ZC3H12A","NINJ1","RHOJ","CNTN1") 

genes <- c("XBP1","CREB3L1","LPCAT3","CLU","HERPUD1")

p <- DotPlot(snRNA, genes, assay = "RNA", group.by = "sub_celltypes") +
    # scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/Visium/public_data/PMID_37127652/downstream/"
dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/ER_Stress_RNA.pdf"), height = 9, width = 9)
print(p)
dev.off()


## Extracting out the GCP HSP and Endothelial cells 
## Since they are very interesting
GCP_HSP_cells  <- rownames(obj_paired@meta.data[grep("GCP_HSP",obj_paired@meta.data$sub_celltypes),])
obj_paired_HSP <- subset(obj_paired, cells = GCP_HSP_cells)

# Recurrent High genes
DefaultAssay(obj_paired_HSP) <- "RNA"
obj_paired_HSP <- NormalizeData(obj_paired_HSP)
obj_paired_HSP <- ScaleData(obj_paired_HSP)

patientid <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

for(i in 1:length(patientid)){
  group1 <- paste(patientid[i],"_Recurrent",sep = "")
  group2 <- paste(patientid[i],"_Primary",sep = "")
  markers <- FindMarkers(obj_paired_HSP, ident.1 = group1, ident.2 = group2, group.by = "patient_tumortype", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
  assign(paste(patientid[i],"markers_MAST",sep = "_"), markers)
  write.table(markers, paste(savedir,"Table/",patientid[i],"_HSP_recurrent_vs_primary_markers_MAST_nocutoff.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

### Performing the similar process with the primary genes
mast_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/",pattern = "_MAST_nocutoff.txt", full.names = TRUE)
markersname <- paste(gsub("_recurrent_vs_primary_markers_MAST_nocutoff.txt","",basename(mast_files)),"_markers_MAST",sep="")

# primary and recurrent samples
filename <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/", pattern = "_HSP_recurrent_vs_primary_markers_MAST_nocutoff.txt", full.names = TRUE)
samplename <- paste(gsub("_recurrent_vs_primary_markers_MAST_nocutoff.txt","",basename(filename)), "_markers_MAST",sep="")

for(i in 1:length(filename)){
  PTmarkers <- read.table(filename[i], header = TRUE, sep = "\t", row.names = 1)
  recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val < 0.05 & PTmarkers$avg_log2FC > 1,]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
  assign(paste(samplename[i],"_recurrent",sep = ""),recurrentmarkers)
  primarymarkers <- rownames(PTmarkers[PTmarkers$p_val < 0.05 & PTmarkers$avg_log2FC < -1,])
  assign(paste(samplename[i],"_primary",sep = ""),primarymarkers)
}


recurrentname <- ls(pattern = "HSP_markers_MAST_recurrent")

H <-list("PT_7WYPEC3Q_HSP_markers_MAST_recurrent"=PT_7WYPEC3Q_HSP_markers_MAST_recurrent,
         "PT_9S6WMQ92_HSP_markers_MAST_recurrent"=PT_9S6WMQ92_HSP_markers_MAST_recurrent,
         "PT_H3WWDMW9_HSP_markers_MAST_recurrent"=PT_H3WWDMW9_HSP_markers_MAST_recurrent,
         "PT_XA98HG1C_HSP_markers_MAST_recurrent"=PT_XA98HG1C_HSP_markers_MAST_recurrent)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
dir.create(paste(savedir,"Table/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"Table/HSP_MAST_recurrent_high_venny.pdf",sep = ""))
a
dev.off()

### Making an upset plot

n <- max(length(PT_7WYPEC3Q_HSP_markers_MAST_recurrent),length(PT_9S6WMQ92_HSP_markers_MAST_recurrent),
         length(PT_H3WWDMW9_HSP_markers_MAST_recurrent), length(PT_XA98HG1C_HSP_markers_MAST_recurrent))

length(PT_7WYPEC3Q_HSP_markers_MAST_recurrent) = n                    
length(PT_9S6WMQ92_HSP_markers_MAST_recurrent) = n
length(PT_H3WWDMW9_HSP_markers_MAST_recurrent) = n
length(PT_XA98HG1C_HSP_markers_MAST_recurrent) = n

df = as.data.frame(cbind((PT_7WYPEC3Q_HSP_markers_MAST_recurrent), (PT_9S6WMQ92_HSP_markers_MAST_recurrent), 
                        (PT_H3WWDMW9_HSP_markers_MAST_recurrent), (PT_XA98HG1C_HSP_markers_MAST_recurrent)))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

pdf(paste(savedir,"Table/HSP_MAST_recurrent_markers_upset_plot.pdf",sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list("PT_7WYPEC3Q"=PT_7WYPEC3Q_HSP_markers_MAST_recurrent,
                    "PT_9S6WMQ92"=PT_9S6WMQ92_HSP_markers_MAST_recurrent,
                    "PT_H3WWDMW9"=PT_H3WWDMW9_HSP_markers_MAST_recurrent,
                    "PT_XA98HG1C"=PT_XA98HG1C_HSP_markers_MAST_recurrent)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

df_int$int2 <- gsub("\\|","__",df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3),1]
morethan_2 <- df_int[(df_int$sample_num >= 2),1]

write.table(df_int, paste(savedir,"Table/HSP_MAST_recurrent_gene_intersect.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"Table/HSP_MAST_recurrent_intersected_gene_number.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir,"Table/HSP_MAST_recurrent_all_samples.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir,"Table/HSP_MAST_recurrent_morethan_3.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")



Endo_cells  <- rownames(obj_impute@meta.data[grep("Endo",obj_impute@meta.data$sub_celltypes),])
obj_Endo <- subset(obj_impute, cells = Endo_cells)

# Recurrent High genes
DefaultAssay(obj_Endo) <- "RNA"
obj_Endo <- NormalizeData(obj_Endo)
obj_Endo <- ScaleData(obj_Endo)


# primary and recurrent samples
patientid <- c("PT_7WYPEC3Q","PT_9S6WMQ92","PT_H3WWDMW9","PT_XA98HG1C")

for(i in 1:length(patientid)){
  group1 <- paste(patientid[i],"_Recurrent",sep = "")
  group2 <- paste(patientid[i],"_Primary",sep = "")
  markers <- FindMarkers(obj_Endo, ident.1 = group1, ident.2 = group2, group.by = "patient_tumortype", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
  assign(paste(patientid[i],"markers_MAST",sep = "_"), markers)
  write.table(markers, 
              paste(savedir,"Table/",patientid[i],"_Endo_recurrent_vs_primary_markers_MAST_nocutoff.txt",sep = ""), 
              quote = F, row.names = T, col.names = T, sep = "\t")
}

endo_paired_cells <- rownames(obj_Endo@meta.data[grep(paste(patientid, collapse="|"), obj_Endo@meta.data$patient),])
obj_Endo_paired <- subset(obj_Endo, cells = endo_paired_cells)

#### Performing GSEA enrichment for differentiated and undifferentiated
### GSEA for R vs P ####
library(fgsea)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(dplyr)

### For the cluster 7 genes we have to take the Cluster 5, 2, 3 and 2+3 for innate and adaptive
### combining cluster 2 and 3
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/"
dir.create(savedir, showWarnings = FALSE)
maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/"
files <- grep("P_vs_R", list.files(maindir, pattern = "DESeq_result.txt", full.names = TRUE, recursive = TRUE), invert = FALSE, value = TRUE)
Ag_name <- gsub(maindir, "", files) %>%
  gsub("_P_vs_R/DESeq_result.txt", "", .) %>%
  gsub("/", "", .)
Ag_name <- Ag_name[1:7]

# j=6
# add j= x number corresponding to number in Ag_name when only running a subset of comparisons, in that case skip the next line "for (....)" and keep running with "Ag_diff..."


# Pathways to be considered
# Load necessary libraries
library(fgsea)
library(msigdbr)
library(dplyr)
library(clusterProfiler)

# Download C5 Biological Process (BP) gene sets
msigdb_c5_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# Define the specific pathways of interest
pathways_of_interest <- c(
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_EIF2_ALPHA_PHOSPHORYLATION",
  "GOBP_NEGATIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_POSITIVE_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_POSITIVE_REGULATION_OF_TRANSLATION_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_REGULATION_OF_TRANSLATION_INITIATION_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_REGULATION_OF_TRANSLATION_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"
)

# Filter the gene sets for pathways of interest
filtered_gene_sets <- msigdb_c5_bp %>% filter(gs_name %in% pathways_of_interest)

# Convert the filtered gene sets to list format
gene_sets <- filtered_gene_sets %>% split(x = .$entrez_gene, f = .$gs_name)

R_vs_P <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/MAST_nocutoff_Fisher_combined_pvalue_adjusted.txt", header = TRUE, sep = "\t")

## GSEA
## extracting the genes required
R_vs_P_Entrez <- bitr(rownames(R_vs_P), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
R_vs_P_Entrez[duplicated(R_vs_P_Entrez$SYMBOL) # Identify the duplicates
| duplicated(R_vs_P_Entrez$SYMBOL, fromLast = TRUE), ]
# Ag_diff_req_Entrez_unique <- tops_all_Entrez[Ag_diff_req_Entrez$ENTREZID != 7006 & Ag_diff_req_Entrez$ENTREZID != 100124696,]
## Adding the fold Change
R_vs_P_Entrez$logFC <- R_vs_P[match(R_vs_P_Entrez$SYMBOL, rownames(R_vs_P)), "mean_logFC"]

write.table(R_vs_P_Entrez, 
            "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/GSEA_R_vs_P_mean.txt", 
            col.names = TRUE, sep = "\t", row.names = F, quote=F)

library(clusterProfiler)
geneList <- R_vs_P_Entrez$logFC ## logfc
names(geneList) <- as.character(R_vs_P_Entrez$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)

library(fgsea)
library(data.table)
library(ggplot2)

# Run fgsea
fgsea_results <- fgsea(
  pathways = gene_sets, stats = geneList,
  minSize = 8,
  eps = 0.0,
  maxSize = 1500
)

# View results
fgsea_results <- fgsea_results %>% arrange(pval)
print(fgsea_results)

# Convert the list column to a string
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ","))

write.table(as.data.frame(fgsea_results), 
            "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/ER_stress_GSEA.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)

pdf(paste(savedir, Ag_name[j], "_P_vs_R/GCP_GN_GSEA_table.pdf", sep = ""))
pdf(paste(savedir,Ag_name[j],"_",filename_pathway,"_all_O_vs_Y_table.pdf",sep = ""))
print(plotGseaTable(pathways[cell_name], geneList, fgseaRes, gseaParam = 0.5))
dev.off()

### It is better to make it separately
R_vs_P_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/", pattern = "_recurrent_vs_primary_markers_MAST_nocutoff.txt", full.names = TRUE, recursive = TRUE)
R_vs_P_files2 <- grep("_HSP_|_ENDO_", R_vs_P_files, invert = TRUE, value = TRUE, ignore.case = T)
R_vs_P_filesname <- gsub("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table//","",R_vs_P_files2) %>% gsub("_recurrent_vs_primary_markers_MAST_nocutoff.txt","",.)

for (j in 1:length(R_vs_P_filesname)) {
  Ag_diff <- read.table(R_vs_P_files2[j], header = TRUE, sep = "\t")

  ## extracting the genes required
  Ag_diff_req <- Ag_diff
  Ag_diff_req_Entrez <- bitr(rownames(Ag_diff_req), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  Ag_diff_req_Entrez[duplicated(Ag_diff_req_Entrez$SYMBOL) # Identify the duplicates
  | duplicated(Ag_diff_req_Entrez$SYMBOL, fromLast = TRUE), ]

  ## Adding the fold Change
  Ag_diff_req_Entrez$logFC <- Ag_diff_req[match(Ag_diff_req_Entrez$SYMBOL, rownames(Ag_diff_req)), "avg_log2FC"]

  ### Saving the file

  # write.csv(Ag_diff_req_Entrez, file = paste(savedir,Ag_name[j],"_all_O_vs_Y_diff_req_Entrez.csv",sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE)

  library(clusterProfiler)
  geneList <- Ag_diff_req_Entrez$logFC ## logfc
  names(geneList) <- as.character(Ag_diff_req_Entrez$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)

  library(fgsea)
  library(data.table)
  library(ggplot2)

  # pathways_combined <- list()
  rm(pathways)
  pathways <- list()
  T_cells <- c("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/ER_stress_entrez_list.txt", 
               "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/PI3K_AKT_signaling.txt")
  cell_name <- gsub("_entrez_list.txt|.txt", "", basename(T_cells))
  for (k in 1:length(T_cells)) {
    Tfh_entrez <- read.table(T_cells[k])[, 1]
    # Tfh_entrez <- bitr(tfh_genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
    # Tfh_entrez[, 2]
    #
    pathways[[cell_name[k]]] <- Tfh_entrez
    geneList

    fgseaRes <- fgsea(
      pathways = pathways,
      stats = geneList,
      minSize = 8,
      eps = 0.0,
      maxSize = 1500
    )

    ## If you get the warning it is due to the  genes in your preranked list having the same numerical/ranking value. For example, if you rank genes by log2 fold change,
    # genes with the same log2 fold change value will be ordered randomly once processed. You can use rank() function and a ties.method parameter to control how genes
    # with the same value are ranked:

    ### for this it
    pdf(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/GSEA/",R_vs_P_filesname[j],"_",cell_name[k], "_GSEA.pdf", sep = ""))
    print(plotEnrichment(pathways[[cell_name[k]]], geneList, gseaParam = 1) + labs(title = cell_name[k]))
    dev.off()
  }

   pdf(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/GSEA/",R_vs_P_filesname[j],"_",cell_name[k], "_GSEA_table.pdf", sep = ""))
  # pdf(paste(savedir,Ag_name[j],"_",filename_pathway,"_all_O_vs_Y_table.pdf",sep = ""))
  print(plotGseaTable(pathways[cell_name], geneList, fgseaRes, gseaParam = 0.5))
  dev.off()

  fgsea_df <- as.matrix(fgseaRes)
  write.table(fgsea_df, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/GSEA/",R_vs_P_filesname[j],"_",cell_name[k], "_GSEA_table.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

write.table(fgsea_df, "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/clus_all_P_vs_R/neuron_differentiation.txt", quote = F, row.names = T, col.names = T, sep = "\t")


## Cell Cyling 
cellcycling <- table(obj_impute@meta.data$Phase,obj_impute@meta.data$sub_celltypes)
ColumnSum <- colSums(cellcycling)
df_CC <- (t(cellcycling)/ColumnSum)*100
df_melted <- melt(df_CC)
colnames(df_melted) <- c("Celltypes","Phase", "Percentage")
p <- ggplot(df_melted, aes(fill=Phase, y=Percentage, x=Celltypes)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('indianred1','green2','steelblue2')) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//cellcycle_barplot.pdf",width =6, height = 5)
p
dev.off()

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table//Endo_phase_NES.pdf")
VlnPlot(obj_Endo,  "NES", pt.size = 0, group.by = "Phase", assay = "MAGIC_RNA", split.by = "tumor_type") + geom_boxplot()
dev.off()


### Cell-Cell Interactions
library(iTALK)
library(dplyr)
## significant ligand-receptor pairs between compare groups
# find DEGenes of regulatory T cells and NK cells between these 2 groups
data <- t(as.data.frame(obj_subset@assays$RNA@counts))
data_df <- as.data.frame(data)
stopifnot(all(rownames(obj_subset@meta.data) == rownames(data_df)))

data_df$cell_type <- obj_subset@meta.data$super_celltype
data_df$compare_group <- obj_subset@meta.data$tissue
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a',
                      "blue","green","yellow","red","brown3","brown2"),
                    names=unique(data_df$cell_type))


deg_t<-DEG(data_df %>% filter(cell_type=='T_cells'),method='Wilcox',contrast=c("Lesional skin","Non-lesional skin"))
deg_mps<-DEG(data_df %>% filter(cell_type=='MPs'),method='Wilcox',contrast=c("Lesional skin","Non-lesional skin"))
# find significant ligand-receptor pairs and do the plotting
par(mfrow=c(1,2))
rm(res)
res <- list()
dir.create(paste(savedir,"LRI",sep = ""), showWarnings = FALSE)
comm_list<-c('growth factor','other','cytokine','checkpoint')
for(comm_type in comm_list){
  res_cat<-FindLR(deg_t,deg_mps,datatype='DEG',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
  pdf(paste(savedir,"LRI/",comm_type,".pdf",sep = ""))
  #plot by ligand category
  if(nrow(res_cat)==0){
    next
  }else if(nrow(res_cat>=20)){
    LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])
  }else{
    LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)
  }
  dev.off()
  pdf(paste(savedir,"LRI/",comm_type,"_NetView.pdf",sep = ""))
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5) 
  dev.off()
  res[[comm_type]] <- res_cat
  # res<-rbind(res,res_cat)
}

res_df <- as.data.frame(do.call(rbind, res))
write.table(res_df, paste(savedir,"LRI/result_table.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

if(is.null(res)){
  print('No significant pairs found')
}else if(nrow(res)>=20){
  res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]
  pdf(paste(savedir,"LRI/",comm_type,"Netview_all.pdf",sep = ""))
  NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  dev.off()
  pdf(paste(savedir,"LRI/",comm_type,"_all.pdf",sep = ""))
  LRPlot(res[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:20],link.arr.width=res$cell_to_logFC[1:20])
  dev.off()
}else{
  pdf(paste(savedir,"LRI/",comm_type,"_NetView_all.pdf",sep = ""))
  NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  dev.off()
  pdf(paste(savedir,"LRI/",comm_type,"_all.pdf",sep = ""))
  LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
  dev.off()
}

### Cell Chat
library(CellChat)
Idents(obj_impute) <- obj_impute@meta.data$sub_celltypes
cellChat <- createCellChat(object = obj_impute, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
# Before users can employ CellChat to infer cell-cell communication, they need to set the ligand-receptor interaction database and identify over-expressed ligands or receptors.

# Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB v2 contains ~3,300 
# validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling. Compared to CellChatDB v1, CellChatDB v2 adds more than 1000 protein and non-protein interactions such as metabolic and synaptic signaling. It should be noted that for molecules that are not directly related to genes measured in scRNA-seq, CellChat v2 estimates the expression of ligands and receptors using those molecules’ key mediators or enzymes for potential communication mediated by non-proteins.

# CellChatDB v2 also adds additional functional annotations of ligand-receptor pairs, such as UniProtKB keywords (including biological process, molecular function, 
# functional class, disease, etc), subcellular location and relevance to neurotransmitter.
# Users can update CellChatDB by adding their own curated ligand-receptor pairs. Please check the tutorial on updating the ligand-receptor interaction database  CellChatDB.
# When analyzing human samples, use the database CellChatDB.human; when analyzing mouse samples, use the database CellChatDB.mouse. CellChatDB categorizes ligand-receptor pairs into different types, including “Secreted Signaling”, “ECM-Receptor”, “Cell-Cell Contact” and “Non-protein Signaling”. By default, the “Non-protein Signaling” are not used.

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/cellchatdb.pdf")
showDatabaseCategory(CellChatDB)
dev.off()


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellChat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis

# To infer the cell state-specific communications, CellChat identifies over-expressed ligands or receptors in one cell group and then identifies over-expressed 
# ligand-receptor interactions if either ligand or receptor are over-expressed.

# We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ 
# expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing 
# single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of
#  subunits of ligands/receptors. One might be concerned about the possible artifact introduced by this diffusion process, however, it will only introduce very weak 
#  communications. By default CellChat uses the raw data (i.e., object@data.signaling) instead of the projected data. To use the projected data, users should run the 
#  function projectData before running computeCommunProb, and then set raw.use = FALSE when running computeCommunProb.

# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

#> [1] 13.20763
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# To determine a proper value of trim, CellChat provides a function computeAveExpr, which can help to check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1). Therefore, if well-known signaling pathways in the studied biological process are not predicted, users can try truncatedMean with lower values of trim to change the method for calculating the average gene expression per cell grou
# computeAveExpr(cellchat, features = c("VEGFA","KDR","FLT1","FLT4"), type =  "truncatedMean", trim = 0.11)
options(future.globals.maxSize = 1000 * 1024^2)
cellChat <- computeCommunProb(cellChat, type = "truncatedMean", trim = 0.11)

# Filtering for minimum number of cells
cellchat <- filterCommunication(cellChat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred
# communications at the level of signaling pathways
df.net <- subsetCommunication(cellchat)

# gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))

# gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net <- subsetCommunication(cellchat, signaling = c("VEGF"))

# Infer the cell-cell communication at a signaling pathway level
# CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions 
# associated with each signaling pathway.
# NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
# CellChat calculates the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. Users can also 
# calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
cellchat <- aggregateNet(cellchat)

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/interaction_strenght.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/interaction_strenght_splitted.pdf", width = 12, height = 12)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("VEGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/VEGF_signaling.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()
# Circle plot
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/VEGF_signaling_circle.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
par(mfrow=c(1,1))
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/VEGF_signaling_chord.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

# Heatmap
par(mfrow=c(1,1))
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/VEGF_signaling_heatmap.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

# Do heatmap based on a single object
pairLR.VEGF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.VEGF # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/VEGF_VEGFR.pdf", width = 12, height = 15)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()

#> [[1]]
# Circle plot
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/VEGF_VEGFR1_cirle.pdf")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()


# Automatically save the plots of the all inferred network for quick exploration

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/"
setwd(savedir)
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  try({
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  pdf(paste0("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/cellcomm/",pathways.show.all[i],"_cirle.pdf"), width = 5, height = 5)
  netVisual_individual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle")
  dev.off()
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(savedir,pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  })
}

## Trying to identify the pathway in the recurrent samples
# Pathway activity inference from scRNA-seq using decoupleR
## We load the required packages
library(Seurat)
library(decoupleR)

# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)


# PROGENy model
# PROGENy is a comprehensive resource containing a curated collection of pathways and their target genes, with weights for each interaction. For this example we will use the human weights (other organisms are available) and we will use the top 500 responsive genes ranked by p-value. Here is a brief description of each pathway:

#     Androgen: involved in the growth and development of the male reproductive organs.
#     EGFR: regulates growth, survival, migration, apoptosis, proliferation, and differentiation in mammalian cells
#     Estrogen: promotes the growth and development of the female reproductive organs.
#     Hypoxia: promotes angiogenesis and metabolic reprogramming when O2 levels are low.
#     JAK-STAT: involved in immunity, cell division, cell death, and tumor formation.
#     MAPK: integrates external signals and promotes cell growth and proliferation.
#     NFkB: regulates immune response, cytokine production and cell survival.
#     p53: regulates cell cycle, apoptosis, DNA repair and tumor suppression.
#     PI3K: promotes growth and proliferation.
#     TGFb: involved in development, homeostasis, and repair of most tissues.
#     TNFa: mediates haematopoiesis, immune surveillance, tumour regression and protection from infection.
#     Trail: induces apoptosis.
#     VEGF: mediates angiogenesis, vascular permeability, and cell migration.
#     WNT: regulates organ morphogenesis during development and tissue repair.

net <- get_progeny(organism = 'human')

# Activity inference with Multivariate Linear Model (MLM)
# To infer pathway enrichment scores we will run the Multivariate Linear Model (mlm) method. For each sample in our dataset (mat), it fits a linear 
# model that predicts the observed gene expression based on all pathways’ Pathway-Gene interactions weights. Once fitted, the obtained t-values of 
# the slopes are the scores. If it is positive, we interpret that the pathway is active and if it is negative we interpret that it is inactive.

# Extract the normalized log-transformed counts
data_RNA <- as.matrix(snRNA@assays$RNA@data)
# mat <- data_RNA
# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
snRNA[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(snRNA) <- "pathwaysmlm"

# Scale the data
snRNA <- ScaleData(snRNA)
snRNA@assays$pathwaysmlm@data <- snRNA@assays$pathwaysmlm@scale.data

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/pathway/"
pdf(paste(savedir,"PI3K_dimplot.pdf"))
(FeaturePlot(snRNA, features = c("PI3K")) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('PI3k activity')

# Extract activities from object as a long dataframe
Idents(snRNA) <- snRNA@meta.data$patient_tumortype
df <- t(as.matrix(snRNA@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(snRNA)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
pdf(paste(savedir,"pathway_heatmap_patient_tumor_type.pdf",sep=""))
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()

### Subsetting only to malignant and paired samples
malignant_cells <- rownames(obj@meta.data[grep("Oligo|Endo|Astro|Pericyte|OPC|Neuron_other|Immune",obj@meta.data$sub_celltypes,invert=TRUE),])
obj_mal <- subset(obj, cells = malignant_cells)

paired_malignant_cells <- rownames(obj_mal@meta.data[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",obj_mal@meta.data$sample2),])
obj_mal_paired <- subset(obj_mal, cells = paired_malignant_cells)

# Change assay
DefaultAssay(obj_mal_paired) <- "pathwaysmlm"

# Scale the data
obj_mal_paired2 <- obj_mal_paired
obj_mal_paired <- ScaleData(obj_mal_paired)
obj_mal_paired@assays$pathwaysmlm@data <- obj_mal_paired@assays$pathwaysmlm@scale.data

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/pathway/"
pdf(paste(savedir,"PI3K_dimplot_mal_paired_2.pdf",sep=""))
(FeaturePlot(obj_mal_paired, features = c("PI3K")) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('PI3k activity')
dev.off()

# Extract activities from object as a long dataframe
Idents(obj_mal_paired) <- obj_mal_paired@meta.data$patient_tumortype
df <- t(as.matrix(obj_mal_paired@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(obj_mal_paired)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
top_acts_mat_ordered <- top_acts_mat[order(rownames(top_acts_mat)),]
pdf(paste(savedir,"pathway_heatmap_patient_tumor_type_mal_paired_ordered.pdf",sep=""))
pheatmap(top_acts_mat_ordered[c(1,3,5,7,2,4,6,8),], border_color = NA, color=my_color, breaks = my_breaks, cluster_rows=FALSE) 
dev.off()

#### Identifying the difference b/w GCP and GN in primary and recurrent
snRNA@meta.data$broad_CT <- snRNA@meta.data$celltypes

 patterns <- c("Oligo","GCP","Astro","Pericyte","OPC","Neuron_other","GN_early",
 "GCP_cycling","GN_late","GN_cycling","Immune","RL_like","Endo")

 replacements <- c("non-malign","GCP","non-malign","non-malign","non-malign","non-malign","GN",
 "GCP","GN","GN","non-malign","non-malign","non-malign")

for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  snRNA@meta.data$broad_CT <- str_replace_all(snRNA@meta.data$broad_CT, pattern, replacements[i])
}

#### Longitudinal analysis
celltypes_ind <- table(snRNA@meta.data$sample,snRNA@meta.data$broad_CT)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_broad_individual.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_broad_individual.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
  df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample","celltype","percentage")
library(ggplot2)

library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header= TRUE,sep="\t")

patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
# patterns<- gsub("-","_",patterns)
replacements = long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

### Patient ID
patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
# patterns<- gsub("-","_",patterns)
replacements = long_sample$Patient.ID
# names(replacement) <- patterns
df_melted$Patient.ID <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$Patient.ID <- str_replace_all(df_melted$Patient.ID, pattern, replacements[i])
}

# df_melted$celltype2 <- gsub("GCP_cycling","GCP",df_melted$celltype) %>% gsub("GNearly|GNlate","GN",.)
df_melted_paired <- df_melted[grep("7316-4529|7316-5881|7316-278|7316-333|7316-2118|7316-737|7316-3023|7316-931",df_melted$sample),]
df_melted_paired$Patient_Timepts <- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts, sep = "_")
# df_wide <- dcast(df_melted_paired, celltype ~ Patient_Timepts, value.var = "percentage", fun.aggregate=mean)
df_melted -> df_melted_paired

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test$p <- (stat.test$p)/2

bxp <- ggboxplot(
  df_melted_paired, x = "celltype", y = "percentage", 
  color = "timepts", palette = c("#00AFBB", "#E7B800")
) 

bxp2 <- bxp +   geom_dotplot(
  aes(fill = timepts, color = timepts), trim = FALSE,
  binaxis='y', stackdir='center', dotsize = 0.15,
  position = position_dodge(0.8)
)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +   theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white"))

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir,"Table/t_test_unpaired_primary_recurrent_samples_bxplot_two_sided_celltype_broad.pdf",sep = ""), width = 5.5, height = 4.5)
bxp4
dev.off()

