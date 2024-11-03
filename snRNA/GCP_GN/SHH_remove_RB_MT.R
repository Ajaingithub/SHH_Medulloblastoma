#### 3rd Time
### Removing the MT and RB gene since it is snRNA and it could be an artifact
library(Seurat)
library(dplyr)
library(glmGamPoi)

obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/obj_integrated_RNA_NN_cluster_0.6.RDS")

DefaultAssay(obj) <- "RNA"
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_and_MT_RB/"
dir.create(savedir, showWarnings=FALSE)

# Assuming 'seurat_obj' is your Seurat object
all_genes <- rownames(obj)
mito_genes <- grep("^MT-", all_genes, value = TRUE)
ribo_genes <- grep("^RPS|^RPL", all_genes, value = TRUE)
genes_to_remove <- c(mito_genes, ribo_genes)

obj_subset <- subset(obj, features = setdiff(all_genes, genes_to_remove))
removed_genes <- intersect(rownames(obj_subset), genes_to_remove) ## Sanity check
print(removed_genes)

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
res = c(0.2,0.4,0.6,0.8,1)
# res = 0.8
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD4_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD4_integrated_RNA, dims = 30, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA","S.Score","G2M.Score"), objname = "obj_subset",
                                                       process = process, col_sel = c("orig.ident","tumor_type","sex","age"))
}

dir.create(paste(savedir,"saveRDS_obj",sep=""),showWarning=FALSE)
saveRDS(CD4_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/obj_integrated_RNA_NN_cluster_1.RDS",sep = ""))

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
GNearly <- c("NHLH2","NHLH1","STMN2","TUBB3","NRG1","GRIN2B","KIF21B","DISP3")
GNlate <- c("NEUROD1","NEUROD2","OXR1","NRXN2","RBFOX3","BRINP1","GRIN2C")
Neuronal <- c("TRPM3","OTX2","ROBO3","UNC5D","ITGA3","ADRA1A","OTX2-AS1","PRDM6","GRIN2D")
Stem_cells <- c("NES","SOX2","MKI67","TOP2A")
HSP <- c("HSPH1","HSPA1B","HSPA1A","DNAJB1","HSPB1")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[,1]
genes <- c(GCP,GNearly, GNlate, Neuronal, Stem_cells, HSP, nonmalignant_cells)

DefaultAssay(obj_impute) <- "MAGIC_RNA"
p <- DotPlot(obj_impute, feature = genes, group.by = "seurat_clusters2", assay = "MAGIC_RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +  coord_flip()  +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/obj_marker_2_genes_0.8.pdf", height = 24, width =11)
p
dev.off()
