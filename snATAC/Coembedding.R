library(Seurat)
library(Signac)
library(future)

options(future.globals.maxSize = 50 * 1024^3)
plan("multicore", workers = 4)

RNA_integrated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed.RDS")
ATAC_integrated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_chromvar.RDS")

DefaultAssay(RNA_integrated) <- "RNA"
DefaultAssay(ATAC_integrated) <- "imputed_rna"

RNA_integrated <- FindVariableFeatures(RNA_integrated)
genes.use <- VariableFeatures(RNA_integrated)

### Since Seurat cannot embed such high number of cells
RNA_random_cells <- sample(rownames(RNA_integrated@meta.data), size = 25000, replace = F)
RNA_integrated_subset <- subset(RNA_integrated, cells = RNA_random_cells)

ATAC_random_cells <- sample(rownames(ATAC_integrated@meta.data), size = 25000, replace = F)
ATAC_integrated_subset <- subset(ATAC_integrated, cells = ATAC_random_cells)

coembed <- merge(x = RNA_integrated_subset, y = ATAC_integrated_subset)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

# DimPlot(coembed, group.by = c("orig.ident", "celltypes"))

saveRDS(coembed, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/scRNA_scATAC_coembedding.RDS")
