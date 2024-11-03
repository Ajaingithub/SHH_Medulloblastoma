### Extracting the Lymphoid Cells
### Extracting out the lymphoid cells to perform the differential
spec_immune_obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/saveRDS_obj/immune_clean.RDS")
lymphoid_cells <- rownames(spec_immune_obj@meta.data[grep("Tcells", spec_immune_obj@meta.data$subcelltypes2, ignore.case = TRUE), ])
lymphoid_obj <- subset(spec_immune_obj, cells = lymphoid_cells)

DefaultAssay(lymphoid_obj) <- "RNA"
lymphoid_obj <- NormalizeData(lymphoid_obj, normalization.method = "LogNormalize", scale.factor = 10000)
lymphoid_obj <- ScaleData(lymphoid_obj, features = rownames(lymphoid_obj))

lymphoid_rec_vs_pri <- FindMarkers(lymphoid_obj, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type", test.use = "MAST")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/lymphoid_cells/"

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(lymphoid_rec_vs_pri, paste(savedir, "table/lymphoid_rec_vs_prim_MAST.txt", sep = ""), quote = F, row.names = T, sep = "\t", col.names = T)

lymphoid_rec_vs_pri <- FindMarkers(lymphoid_obj, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/lymphoid_cells/"

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(lymphoid_rec_vs_pri, paste(savedir, "table/lymphoid_rec_vs_prim.txt", sep = ""), quote = F, row.names = T, sep = "\t", col.names = T)


### Natural Killer Cells high expression
genes <- c("KLRK1", "KLRD1", "NKG7", "RAB27A", "PRDM1", "TOX", "RASGRP1", "STAT5B", "ZNF683", "CD2", "IRF1")

### Extracting immune cells
library(Seurat)
library(dplyr)
DefaultAssay(lymphoid_obj) <- "RNA"

lymphoid_obj <- FindVariableFeatures(lymphoid_obj)
lymphoid_obj <- lymphoid_obj %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/lymphoid_cells/"
# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/EndoDimPlot_patient.pdf")
dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/ImmunePlot_patient.pdf", sep = ""), width = 6, height = 3)
DimPlot(lymphoid_obj, reduction = "umap", group.by = "patient_tumortype")
dev.off()

### It require integration
lymphoid_obj <- RunHarmony(
    lymphoid_obj,
    group.by.vars = "batch",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

lymphoid_obj <- RunUMAP(
    lymphoid_obj,
    assay = "RNA",
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

p1 <- DimPlot(lymphoid_obj, reduction = "harmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(lymphoid_obj, reduction = "harmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(lymphoid_obj, reduction = "harmonyumap", label = TRUE, group.by = "batch")

dir.create(paste0(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste0(savedir, "UMAP/batches_chromharmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

### Performing integration based on the samples
lymphoid_obj <- RunHarmony(
    lymphoid_obj,
    group.by.vars = "SampleID",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "sample_harmony"
)

lymphoid_obj <- RunUMAP(
    lymphoid_obj,
    assay = "RNA",
    reduction.key = "sampleharmonyUMAP_",
    reduction = "sample_harmony",
    reduction.name = "sampleharmonyumap",
    dims = 1:30
)

p1 <- DimPlot(lymphoid_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(lymphoid_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(lymphoid_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/sample_harmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

# Find neighbors based on the UMAP reduction
lymphoid_obj <- FindNeighbors(
    lymphoid_obj,
    reduction = "sample_harmony",
    dims = 1:30
)

# Perform clustering
lymphoid_obj <- FindClusters(
    lymphoid_obj,
    resolution = 1.1 # You can adjust the resolution parameter based on your needs
)

# Optional: Visualize the clusters using UMAP
pdf(paste0(savedir, "UMAP/sampleharmony_dinmplot.pdf"), width = 3, height = 3)
DimPlot(lymphoid_obj, reduction = "sampleharmonyumap", label = TRUE)
dev.off()

pdf(paste0(savedir, "UMAP/sampleharmony_dimplot_splitted.pdf"), width = 15, height = 15)
DimPlot(lymphoid_obj, reduction = "sampleharmonyumap", group.by = "celltype_pbmc_predicted", split.by = "orig.ident", ncol = 4)
dev.off()

pdf(paste0(savedir, "UMAP/celltypes_dimplot_splitted.pdf"), width = 15, height = 15)
DimPlot(lymphoid_obj, reduction = "sampleharmonyumap", group.by = "celltype_pbmc_predicted", split.by = "celltype_pbmc_predicted", ncol = 4)
dev.off()

lymphoid_obj_markers <- FindAllMarkers(lymphoid_obj)
dir.create(paste(savedir, "Table", sep = ""), showWarnings = FALSE)
write.table(lymphoid_obj_markers, paste(savedir, "table/lymphoid_obj_cell_markers.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

saveRDS(lymphoid_obj@assays$RNA@counts, paste(savedir, "table/lymhoid_counts.RDS", sep = ""))


genes <- c("CD70", "GNLY", "CCL4", "KLRK1", "KLRD1", "CD8A", "NKG7") ### NK cells gene high
genes <- c("ATG16L1", "ATG7", "ATG4D", "UBE2D1", "UBE2G1", "UBE2L3") ## Autophagy and Ubiquitination
genes <- c("SOX6", "SOX5")

lymphoid_paired_cells <- rownames(lymphoid_obj@meta.data[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", lymphoid_obj@meta.data$sample2), ])
lymphoid_obj@meta.data$patient_tumortype <- paste(lymphoid_obj@meta.data$patient, lymphoid_obj@meta.data$tumor_type, sep = "_")
lymphoid_obj_paired <- subset(lymphoid_obj, cells = lymphoid_paired_cells)
DefaultAssay(lymphoid_obj_paired) <- "RNA"
lymphoid_obj_paired <- NormalizeData(lymphoid_obj_paired)

DefaultAssay(lymphoid_obj) <- "MAGIC_RNA"
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(lymphoid_obj_paired, genes[i],
        pt.size = 0,
        group.by = "patient_tumortype",
        assay = "MAGIC_RNA",
        cols = c(
            "lightblue3", "lightcoral",
            "lightblue3", "lightcoral",
            "lightblue3", "lightblue3", "lightcoral",
            "lightblue3"
        )
    ) + geom_boxplot()
}

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/lymphoid_genes_paired_high_primary_SOX.pdf", sep = ""))
plot_list
dev.off()

DefaultAssay(lymphoid_obj) <- "MAGIC_RNA"
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(lymphoid_obj, genes[i],
        pt.size = 0,
        group.by = "tumor_type",
        assay = "MAGIC_RNA",
        # cols = c("lightblue3", "lightcoral", "lightblue3", "lightcoral", "lightblue3", "lightcoral", "lightblue3", "lightcoral")
    ) + geom_boxplot()
}

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/lymphoid_genes_high_recurrent_tumor_type.pdf", sep = ""))
plot_list
dev.off()

genes <- c(
    "CD70", "GNLY", "CCL4", "KLRK1", "KLRD1", "CD8A", "NKG7",
    "HAVCR2", "LAG3", "PDCD1", "CTLA4", "GZMB", "PRF1", "TIGIT"
)


dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/lymphoid_genes_celltypes_MAGIC_RNA.pdf", sep = ""), width = 5, height = 6)
DotPlot(lymphoid_obj, features = genes, group.by = "seurat_clusters", assay = "MAGIC_RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()
