### Harmony Integration
library(Signac)
library(Seurat)
library(harmony)

combined <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/combined.RDS")
combined@meta.data$Sample <- gsub("_128.*.", "", rownames(combined@meta.data))

DefaultAssay(combined) <- "ATAC"

scATAC <- RunHarmony(
    combined,
    group.by.vars = "Sample",
    reduction.use = "lsi",
    assay.use = "ATAC",
    reduction.save = "ATACharmony",
    project.dim = FALSE
)

scATAC <- RunUMAP(
    scATAC,
    assay = "ATAC",
    reduction.key = "ATACharmonyUMAP_",
    reduction = "ATACharmony",
    reduction.name = "ATACharmonyumap",
    dims = 2:30
)

p1 <- DimPlot(scATAC, reduction = "ATACharmonyumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(scATAC, reduction = "ATACharmonyumap", label = FALSE, group.by = "Sample", split.by = "Sample")

p2 <- DimPlot(scATAC, reduction = "ATACharmonyumap", label = TRUE, group.by = "celltypes")
p3 <- DimPlot(scATAC, reduction = "ATACharmonyumap", label = TRUE, group.by = "Batch")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste0(savedir, "motif/chromvar/UMAP/batches_ATACharmony_sample_celltypes_batch.pdf"), width = 15, height = 5.5)
p1 + p2
dev.off()

# saveRDS(scATAC, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/ATAC_harmony.RDS")
