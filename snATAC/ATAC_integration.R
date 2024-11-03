library("Seurat")
library("dplyr")
library("Signac")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
library("ggplot2")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
integrated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_with_RNA_with_imputation_and_motifs.RDS")
snRNA_obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed.RDS")
DefaultAssay(snRNA_obj) <- "RNA"
DefaultAssay(integrated) <- "RNA"
snRNA_obj <- FindVariableFeatures(snRNA_obj, nfeatures = 3000)

DefaultAssay(integrated) <- "RNA"
integrated <- FindVariableFeatures(integrated, nfeatures = 4000)

transfer.anchors <- FindTransferAnchors(
    reference = snRNA_obj,
    query = integrated,
    reduction = "cca"
)

saveRDS(transfer.anchors, paste0(savedir, "anchors_snRNA_integrated_cca.RDS"))

predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = snRNA_obj$sub_celltypes,
    weight.reduction = integrated[["integrated_lsi"]],
    dims = 2:30
)

predicted.labels2 <- TransferData(
    anchorset = transfer.anchors,
    refdata = snRNA_obj$celltypes,
    weight.reduction = integrated[["integrated_lsi"]],
    dims = 2:30
)

integrated <- AddMetaData(object = integrated, metadata = predicted.labels)

integrated@meta.data$predicted_celltypes <- predicted.labels2$predicted.id

plot1 <- DimPlot(
    object = integrated,
    group.by = "predicted_celltypes",
    label = TRUE,
    repel = TRUE
)

plot2 <- DimPlot(
    object = integrated,
    group.by = "predicted.id",
    label = TRUE,
    repel = TRUE
)

pdf(paste0(savedir, "/UMAP/integrated_projection_CCA.pdf"), width = 12, height = 5.5)
print(plot1 + plot2)
dev.off()

saveRDS(integrated, paste0(savedir, "/saveRDS_obj/integrated_QC_snRNA_ref_mapped.RDS"))
