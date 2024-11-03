library(Signac)
library(Seurat)

Endo_peri_cellnames <- integrated@meta.data[grep("Endo-Pericytes", integrated@meta.data$celltypes), ] %>% row.names()
Endo_peri_obj <- subset(integrated, cells = Endo_peri_cellnames)

DefaultAssay(Endo_peri_obj) <- "ATAC"

Endo_peri_obj <- RunTFIDF(Endo_peri_obj)
Endo_peri_obj <- FindTopFeatures(Endo_peri_obj, min.cutoff = 20)
Endo_peri_obj <- RunSVD(Endo_peri_obj)
Endo_peri_obj <- RunUMAP(Endo_peri_obj, dims = 2:50, reduction = "lsi")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Endo_pericytes/"
dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/combined_sample.pdf"))
DimPlot(Endo_peri_obj, group.by = "dataset")
dev.off()

pdf(paste0(savedir, "UMAP/combined_sample_splitted.pdf"), width = 15, height = 12)
DimPlot(Endo_peri_obj, group.by = "dataset", split.by = "dataset", ncol = 4) + NoLegend()
dev.off()

DefaultAssay(Endo_peri_obj) <- "ATAC"

Endo_peri_obj <- RunHarmony(
    Endo_peri_obj,
    group.by.vars = "Sample",
    reduction.use = "lsi",
    assay.use = "ATAC",
    reduction.save = "ATACharmony",
    project.dim = FALSE
)

Endo_peri_obj <- RunUMAP(
    Endo_peri_obj,
    assay = "ATAC",
    reduction.key = "ATACharmonyUMAP_",
    reduction = "ATACharmony",
    reduction.name = "ATACharmonyumap",
    dims = 2:30
)

p1 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = FALSE, group.by = "Sample", split.by = "Sample")

# p2 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = TRUE, group.by = "celltypes")
# p3 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = TRUE, group.by = "Batch")

dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
pdf(paste0(savedir, "UMAP/batches_ATACharmony_sample_celltypes_batch.pdf"), width = 15, height = 5.5)
p1 + p2
dev.off()

rec_vs_pri <- FindMarkers(
    object = Endo_peri_obj,
    ident.1 = "Recurrent",
    ident.2 = "Primary",
    only.pos = FALSE,
    test.use = "LR",
    min.pct = 0.05,
    group.by = "tumor_type",
    latent.vars = "nCount_ATAC"
)

dir.create(paste0(savedir, "scdifferential/", sep = ""), showWarnings = FALSE)
write.table(rec_vs_pri, paste0(savedir, "scdifferential/Endo_pericytes_rec_vs_pri_ATAC_both.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

da_peaks <- rec_vs_pri
da_peaks_rec <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5 & da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
da_peaks_pri <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5 & da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2, ])

closest_genes_rec <- ClosestFeature(Endo_peri_obj, regions = da_peaks_rec)
closest_genes_prim <- ClosestFeature(Endo_peri_obj, regions = da_peaks_pri)

write.table(closest_genes_rec, paste0(savedir, "scdifferential/close_features/Endo_pericytes_close_gene_rec.txt"),
    quote = F, row.names = T, col.names = T, sep = "\t"
)
write.table(closest_genes_prim, paste0(savedir, "scdifferential/close_features/Endo_pericytes_close_gene_prim.txt"),
    quote = F, row.names = T, col.names = T, sep = "\t"
)

### Motif Enrichment
top.da.peak <- row.names(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2 & da_peaks$avg_log2FC < 0, ])
top.da.peak <- grep("^GL|^KI", top.da.peak, invert = TRUE, value = TRUE)
# test enrichment
enriched.motifs <- FindMotifs(
    object = Endo_peri_obj,
    features = top.da.peak
)

primary_motif_significant_fold_enrichment <- enriched.motifs[enriched.motifs$p.adjust < 0.01 & enriched.motifs$fold.enrichment > 3, ]
dir.create(paste(savedir, "motif", sep = ""))
write.table(enriched.motifs, paste0(savedir, "motif/Endothelial_primary_high_ATAC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(primary_motif_significant_fold_enrichment, paste0(savedir, "motif/Endothelial_primary_high_ATAC_significant.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

top.da.peak <- row.names(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2 & da_peaks$avg_log2FC > 0, ])
top.da.peak <- grep("^GL|^KI", top.da.peak, invert = TRUE, value = TRUE)
# test enrichment
enriched.motifs <- FindMotifs(
    object = Endo_peri_obj,
    features = top.da.peak
)

recurrent_motif_significant_fold_enrichment <- enriched.motifs[enriched.motifs$p.adjust < 0.01 & enriched.motifs$fold.enrichment > 3, ]

write.table(enriched.motifs, paste0(savedir, "motif/Endothelial_recurrent_high_ATAC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(recurrent_motif_significant_fold_enrichment, paste0(savedir, "motif/Endothelial_recurrent_high_ATAC_significant.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

p1 <- CoveragePlot(
    object = Endo_peri_obj,
    region = "KDR",
    features = "KDR",
    assay = "ATAC",
    expression.assay = "MAGIC_RNA",
    group.by = "Sample_pt_tumortype"
    # extend.upstream = 200000,
    # extend.downstream = 200000
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Endo_pericytes/"
dir.create(savedir)
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_MAGIC_RNA.pdf", sep = ""))
p1
dev.off()


dir.create(paste0(savedir, "saveRDS_obj"), showWarnings = FALSE, recursive = TRUE)
saveRDS(Endo_peri_obj, paste0(savedir, "saveRDS_obj/combined.RDS"))

### Trying with the predicted Endothelial cells
Endocells <- row.names(integrated@meta.data[grep("Endo", integrated@meta.data$predicted.id), ])
Endo_peri_obj <- subset(integrated, cells = Endocells)

DefaultAssay(Endo_peri_obj) <- "ATAC"

Endo_peri_obj <- RunTFIDF(Endo_peri_obj)
Endo_peri_obj <- FindTopFeatures(Endo_peri_obj, min.cutoff = 20)
Endo_peri_obj <- RunSVD(Endo_peri_obj)
Endo_peri_obj <- RunUMAP(Endo_peri_obj, dims = 2:50, reduction = "lsi")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Endo_pericytes/"
dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/combined_sample.pdf"))
DimPlot(Endo_peri_obj, group.by = "dataset")
dev.off()

pdf(paste0(savedir, "UMAP/combined_sample_splitted.pdf"), width = 15, height = 12)
DimPlot(Endo_peri_obj, group.by = "dataset", split.by = "dataset", ncol = 4) + NoLegend()
dev.off()

DefaultAssay(Endo_peri_obj) <- "ATAC"

Endo_peri_obj <- RunHarmony(
    Endo_peri_obj,
    group.by.vars = "Sample",
    reduction.use = "lsi",
    assay.use = "ATAC",
    reduction.save = "ATACharmony",
    project.dim = FALSE
)

Endo_peri_obj <- RunUMAP(
    Endo_peri_obj,
    assay = "ATAC",
    reduction.key = "ATACharmonyUMAP_",
    reduction = "ATACharmony",
    reduction.name = "ATACharmonyumap",
    dims = 2:30
)

p1 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = FALSE, group.by = "Sample", split.by = "Sample")

# p2 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = TRUE, group.by = "celltypes")
# p3 <- DimPlot(Endo_peri_obj, reduction = "ATACharmonyumap", label = TRUE, group.by = "Batch")

dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
pdf(paste0(savedir, "UMAP/batches_ATACharmony_sample_celltypes_batch.pdf"), width = 15, height = 5.5)
p1 + p2
dev.off()

rec_vs_pri <- FindMarkers(
    object = Endo_peri_obj,
    ident.1 = "Recurrent",
    ident.2 = "Primary",
    only.pos = FALSE,
    test.use = "LR",
    min.pct = 0.05,
    group.by = "tumor_type",
    latent.vars = "nCount_ATAC"
)

dir.create(paste0(savedir, "scdifferential/", sep = ""), showWarnings = FALSE)
write.table(rec_vs_pri, paste0(savedir, "scdifferential/Endo_pericytes_rec_vs_pri_ATAC_both.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

da_peaks <- rec_vs_pri
da_peaks_rec <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5 & da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
da_peaks_pri <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5 & da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2, ])

closest_genes_rec <- ClosestFeature(Endo_peri_obj, regions = da_peaks_rec)
closest_genes_prim <- ClosestFeature(Endo_peri_obj, regions = da_peaks_pri)

write.table(closest_genes_rec, paste0(savedir, "scdifferential/close_features/Endo_pericytes_close_gene_rec.txt"),
    quote = F, row.names = T, col.names = T, sep = "\t"
)
write.table(closest_genes_prim, paste0(savedir, "scdifferential/close_features/Endo_pericytes_close_gene_prim.txt"),
    quote = F, row.names = T, col.names = T, sep = "\t"
)

### Motif Enrichment
top.da.peak <- row.names(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2 & da_peaks$avg_log2FC < 0, ])
top.da.peak <- grep("^GL|^KI", top.da.peak, invert = TRUE, value = TRUE)
# test enrichment
enriched.motifs <- FindMotifs(
    object = Endo_peri_obj,
    features = top.da.peak
)

primary_motif_significant_fold_enrichment <- enriched.motifs[enriched.motifs$p.adjust < 0.01 & enriched.motifs$fold.enrichment > 3, ]
dir.create(paste(savedir, "motif", sep = ""))
write.table(enriched.motifs, paste0(savedir, "motif/Endothelial_primary_high_ATAC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(primary_motif_significant_fold_enrichment, paste0(savedir, "motif/Endothelial_primary_high_ATAC_significant.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

top.da.peak <- row.names(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2 & da_peaks$avg_log2FC > 0, ])
top.da.peak <- grep("^GL|^KI", top.da.peak, invert = TRUE, value = TRUE)
# test enrichment
enriched.motifs <- FindMotifs(
    object = Endo_peri_obj,
    features = top.da.peak
)

recurrent_motif_significant_fold_enrichment <- enriched.motifs[enriched.motifs$p.adjust < 0.01 & enriched.motifs$fold.enrichment > 1.5, ]

write.table(enriched.motifs, paste0(savedir, "motif/Endothelial_recurrent_high_ATAC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
write.table(recurrent_motif_significant_fold_enrichment, paste0(savedir, "motif/Endothelial_recurrent_high_ATAC_significant.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

p1 <- CoveragePlot(
    object = Endo_peri_obj,
    region = "KDR",
    features = "KDR",
    assay = "ATAC",
    expression.assay = "MAGIC_RNA",
    group.by = "Sample_pt_tumortype"
    # extend.upstream = 200000,
    # extend.downstream = 200000
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Endo_pericytes/"
dir.create(savedir)
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_MAGIC_RNA.pdf", sep = ""))
p1
dev.off()


dir.create(paste0(savedir, "saveRDS_obj"), showWarnings = FALSE, recursive = TRUE)
saveRDS(Endo_peri_obj, paste0(savedir, "saveRDS_obj/combined.RDS"))
