#### Extracting out the immune cells
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"

immune_cells <- rownames(snRNA@meta.data[grep("Immune", snRNA@meta.data$celltypes), ])

p <- DimPlot(snRNA,
    cells.highlight = immune_cells,
    reduction = "umap",
    label = FALSE, cols.highlight = "deeppink2",
    sizes.highlight = 0.5,
    cols = "gray92"
) +
    ggtitle("Immune Cells") +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.2))

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/snRNA_immune_cells.pdf", sep = ""))
p
dev.off()

genes <- c("GLI2", "CUX2", "GRIN2A", "SMO")
source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    p <- featureplot_front(snRNA, genes[i],
        reduction = "umap", x = "umap_1", y = "umap_2", size = 0.
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

dir.create(paste(savedir, "featureplot/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/Immune_GCP_genes.pdf", sep = ""), width = 5.5, height = 5)
print(plot_list)
dev.off()


### Extracting immune cells
library(Seurat)
library(Signac)

immune_cells <- rownames(snRNA@meta.data[grep("Immune", snRNA@meta.data$celltypes), ])
immune <- subset(snRNA, cells = immune_cells)

DefaultAssay(immune) <- "RNA"

immune <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)
immune <- ScaleData(immune, features = rownames(immune))

immune <- FindVariableFeatures(immune)
immune <- immune %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/EndoDimPlot_patient.pdf")
dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/ImmunePlot_patient.pdf", sep = ""))
DimPlot(immune, reduction = "umap", group.by = "patient_tumortype")
dev.off()

### It require integration
immune <- RunHarmony(
    immune,
    group.by.vars = "batch",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

immune <- RunUMAP(
    immune,
    assay = "RNA",
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

p1 <- DimPlot(immune, reduction = "harmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(immune, reduction = "harmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(immune, reduction = "harmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/batches_chromharmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

### Performing integration based on the samples
immune <- RunHarmony(
    immune,
    group.by.vars = "SampleID",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "sample_harmony"
)

immune <- RunUMAP(
    immune,
    assay = "RNA",
    reduction.key = "sampleharmonyUMAP_",
    reduction = "sample_harmony",
    reduction.name = "sampleharmonyumap",
    dims = 1:30
)

p1 <- DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/sample_harmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

# Find neighbors based on the UMAP reduction
immune <- FindNeighbors(
    immune,
    reduction = "sample_harmony",
    dims = 1:30
)

# Perform clustering
immune <- FindClusters(
    immune,
    resolution = 0.6 # You can adjust the resolution parameter based on your needs
)

# Optional: Visualize the clusters using UMAP
pdf(paste0(savedir, "UMAP/sampleharmony_dinmplot.pdf"), width = 5.5, height = 5.5)
DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE)
dev.off()

pdf(paste0(savedir, "UMAP/sampleharmony_dimplot_splitted.pdf"), width = 15, height = 15)
DimPlot(immune, reduction = "sampleharmonyumap", group.by = "celltype_pbmc_predicted", split.by = "orig.ident", ncol = 4)
dev.off()

pdf(paste0(savedir, "UMAP/celltypes_dimplot_splitted.pdf"), width = 15, height = 15)
DimPlot(immune, reduction = "sampleharmonyumap", group.by = "celltype_pbmc_predicted", split.by = "celltype_pbmc_predicted", ncol = 4)
dev.off()


immune_markers <- FindAllMarkers(immune)
dir.create(paste(savedir, "Table", sep = ""), showWarnings = FALSE)
write.table(immune_markers, paste(savedir, "Table/immune_cell_markers.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


### Reference mapping to the multimodal data
immune <- SCTransform(immune, verbose = FALSE)
anchors <- FindTransferAnchors(
    reference = reference,
    query = immune,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
)

# We then transfer cell type labels and protein data from the reference to the query. Additionally, we project the query data onto the UMAP structure of the reference.
immune <- MapQuery(
    anchorset = anchors,
    query = immune,
    reference = reference,
    refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
)

#### Generate the prediction from the Azimuth
az_pred <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/Table/azimuth_pred.tsv", sep = "\t", header = TRUE)
stopifnot(all(rownames(immune@meta.data) == az_pred$cell))
immune@meta.data$celltype_pbmc_predicted <- az_pred$predicted.celltype.l2

pdf(paste0(savedir, "UMAP/sampleharmony_dinmplot_celltypes.pdf"), width = 5.5, height = 5.5)
DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "celltype_pbmc_predicted")
dev.off()

#### Longitudinal analysis
celltypes_ind <- table(immune@meta.data$SampleID, immune@meta.data$seurat_clusters)

write.table(celltypes_ind,
    paste(savedir, "Table/celltype_individual.txt", sep = ""),
    row.names = T,
    col.names = T,
    sep = "\t",
    quote = F
)

celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/Table/celltype_individual.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

### Sanity Check
rowSums(df)

### Making some change in the naming of the samples
rownames(df) <- gsub("-", "_", rownames(df))

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample", "celltype", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c(
    "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
    "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
    "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994", "SF8368", "SF8539", "SF9232"
) -> sample_levels

df_melted$sample <- gsub("-", "_", df_melted$sample) %>% gsub("SF7994R", "SF7994", .)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = sample)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"
pdf(paste(savedir, "Table/celltype_individual_barplot.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.delim("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
replacements <- long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

# df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]

df_melted -> df_melted_paired

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    # t_test(percentage ~ timepts, paired = TRUE, alternative = "greater") %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2

bxp <- ggboxplot(
    df_melted_paired,
    x = "celltype", y = "percentage",
    color = "timepts", palette = c("#00AFBB", "#E7B800")
)

bxp2 <- bxp + geom_dotplot(
    aes(fill = timepts, color = timepts),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.grid.major = element_line(colour = "white")
)


# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"
pdf(paste(savedir, "Table/t_test_all_primary_recurrent_bxplot_one_sided_greater.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()

### barplot
gbr <- ggbarplot(df_melted,
    x = "celltype", y = "percentage",
    add = c("mean_se", "jitter"),
    color = "timepts", palette = c("#00AFBB", "#E7B800"),
    position = position_dodge(0.8)
)

stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
gbr2 <- gbr + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")
    )

pdf(paste(savedir, "Table/t_test_all_primary_recurrent_barplot_t_test.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]
gbr <- ggbarplot(df_melted_paired,
    x = "celltype", y = "percentage",
    add = c("mean_se", "jitter"),
    color = "timepts", palette = c("#00AFBB", "#E7B800"),
    position = position_dodge(0.8)
)

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    # t_test(percentage ~ timepts, paired = TRUE, alternative = "greater") %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2


stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
gbr2 <- gbr + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")
    )
d
pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_barplot_t_test.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

### Identified the GCP genes in the immune cells
genes <- c("GLI2", "CUX2", "GRIN2A", "SMO")
genes <- c("PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L", "ATOH1", "GLI2", "GLI1", "CD163", "CD4", "CD8A", "CD86", "CD74")
source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    p <- featureplot_front(snRNA, genes[i],
        reduction = "umap", x = "umap_1", y = "umap_2", size = 0.
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

dir.create(paste(savedir, "featureplot/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/Immune_GCP_and_immune_genes.pdf", sep = ""), width = 5.5, height = 5)
print(plot_list)
dev.off()

DefaultAssay(immune) <- "MAGIC_RNA"
# genes <- c("GLI2", "CUX2", "GRIN2A", "SMO", "CD4", "CD14", "CD163", "CD69")
genes <- c(
    "PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L",
    "ATOH1", "GLI2", "GLI1"
)

genes <- c(
    "CD163", "CD4", "CD8A", "CD86", "CD74", "CD14", "CLEC7A",
    "ITGAX", "LY86", "IL2RA", "PTPRC", "NKG7", "CD69"
)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    p <- featureplot_front(immune, genes[i],
        reduction = "sampleharmonyumap", x = "sampleharmonyUMAP_1", y = "sampleharmonyUMAP_2", size = 0.4
    ) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

genes <- c(
    "PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L",
    "ATOH1", "GLI2", "GLI1", "CD163", "CD4", "CD8A", "CD86", "CD74", "CD14", "CLEC7A",
    "ITGAX", "LY86", "IL2RA", "PTPRC", "NKG7", "CD69"
)

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/GCP_and_immune_genes_celltypes.pdf", sep = ""), width = 9, height = 8)
DotPlot(immune, features = genes, group.by = "subcelltypes") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

### Since I am getting the GCP cells also in the Immune cells. Need to check why they are coming.
p <- VlnPlot(immune, features = c("nCount_RNA", "nFeature_RNA", "doublet_score"), group.by = "subcelltypes")
dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/VlnPlot/")
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/immune_featureplot_subcelltypes.pdf")
p
dev.off()

immune_GCP_cells <- rownames(immune@meta.data[grep("GCP", immune@meta.data$subcelltypes), ])

p <- DimPlot(snRNA,
    cells.highlight = immune_GCP_cells,
    reduction = "umap",
    label = FALSE, cols.highlight = "deeppink2",
    sizes.highlight = 0.5,
    cols = "gray92"
) +
    ggtitle("Immune Cells") +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.2))

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/snRNA_immune_GCP_cells.pdf", sep = ""))
p
dev.off()

snRNA@meta.data$celltypes2 <- snRNA@meta.data$celltypes
snRNA@meta.data[match(immune_GCP_cells, rownames(snRNA@meta.data)), "celltypes2"] <- "Immune_GCP"

DefaultAssay(snRNA) <- "RNA"
immune_GCP_vs_GCP <- FindMarkers(snRNA, ident.1 = "Immune_GCP", ident.2 = "GCP", group.by = "celltypes2")
write.table(immune_GCP_vs_GCP, paste(savedir, "Table/immune_GCP_vs_GCP.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Extracting randomly 1000 GCP cells from the main object
GCP_random_cells <- sample(rownames(snRNA@meta.data[grep("^GCP$", snRNA@meta.data$celltypes), ]), 1000, replace = FALSE)
immune_cells <- rownames(snRNA@meta.data[grep("Immune", snRNA@meta.data$celltypes), ])
GCP_immune_cells <- c(GCP_random_cells, immune_cells)

#### Extracting out the immune cells
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/immune_random_GCP/"
dir.create(savedir, showWarnings = FALSE)

p <- DimPlot(snRNA,
    cells.highlight = GCP_immune_cells,
    reduction = "umap",
    label = FALSE, cols.highlight = "deeppink2",
    sizes.highlight = 0.5,
    cols = "gray92"
) +
    ggtitle("Immune Cells") +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.2))

dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/snRNA_GCP_immune_cells.pdf", sep = ""))
p
dev.off()

immune <- subset(snRNA, cells = GCP_immune_cells)

DefaultAssay(immune) <- "RNA"

immune <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)
immune <- ScaleData(immune, features = rownames(immune))

immune <- FindVariableFeatures(immune)
immune <- immune %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/EndoDimPlot_patient.pdf")
dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/ImmunePlot_patient.pdf", sep = ""))
DimPlot(immune, reduction = "umap", group.by = "patient_tumortype")
dev.off()

### It require integration
immune <- RunHarmony(
    immune,
    group.by.vars = "batch",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

immune <- RunUMAP(
    immune,
    assay = "RNA",
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

p1 <- DimPlot(immune, reduction = "harmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(immune, reduction = "harmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(immune, reduction = "harmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/batches_chromharmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

### Performing integration based on the samples
immune <- RunHarmony(
    immune,
    group.by.vars = "SampleID",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "sample_harmony"
)

immune <- RunUMAP(
    immune,
    assay = "RNA",
    reduction.key = "sampleharmonyUMAP_",
    reduction = "sample_harmony",
    reduction.name = "sampleharmonyumap",
    dims = 1:30
)

p1 <- DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/sample_harmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

# Find neighbors based on the UMAP reduction
immune <- FindNeighbors(
    immune,
    reduction = "sample_harmony",
    dims = 1:30
)

# Perform clustering
immune <- FindClusters(
    immune,
    resolution = 0.6 # You can adjust the resolution parameter based on your needs
)

# Optional: Visualize the clusters using UMAP
pdf(paste0(savedir, "UMAP/sampleharmony_dinmplot.pdf"), width = 5.5, height = 5.5)
DimPlot(immune, reduction = "sampleharmonyumap", label = TRUE)
dev.off()

pdf(paste0(savedir, "UMAP/celltypes_dimplot.pdf"))
DimPlot(immune, reduction = "sampleharmonyumap", group.by = "celltypes2")
dev.off()

pdf(paste0(savedir, "UMAP/celltypes_dimplot_splitted.pdf"), width = 15, height = 15)
DimPlot(immune, reduction = "sampleharmonyumap", group.by = "celltypes2", split.by = "celltypes2", ncol = 4)
dev.off()

immune_markers <- FindAllMarkers(immune)
dir.create(paste(savedir, "Table", sep = ""), showWarnings = FALSE)
write.table(immune_markers, paste(savedir, "Table/immune_cell_markers.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

p <- VlnPlot(immune, features = c("nCount_RNA", "nFeature_RNA", "doublet_score"), group.by = "celltypes2")
dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/VlnPlot/")
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/VlnPlot/immune_featureplot_subcelltypes.pdf")
p
dev.off()

genes <- c(
    "PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L",
    "ATOH1", "GLI2", "GLI1", "CD163", "CD4", "CD8A", "CD86", "CD74", "CD14", "CLEC7A",
    "ITGAX", "LY86", "IL2RA", "PTPRC", "NKG7", "CD69"
)

genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/Table/genes_markers.txt")
genes <- unique(genes[, 1])

spec_immune_obj@meta.data$subcelltypes <- factor(spec_immune_obj@meta.data$subcelltypes,
    levels = c(
        "M1 macrophages", "M2 macrophages", "Microglia", "Tumor Associated Macrophages",
        "Neutrophils", "Macrophages", "Lymphoid", "Proliferating Cells", "Unknown", "Neuronal Cells"
    )
)

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/immune_genes_celltypes_MAGIC_RNA.pdf", sep = ""), width = 9, height = 16)
DotPlot(spec_immune_obj, features = genes, group.by = "subcelltypes", assay = "MAGIC_RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/immune_genes_celltypes_MAGIC_RNA_cluster_2.pdf", sep = ""), width = 9, height = 14)
DotPlot(spec_immune_obj, features = genes, group.by = "seurat_clusters", assay = "MAGIC_RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

## Removing GCP cluster from the immune cells
immune_spec <- rownames(immune@meta.data[grep("^1$|^3$|^6$", immune@meta.data$seurat_clusters, invert = TRUE), ])
spec_immune_obj <- subset(immune, cells = immune_spec)

DefaultAssay(spec_immune_obj) <- "RNA"
spec_immune_obj <- NormalizeData(spec_immune_obj, normalization.method = "LogNormalize", scale.factor = 10000)
spec_immune_obj <- ScaleData(spec_immune_obj, features = rownames(spec_immune_obj))

spec_immune_obj <- FindVariableFeatures(spec_immune_obj)
spec_immune_obj <- spec_immune_obj %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
# pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/EndoDimPlot_patient.pdf")
dir.create(paste(savedir, "UMAP", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "UMAP/ImmunePlot_patient.pdf", sep = ""))
DimPlot(spec_immune_obj, reduction = "umap", group.by = "patient_tumortype")
dev.off()

##
spec_immune_obj <- RunHarmony(
    spec_immune_obj,
    group.by.vars = "batch",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

### It require integration
spec_immune_obj <- RunHarmony(
    spec_immune_obj,
    group.by.vars = "batch",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "harmony"
)

spec_immune_obj <- RunUMAP(
    spec_immune_obj,
    assay = "RNA",
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

p1 <- DimPlot(spec_immune_obj, reduction = "harmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(spec_immune_obj, reduction = "harmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(spec_immune_obj, reduction = "harmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/batches_chromharmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

### Performing integration based on the samples
spec_immune_obj <- RunHarmony(
    spec_immune_obj,
    group.by.vars = "SampleID",
    reduction.use = "pca",
    assay.use = "RNA",
    reduction.save = "sample_harmony"
)

spec_immune_obj <- RunUMAP(
    spec_immune_obj,
    assay = "RNA",
    reduction.key = "sampleharmonyUMAP_",
    reduction = "sample_harmony",
    reduction.name = "sampleharmonyumap",
    dims = 1:30
)

p1 <- DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "SampleID")
p2 <- DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "patient")
p3 <- DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "batch")

pdf(paste0(savedir, "UMAP/sample_harmony_sample_age_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

# Find neighbors based on the UMAP reduction
spec_immune_obj <- FindNeighbors(
    spec_immune_obj,
    reduction = "sample_harmony",
    dims = 1:30
)

# Perform clustering
# res <- c(0.4, 0.5, 0.6, 0.7, 0.8)
# res <- c(0.9, 1, 1.1, 1.2, 1.3)
res <- 1.3
for (i in 1:length(res)) {
    spec_immune_obj <- FindClusters(
        spec_immune_obj,
        resolution = res[i]
    )
    pdf(paste0(savedir, "UMAP/sampleharmony_dinmplot", res[i], ".pdf"), width = 5.5, height = 5.5)
    print(DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", label = TRUE))
    dev.off()
}

pdf(paste0(savedir, "UMAP/celltypes_dimplot_splitted.pdf"), width = 15, height = 15)
DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", group.by = "celltype_pbmc_predicted", split.by = "celltype_pbmc_predicted", ncol = 4)
dev.off()

pdf(paste0(savedir, "UMAP/celltypes_dimplot.pdf"), width = 5, height = 5)
DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", group.by = "celltype_pbmc_predicted")
dev.off()

immune_markers <- FindAllMarkers(spec_immune_obj)
dir.create(paste(savedir, "Table", sep = ""), showWarnings = FALSE)
write.table(immune_markers, paste(savedir, "Table/immune_cell_markers_0.6.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

### Imputing the data
DefaultAssay(obj) <- "RNA"
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check
obj <- magic(obj, npca = 20) ## imputing the RNA data as for RNA PCs are 20
DefaultAssay(obj) <- "MAGIC_RNA"
obj <- ScaleData(obj, features = rownames(obj)) # Sacling data on all the features as we need it during the feature plot scaled data for gene that are not variable

#### Performing Add Module Score to identify the M0, M1, and M2 based on the Aaron Nature Cancer paper
library(ArchR)
library(UCell)
library(ggplot2)

DefaultAssay(obj) <- "MAGIC_RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/xCell_genes/", pattern = ".txt", full.names = TRUE)
filename <- paste(gsub(".txt|_cell|_genes", "", basename(files)), "_xCell", sep = "")
# gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i], sep = "\t")[, 1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset, rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

names(markers) <- paste(names(markers), "data", sep = "_")
spec_immune_obj <- AddModuleScore_UCell(spec_immune_obj, features = markers, ncores = 4, slot = "data", assay = "RNA")
featurename <- grep("_xCell_data_UCell", colnames(spec_immune_obj@meta.data), value = TRUE)

p1 <- VlnPlot(spec_immune_obj, features = featurename, group.by = "seurat_clusters")
# p3 <- VlnPlot(obj, features = featurename[1], group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()
# p4 <- VlnPlot(obj, features = featurename[2], group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(spec_immune_obj, features = featurename[i], reduction = "sampleharmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/myeloid_Ucell_xcell_data_RNA.pdf", sep = ""), width = 15, height = 11)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/myeloid_Ucell_xcell_data_RNA.pdf", sep = ""))
print(plot_list)
dev.off()

### Add Module score Seurat
spec_immune_obj <- AddModuleScore(spec_immune_obj, features = markers, slot = "scale.data")
colnames(spec_immune_obj@meta.data)[69:79] <- paste(featurename, "_score", sep = "")

p1 <- VlnPlot(spec_immune_obj, features = paste(featurename, "_score", sep = ""), group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(spec_immune_obj,
        features = paste(featurename, "_score", sep = "")[i],
        reduction = "sampleharmonyumap"
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/myeloid_Ucell_Abhinav_Addmod.pdf", sep = ""), width = 15, height = 11)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/myeloid_Ucell_Abhinav_Addmod.pdf", sep = ""))
print(plot_list)
dev.off()

spec_immune_obj <- AddModuleScore(spec_immune_obj, features = markers)
colnames(spec_immune_obj@meta.data)[80:90] <- paste(featurename, "_log_score", sep = "")

p1 <- VlnPlot(spec_immune_obj, features = paste(featurename, "_log_score", sep = ""), group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(spec_immune_obj,
        features = paste(featurename, "_log_score", sep = "")[i],
        reduction = "sampleharmonyumap"
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/myeloid_Ucell_Abhinav_log_Addmod.pdf", sep = ""), width = 15, height = 11)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/myeloid_Ucell_Abhinav_log_Addmod.pdf", sep = ""))
print(plot_list)
dev.off()

#### Abhinav including Nature Cancer Aaron markers
library(ArchR)
library(UCell)
library(ggplot2)

obj <- spec_immune_obj
obj@meta.data <- obj@meta.data[, 1:57]

DefaultAssay(obj) <- "MAGIC_RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/", pattern = ".txt", full.names = TRUE)
filename <- paste(gsub(".txt|_cell|_genes", "", basename(files)), "_Abhinav", sep = "")
# gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i], sep = "\t")[, 1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset, rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

# names(markers) <- paste(names(markers),"data",sep="_")
obj <- AddModuleScore_UCell(obj, features = markers, ncores = 4, slot = "scale.data", assay = "MAGIC_RNA")
featurename <- c("M1_Abhinav_UCell", "M2_Abhinav_UCell")

p1 <- VlnPlot(obj, features = featurename, group.by = "seurat_clusters")
# p3 <- VlnPlot(obj, features = featurename[1], group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()
# p4 <- VlnPlot(obj, features = featurename[2], group.by = "seurat_clusters", pt.size = 0) + geom_boxplot()

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(obj, features = featurename[i], reduction = "sampleharmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/Abhinav_M1_M2_scale_data_MAGIC_RNA.pdf", sep = ""), width = 8, height = 5)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/Abhinav_M1_M2_scale_data_MAGIC_RNA.pdf", sep = ""))
print(plot_list)
dev.off()

### Add Module score Seurat
names(markers) <- paste(names(markers), "scale", sep = "_")
obj <- spec_immune_obj
obj@meta.data <- obj@meta.data[, 1:57]
obj <- AddModuleScore(obj, features = markers, slot = "data")
colnames(obj@meta.data)[58:59] <- paste(featurename, "_score", sep = "")

p1 <- VlnPlot(obj, features = paste(featurename, "_score", sep = ""), group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(obj,
        features = paste(featurename, "_score", sep = "")[i],
        reduction = "sampleharmonyumap"
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/Abhinav_M1_M2_Ucell_Addmod_2.pdf", sep = ""), width = 15, height = 11)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/Abhinav_M1_M2_Ucell_Addmod_2.pdf", sep = ""))
print(plot_list)
dev.off()

### Now for the data score for AddModuleScore making the M1 and M2
obj@meta.data$macrophage_subtype <- "M0"
obj@meta.data[obj@meta.data$M1_Abhinav_UCell_score > obj@meta.data$M2_Abhinav_UCell_score, "macrophage_subtype"] <- "M1"
obj@meta.data[obj@meta.data$M1_Abhinav_UCell_score < obj@meta.data$M2_Abhinav_UCell_score, "macrophage_subtype"] <- "M2"
obj@meta.data[obj@meta.data$M2_Abhinav_UCell_score < 0 & obj@meta.data$M1_Abhinav_UCell_score < 0, "macrophage_subtype"] <- "M0"
obj@meta.data[grep("^1$", obj@meta.data$seurat_clusters), "macrophage_subtype"] <- "T_cells"
obj@meta.data[grep("^9$", obj@meta.data$seurat_clusters), "macrophage_subtype"] <- "Neuronal"

pdf(paste(savedir, "UMAP/macrophage_subtype.pdf", sep = ""))
DimPlot(obj, group.by = "macrophage_subtype", reduction = "sampleharmonyumap") + scale_color_manual(values = c(
    "#279e68", "#d62728", "orange", "#1f77b4", "#aa40fc",
    "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white"
))
dev.off()

pdf(paste(savedir, "UMAP/seurat_cluster.pdf", sep = ""))
DimPlot(obj, reduction = "sampleharmonyumap", label = TRUE)
dev.off()

## Finding the cluster distribution
seurat_macrophage <- table(obj@meta.data$seurat_clusters, obj@meta.data$macrophage_subtype)

write.table(seurat_macrophage,
    "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/Table/seurat_macrophages.txt",
    quote = F, col.names = T, row.names = T, sep = "\t"
)

seurat_macrophage <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/Table/seurat_macrophages.txt", header = TRUE, row.names = 1)

n_cells <- t(seurat_macrophage)
n_cells_sum <- as.vector(rowSums(n_cells))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

### Sanity Check
rowSums(df)

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("celltypes", "seurat_cluster", "percentage")

p <- ggplot(df_melted, aes(fill = celltypes, y = percentage, x = seurat_cluster)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    theme(axis.text = element_text(angle = 0, vjust = 0, hjust = 0, size = 15))

pdf(paste(savedir, "Table/macrophage_seurat_clusters_2.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

df_melted_noT_noNeuronal <- df_melted[grep("^1$|^9$", df_melted$seurat_cluster, invert = TRUE), ]
p <- ggplot(df_melted_noT_noNeuronal, aes(fill = celltypes, y = percentage, x = seurat_cluster)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    theme(axis.text = element_text(angle = 0, vjust = 0, hjust = 0, size = 15))

pdf(paste(savedir, "Table/macrophage_seurat_clusters_noT_Neuronal.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

clus_CT <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/Table/cluster_celltypes.txt", header = TRUE)
patterns <- clus_CT$Cluster
replacements <- clus_CT$Celltypes
# names(replacement) <- patterns
obj@meta.data$subcelltypes2 <- obj@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    obj@meta.data$subcelltypes2 <- str_replace_all(obj@meta.data$subcelltypes2, pattern, replacements[i])
}

pdf(paste(savedir, "UMAP/subcelltypes2.pdf", sep = ""))
DimPlot(obj, group.by = "subcelltypes2", reduction = "sampleharmonyumap", label = TRUE)
dev.off()

# obj@meta.data$mac_subcelltypes <- paste(obj@meta.data$macrophage_subtype, obj@meta.data$subcelltypes2, sep = "_")

# pdf(paste(savedir, "UMAP/mac_subcelltypes2.pdf", sep = ""))
# DimPlot(obj, group.by = "mac_subcelltypes", reduction = "sampleharmonyumap")
# dev.off()

obj@meta.data$SampleID <- gsub("-", "_", obj@meta.data$SampleID)
celltypes_ind <- table(obj@meta.data$SampleID, obj@meta.data$subcelltypes2)

n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

### Sanity Check
rowSums(df)

### Making some change in the naming of the samples
rownames(df) <- gsub("-", "_", rownames(df))

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample", "celltype", "percentage")
library(ggplot2)
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.delim("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

patterns <- long_sample$Sample.ID
# patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
patterns <- gsub("-", "_", patterns)
replacements <- long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]

df_melted -> df_melted_paired

# df_melted_paired$celltype <- factor(df_melted_paired$celltype,
#     levels = c(
#         "M1 macrophages", "M2 macrophages", "Microglia", "Tumor Associated Macrophages",
#         "Neutrophils", "Macrophages", "Lymphoid", "Proliferating Cells", "Unknown", "Neuronal Cells"
#     )
# )

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    # t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2

bxp <- ggboxplot(
    df_melted_paired,
    x = "celltype", y = "percentage",
    color = "timepts", palette = c("#00AFBB", "#E7B800")
)

bxp2 <- bxp + geom_dotplot(
    aes(fill = timepts, color = timepts),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.grid.major = element_line(colour = "white")
)

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"
pdf(paste(savedir, "Table/t_test_all_celltypes_primary_recurrent_bxplot_one_sided_macrophage_subtypes.pdf", sep = ""), width = 5, height = 4)
bxp4
dev.off()

### Dotplot for the markers
M1_genes <- read.table("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/M1_genes.txt")[, 1]
M2_genes <- read.table("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/M2_genes.txt")[, 1]
genes <- c(M1_genes, M2_genes)

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/M1_M2_seurat_clus_RNA.pdf", sep = ""), width = 9, height = 10)
DotPlot(obj, features = genes, group.by = "seurat_clusters", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

### Since I am getting the GCP cells also in the Immune cells. Need to check why they are coming.
p <- VlnPlot(immune, features = c("nCount_RNA", "nFeature_RNA", "doublet_score"), group.by = "subcelltypes")
dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/VlnPlot/")
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/immune_featureplot_subcelltypes.pdf")
p
dev.off()

### Making a gene correlation and graph
source("/diazlab/data3/.abhinav/resources/all_scripts/R/gene_corr_and_graph.R")
p <- graph_correaltion(obj, t(obj@assays$MAGIC_RNA@data), savedir, "CD86", "CD209", "immune")

### Reference mapping to the multimodal data
spec_immune_obj <- SCTransform(spec_immune_obj, verbose = FALSE)
anchors <- FindTransferAnchors(
    reference = reference,
    query = spec_immune_obj,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
)

# We then transfer cell type labels and protein data from the reference to the query. Additionally, we project the query data onto the UMAP structure of the reference.
spec_immune_obj <- MapQuery(
    anchorset = anchors,
    query = spec_immune_obj,
    reference = reference,
    refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
)

#### Generate the prediction from the Azimuth
saveRDS(spec_immune_obj@assays$RNA@counts, paste(savedir, "Table/spec_immune_obj_count.RDS", sep = ""))
az_pred <- read.table(paste(savedir, "Table/azimuth_pred_spec_immune.tsv", sep = ""), sep = "\t", header = TRUE)
stopifnot(all(rownames(spec_immune_obj@meta.data) == az_pred$cell))
spec_immune_obj@meta.data$celltype_pbmc_predicted_spec <- az_pred$predicted.celltype.l2

pdf(paste0(savedir, "UMAP/sampleharmony_dinmplot_celltypes_spec.pdf"), width = 5.5, height = 5.5)
DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", label = TRUE, group.by = "celltype_pbmc_predicted_spec")
dev.off()

#### Immune cell with GCP has more clear marker in the differential adding that
clus_CT <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/Table/clus_celltype_markers.txt",
    header = TRUE,
    sep = "\t"
)

patterns <- clus_CT$Clus
replacements <- clus_CT$Celltypes

spec_immune_obj@meta.data$seurat_clusters -> spec_immune_obj@meta.data$subcelltypes
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    spec_immune_obj@meta.data$subcelltypes <- str_replace_all(spec_immune_obj@meta.data$subcelltypes, pattern, replacements[i])
}

pdf(paste(savedir, "UMAP/celltypes.pdf", sep = ""))
DimPlot(spec_immune_obj, reduction = "sampleharmonyumap", group.by = "subcelltypes", label = TRUE, label.size = 5) + NoLegend()
dev.off()

#### Longitudinal analysis
celltypes_ind <- table(spec_immune_obj@meta.data$SampleID, spec_immune_obj@meta.data$subcelltypes)

write.table(celltypes_ind,
    paste(savedir, "Table/celltype_individual_2.txt", sep = ""),
    row.names = T,
    col.names = T,
    sep = "\t",
    quote = F
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
celltypes_ind <- read.table(paste(savedir, "Table/celltype_individual_2.txt", sep = ""), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

### Sanity Check
rowSums(df)

### Making some change in the naming of the samples
rownames(df) <- gsub("-", "_", rownames(df))

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample", "celltype", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c(
    "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
    "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
    "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994", "SF8368", "SF8539", "SF9232"
) -> sample_levels

df_melted$sample <- gsub("-", "_", df_melted$sample) %>% gsub("SF7994R", "SF7994", .)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = sample)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste(savedir, "Table/celltype_individual_barplot_2.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.delim("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
replacements <- long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]

df_melted -> df_melted_paired
df_melted_paired$celltype <- factor(df_melted_paired$celltype,
    levels = c(
        "M1 macrophages", "M2 macrophages", "Microglia", "Tumor Associated Macrophages",
        "Neutrophils", "Macrophages", "Lymphoid", "Proliferating Cells", "Unknown", "Neuronal Cells"
    )
)

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    # t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2

bxp <- ggboxplot(
    df_melted_paired,
    x = "celltype", y = "percentage",
    color = "timepts", palette = c("#00AFBB", "#E7B800")
)

bxp2 <- bxp + geom_dotplot(
    aes(fill = timepts, color = timepts),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.grid.major = element_line(colour = "white")
)

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"
pdf(paste(savedir, "Table/t_test_all_celltypes_primary_recurrent_bxplot_one_sided.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()

### barplot
gbr <- ggbarplot(df_melted_paired,
    x = "celltype", y = "percentage",
    add = c("mean_se", "jitter"),
    color = "timepts", palette = c("#00AFBB", "#E7B800"),
    position = position_dodge(0.8)
)

stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
gbr2 <- gbr + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")
    )

pdf(paste(savedir, "Table/t_test_celltypes_paired_primary_recurrent_barplot_t_test.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]
gbr <- ggbarplot(df_melted_paired,
    x = "celltype", y = "percentage",
    add = c("mean_se", "jitter"),
    color = "timepts", palette = c("#00AFBB", "#E7B800"),
    position = position_dodge(0.8)
)

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    # t_test(percentage ~ timepts, paired = TRUE, alternative = "greater") %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2


stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
gbr2 <- gbr + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")
    )

pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_barplot_t_test.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

## FeaturePlot markers
mps_marker <- c("CD14", "CD163", "CD68", "CCR5", "CD19", "FCRL2", "TMEM179")

rm(plot_list)
plot_list <- list()
DefaultAssay(spec_immune_obj) <- "MAGIC_RNA"
for (i in 1:length(mps_marker)) {
    plot_list[[i]] <- FeaturePlot(spec_immune_obj,
        mps_marker[i],
        reduction = "sampleharmonyumap"
    ) + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
}

dir.create(paste(savedir, "featureplot", sep = ""))
pdf(paste(savedir, "featureplot/mps_markers.pdf", sep = ""))
plot_list
dev.off()

#### Immune cell with GCP has more clear marker in the differential adding that
clus_CT <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/Table/clus_celltype.txt",
    header = TRUE,
    sep = "\t"
)

patterns <- clus_CT$Clus
replacements <- clus_CT$Celltype

immune@meta.data$seurat_clusters -> immune@meta.data$subcelltypes
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    immune@meta.data$subcelltypes <- str_replace_all(immune@meta.data$subcelltypes, pattern, replacements[i])
}

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/UMAP/subcelltype_dimplot.pdf")
DimPlot(immune, group.by = "subcelltypes", reduction = "sampleharmonyumap")
dev.off()

### Adding to the removed GCP cells object
spec_immune_obj@meta.data$subcelltypes <- immune@meta.data[match(rownames(spec_immune_obj@meta.data), rownames(immune@meta.data), nomatch = 0), "subcelltypes"]

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/UMAP/subcelltype_dimplot.pdf")
DimPlot(spec_immune_obj, group.by = "subcelltypes", reduction = "sampleharmonyumap", label = TRUE)
dev.off()

#### Longitudinal analysis
celltypes_ind <- table(spec_immune_obj@meta.data$SampleID, spec_immune_obj@meta.data$subcelltypes)

write.table(celltypes_ind,
    paste(savedir, "Table/subcelltype_individual.txt", sep = ""),
    row.names = T,
    col.names = T,
    sep = "\t",
    quote = F
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/"
celltypes_ind <- read.table(paste(savedir, "Table/subcelltype_individual.txt", sep = ""), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

### Sanity Check
rowSums(df)

### Making some change in the naming of the samples
rownames(df) <- gsub("-", "_", rownames(df))

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample", "celltype", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c(
    "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
    "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
    "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994", "SF8368", "SF8539", "SF9232"
) -> sample_levels

df_melted$sample <- gsub("-", "_", df_melted$sample) %>% gsub("SF7994R", "SF7994", .)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = sample)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste(savedir, "Table/subcelltype_individual_barplot.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.delim("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
replacements <- long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]

df_melted -> df_melted_paired

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    # t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2

bxp <- ggboxplot(
    df_melted_paired,
    x = "celltype", y = "percentage",
    color = "timepts", palette = c("#00AFBB", "#E7B800")
)

bxp2 <- bxp + geom_dotplot(
    aes(fill = timepts, color = timepts),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.grid.major = element_line(colour = "white")
)


# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"
pdf(paste(savedir, "Table/t_test_all_primary_recurrent_bxplot_one_sided_subcelltypes.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()

### barplot
gbr <- ggbarplot(df_melted,
    x = "celltype", y = "percentage",
    add = c("mean_se", "jitter"),
    color = "timepts", palette = c("#00AFBB", "#E7B800"),
    position = position_dodge(0.8)
)

stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
gbr2 <- gbr + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")
    )

pdf(paste(savedir, "Table/t_test_all_primary_recurrent_barplot_t_test_subcelltype.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]
gbr <- ggbarplot(df_melted_paired,
    x = "celltype", y = "percentage",
    add = c("mean_se", "jitter"),
    color = "timepts", palette = c("#00AFBB", "#E7B800"),
    position = position_dodge(0.8)
)

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    # t_test(percentage ~ timepts, paired = TRUE, alternative = "greater") %>%
    t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test$p <- stat.test$p / 2


stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
gbr2 <- gbr + stat_pvalue_manual(
    stat.test,
    label = "p", tip.length = 0
) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")
    )

pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_barplot_t_test.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

### Extracting out the Myeloid cells to perform the differential
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"
myeloid_cells <- rownames(obj@meta.data[grep("M1_Microglia|M0_Macrophage|M1_Macrophage|M2_Macrophage|M0_TAMs", obj@meta.data$subcelltypes2, ignore.case = TRUE), ])
myeloid_obj <- subset(spec_immune_obj, cells = myeloid_cells)

DefaultAssay(myeloid_obj) <- "RNA"
myeloid_obj <- NormalizeData(myeloid_obj, normalization.method = "LogNormalize", scale.factor = 10000)
myeloid_obj <- ScaleData(myeloid_obj, features = rownames(myeloid_obj))

# primary and recurrent samples
patientid <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")

for (i in 1:length(patientid)) {
    group1 <- paste(patientid[i], "_Recurrent", sep = "")
    group2 <- paste(patientid[i], "_Primary", sep = "")
    markers <- FindMarkers(myeloid_obj, ident.1 = group1, ident.2 = group2, group.by = "patient_tumortype", test.use = "MAST", min.pct = 0, logfc.threshold = 0)
    assign(paste(patientid[i], "markers_MAST", sep = "_"), markers)
    write.table(markers, paste(savedir, "table/", patientid[i], "_recurrent_vs_primary_markers_MAST_nocutoff.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

#### Identifying the common genes in primary and recurrent samples
#### Making an upset plot for the GCP markers
### Making an Upset plot ####
markers2 <- ls(pattern = "_markers_MAST")

### First  putting the cutoff for padj < 0.05 and logFC > 0.5
for (i in 1:length(markers2)) {
    diff_table <- get(markers2[i])
    diff_table_sig <- diff_table[diff_table$p_val < 0.05 & diff_table$avg_log2FC > 0, ]
    assign(paste(markers[i], "_sig", sep = ""), diff_table_sig)
}

markers <- ls(pattern = "_markers_MAST_sig")
library(data.table)
PT_7WYPEC3Q <- rownames(get(markers[1]))
PT_9S6WMQ92 <- rownames(get(markers[2]))
PT_H3WWDMW9 <- rownames(get(markers[3]))
PT_XA98HG1C <- rownames(get(markers[4]))

n <- max(
    length(PT_7WYPEC3Q), length(PT_9S6WMQ92),
    length(PT_H3WWDMW9), length(PT_XA98HG1C)
)

length(PT_7WYPEC3Q) <- n
length(PT_9S6WMQ92) <- n
length(PT_H3WWDMW9) <- n
length(PT_XA98HG1C) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q), (PT_9S6WMQ92),
    (PT_H3WWDMW9), (PT_XA98HG1C)
))

library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_markers_MAST_sig", "", markers)

dir.create(paste(savedir, "plots", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "plots/GCP_marker_upset_plot.pdf", sep = ""), width = 12, height = 8)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q, "PT_9S6WMQ92" = PT_9S6WMQ92,
    "PT_H3WWDMW9" = PT_H3WWDMW9, "PT_XA98HG1C" = PT_XA98HG1C
)

df2 <- data.frame(gene = unique(unlist(list_filter)))

df1 <- lapply(list_filter, function(x) {
    data.frame(gene = x)
}) %>%
    bind_rows(.id = "path")

df_int <- lapply(df2$gene, function(x) {
    # pull the name of the intersections
    intersection <- df1 %>%
        dplyr::filter(gene == x) %>%
        arrange(path) %>%
        pull("path") %>%
        paste0(collapse = "|")

    # build the dataframe
    data.frame(gene = x, int = intersection)
}) %>%
    bind_rows()

library(dplyr)
final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

df_int$int2 <- gsub("\\|", "__", df_int$int)
# GCP_common <- df_int[grep("mar_7316_2118__mar_7316_278__mar_7316_2978__mar_7316_3023__mar_7316_311__mar_7316_333__mar_7316_4529__mar_7316_5881__mar_7316_737__mar_7316_931__mar_DOD4182__mar_SF10961__mar_SF12930__mar_SF7994__mar_SF8368__mar_SF8539",df_int$int2),1]

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
recurrent_common <- df_int[(df_int$sample_num > 2), 1]

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(df_int, paste(savedir, "table/Recurrent_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(recurrent_common, paste(savedir, "table/Recurrent_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

#### Dotplot for the specific genes
genes <- c("MAPK1", "AKT3", "PIK3CB", "HIF1A", "VEGFA", "SPP1", "HK2", "PGK1", "EGLN3", "EGLN1", "SLC2A3", "GRB2", "SOS1", "MDM2", "GNB4")

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/PI3K_AKT_pathway_recurrent_RNA.pdf", sep = ""))
DotPlot(obj, genes, group.by = "subcelltypes2", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

myeloid_rec_vs_pri <- FindMarkers(myeloid_obj, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type", test.use = "MAST")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(myeloid_rec_vs_pri, paste(savedir, "table/myeloid_rec_vs_prim_MAST.txt", sep = ""), quote = F, row.names = T, sep = "\t", col.names = T)

### Primary
#### Identifying the common genes in primary and recurrent samples
#### Making an upset plot for the GCP markers
### Making an Upset plot ####
markers2 <- ls(pattern = "_markers_MAST$")

### First  putting the cutoff for padj < 0.05 and logFC > 0.5
for (i in 1:length(markers2)) {
    diff_table <- get(markers2[i])
    diff_table_sig <- diff_table[diff_table$p_val < 0.05 & diff_table$avg_log2FC < 0, ]
    assign(paste(markers2[i], "_primary_sig", sep = ""), diff_table_sig)
}

markers <- ls(pattern = "_markers_MAST_primary_sig")
library(data.table)
PT_7WYPEC3Q <- rownames(get(markers[1]))
PT_9S6WMQ92 <- rownames(get(markers[2]))
PT_H3WWDMW9 <- rownames(get(markers[3]))
PT_XA98HG1C <- rownames(get(markers[4]))

n <- max(
    length(PT_7WYPEC3Q), length(PT_9S6WMQ92),
    length(PT_H3WWDMW9), length(PT_XA98HG1C)
)

length(PT_7WYPEC3Q) <- n
length(PT_9S6WMQ92) <- n
length(PT_H3WWDMW9) <- n
length(PT_XA98HG1C) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q), (PT_9S6WMQ92),
    (PT_H3WWDMW9), (PT_XA98HG1C)
))

library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_markers_MAST_primary_sig", "", markers)

dir.create(paste(savedir, "plots", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "plots/primary_marker_upset_plot.pdf", sep = ""), width = 12, height = 8)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q, "PT_9S6WMQ92" = PT_9S6WMQ92,
    "PT_H3WWDMW9" = PT_H3WWDMW9, "PT_XA98HG1C" = PT_XA98HG1C
)

df2 <- data.frame(gene = unique(unlist(list_filter)))

df1 <- lapply(list_filter, function(x) {
    data.frame(gene = x)
}) %>%
    bind_rows(.id = "path")

df_int <- lapply(df2$gene, function(x) {
    # pull the name of the intersections
    intersection <- df1 %>%
        dplyr::filter(gene == x) %>%
        arrange(path) %>%
        pull("path") %>%
        paste0(collapse = "|")

    # build the dataframe
    data.frame(gene = x, int = intersection)
}) %>%
    bind_rows()

library(dplyr)
final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

df_int$int2 <- gsub("\\|", "__", df_int$int)
# GCP_common <- df_int[grep("mar_7316_2118__mar_7316_278__mar_7316_2978__mar_7316_3023__mar_7316_311__mar_7316_333__mar_7316_4529__mar_7316_5881__mar_7316_737__mar_7316_931__mar_DOD4182__mar_SF10961__mar_SF12930__mar_SF7994__mar_SF8368__mar_SF8539",df_int$int2),1]

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
recurrent_common <- df_int[(df_int$sample_num > 2), 1]

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(df_int, paste(savedir, "table/Primary_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(recurrent_common, paste(savedir, "table/Primary_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

#### Dotplot for the specific genes
genes <- c("MAPK1", "AKT3", "PIK3CB", "HIF1A", "VEGFA", "SPP1", "HK2", "PGK1", "EGLN3", "EGLN1", "SLC2A3", "GRB2", "SOS1", "MDM2", "GNB4")

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/PI3K_AKT_pathway_recurrent_RNA.pdf", sep = ""))
DotPlot(obj, genes, group.by = "subcelltypes2", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

myeloid_rec_vs_pri <- FindMarkers(myeloid_obj, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type", test.use = "MAST")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(myeloid_rec_vs_pri, paste(savedir, "table/myeloid_rec_vs_prim_MAST.txt", sep = ""), quote = F, row.names = T, sep = "\t", col.names = T)


myeloid_rec_vs_pri <- FindMarkers(myeloid_obj, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"

dir.create(paste(savedir, "table", sep = ""), showWarnings = FALSE)
write.table(myeloid_rec_vs_pri, paste(savedir, "table/myeloid_rec_vs_prim.txt", sep = ""), quote = F, row.names = T, sep = "\t", col.names = T)

# genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/Table/recurrent_important_genes", header = FALSE)[, 1]

myeloid_paired_cells <- rownames(myeloid_obj@meta.data[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", myeloid_obj@meta.data$sample2), ])

myeloid_obj@meta.data$patient_tumortype <- paste(myeloid_obj@meta.data$patient, myeloid_obj@meta.data$tumor_type, sep = "_")
myeloid_obj_paired <- subset(myeloid_obj, cells = myeloid_paired_cells)
DefaultAssay(myeloid_obj_paired) <- "RNA"
myeloid_obj_paired <- NormalizeData(myeloid_obj_paired)
genes <- c("EGLN1", "EGLN3", "HIF1A", "VEGFA", "PGK1", "BNIP3")
genes <- c("HK2", "SLC2A1", "SLC2A3", "SPP1", "MMP9")
genes <- c("B2M", "HLA-DRA", "HLA-DPB1")
genes <- c("CD163", "CTSB")

rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(myeloid_obj_paired, genes[i],
        pt.size = 0,
        group.by = "patient_tumortype",
        assay = "MAGIC_RNA",
        cols = c("lightblue3", "lightcoral", "lightblue3", "lightcoral", "lightblue3", "lightcoral", "lightblue3", "lightcoral")
    ) + geom_boxplot()
}

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/macrophage_activation_genes_paired_2.pdf", sep = ""))
plot_list
dev.off()
