library(Seurat)
library(dplyr)
library(UCell)
spec_immune_obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/saveRDS_obj/immune_clean.RDS")
mac_clus <- c("0", "3", "8", "6", "2", "5", "11") ### based on the cell markers and add MOdule score

myeloid_cells <- rownames(spec_immune_obj@meta.data[grep(paste("^", mac_clus, "$", sep = "", collapse = "|"), spec_immune_obj@meta.data$seurat_clusters), ])

myeloid_obj <- subset(spec_immune_obj, cells = myeloid_cells)

DefaultAssay(myeloid_obj) <- "RNA"
myeloid_obj <- NormalizeData(myeloid_obj)

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
myeloid_obj <- AddModuleScore_UCell(myeloid_obj, features = markers, ncores = 4, slot = "data", assay = "RNA")
featurename <- c("M1_Abhinav_UCell", "M2_Abhinav_UCell")

p1 <- VlnPlot(myeloid_obj, features = featurename, group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(myeloid_obj, features = featurename[i], reduction = "sampleharmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/Abhinav_M1_M2_data_RNA.pdf", sep = ""), width = 8, height = 5)
print(p1)
dev.off()

dir.create(paste(savedir, "featureplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/Abhinav_M1_M2_RNA.pdf", sep = ""))
print(plot_list)
dev.off()

names(markers) <- paste(names(markers), "_imputed", sep = "")
myeloid_obj <- AddModuleScore_UCell(myeloid_obj, features = markers, ncores = 4, slot = "data", assay = "MAGIC_RNA")
featurename <- c("M1_Abhinav_imputed_UCell", "M2_Abhinav_imputed_UCell")

p1 <- VlnPlot(myeloid_obj, features = featurename, group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(myeloid_obj, features = featurename[i], reduction = "sampleharmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/Abhinav_M1_M2_data_MAGIC_RNA.pdf", sep = ""), width = 8, height = 5)
print(p1)
dev.off()

dir.create(paste(savedir, "featureplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/Abhinav_M1_M2_MAGIC_RNA.pdf", sep = ""))
print(plot_list)
dev.off()

### Using Seurat AddModule Score
DefaultAssay(myeloid_obj) <- "MAGIC_RNA"
myeloid_obj <- AddModuleScore(myeloid_obj, features = markers, slot = "data")
featurename <- c("M1_Abhinav", "M2_Abhinav")
colnames(myeloid_obj@meta.data)[69:70] <- paste(featurename, "_imputed_seurat", sep = "")

p1 <- VlnPlot(myeloid_obj, features = paste(featurename, "_imputed_seurat", sep = ""), group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(myeloid_obj,
        features = paste(featurename, "_imputed_seurat", sep = "")[i],
        reduction = "sampleharmonyumap"
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/Abhinav_M1_M2_Ucell_Addmod_imputed_seurat.pdf", sep = ""), width = 15, height = 11)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/Abhinav_M1_M2_Ucell_imputed_Addmod_2.pdf", sep = ""))
print(plot_list)
dev.off()

### M0 M1 and M2
### Now for the data score for AddModuleScore making the M1 and M2
myeloid_obj@meta.data$macrophage_subtype <- "M0"
myeloid_obj@meta.data[myeloid_obj@meta.data$M1_Abhinav_imputed_seurat > myeloid_obj@meta.data$M2_Abhinav_imputed_seurat, "macrophage_subtype"] <- "M1"
myeloid_obj@meta.data[myeloid_obj@meta.data$M1_Abhinav_imputed_seurat < myeloid_obj@meta.data$M2_Abhinav_imputed_seurat, "macrophage_subtype"] <- "M2"
myeloid_obj@meta.data[myeloid_obj@meta.data$M1_Abhinav_imputed_seurat < 0 & myeloid_obj@meta.data$M2_Abhinav_imputed_seurat < 0, "macrophage_subtype"] <- "M0"

pdf(paste(savedir, "UMAP/macrophage_subtype.pdf", sep = ""))
DimPlot(myeloid_obj, group.by = "macrophage_subtype", reduction = "sampleharmonyumap") + scale_color_manual(values = c(
    "#279e68", "#d62728", "orange", "#1f77b4", "#aa40fc",
    "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white"
))
dev.off()

pdf(paste(savedir, "UMAP/seurat_cluster.pdf", sep = ""))
DimPlot(myeloid_obj, reduction = "sampleharmonyumap", label = TRUE)
dev.off()

## Finding the cluster distribution
seurat_macrophage <- table(myeloid_obj@meta.data$seurat_clusters, myeloid_obj@meta.data$macrophage_subtype)
seurat_macrophage_req <- seurat_macrophage[rowSums(seurat_macrophage) != 0, ]

dir.create(paste(savedir, "Table", sep = ""), showWarnings = FALSE)
write.table(seurat_macrophage_req,
    paste(savedir, "Table/seurat_macrophages.txt", sep = ""),
    quote = F, col.names = T, row.names = T, sep = "\t"
)

seurat_macrophage <- read.table("//diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/Table/seurat_macrophages.txt", header = TRUE, row.names = 1)

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

pdf(paste(savedir, "UMAP/dimplot.pdf", sep = ""))
DimPlot(myeloid_obj, label = TRUE, reduction = "sampleharmonyumap")
dev.off()

myeloid_obj@meta.data$SampleID <- gsub("-", "_", myeloid_obj@meta.data$SampleID)
celltypes_ind <- table(myeloid_obj@meta.data$SampleID, myeloid_obj@meta.data$macrophage_subtype)

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

# df_melted -> df_melted_paired

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    t_test(percentage ~ timepts, paired = TRUE, alternative = "two.sided") %>%
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
pdf(paste(savedir, "Table/t_Test_macrophage_subtypes_paired.pdf", sep = ""), width = 5, height = 4)
bxp4
dev.off()

### Dotplot for the markers
M1_genes <- read.table("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/M1_genes.txt")[, 1]
M2_genes <- read.table("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/M2_genes.txt")[, 1]
genes <- unique(c(M1_genes, M2_genes))

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/M1_M2_seurat_clus_RNA.pdf", sep = ""), width = 9, height = 10)
DotPlot(myeloid_obj, features = genes, group.by = "seurat_clusters", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

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

pdf(paste(savedir, "Table/celltype_individual_M0_M1.pdf", sep = ""), width = 10, height = 7)
p
dev.off()


### Microglia
library(Seurat)
library(dplyr)
library(UCell)
spec_immune_obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/saveRDS_obj/immune_clean.RDS")
mic_clus <- c("4") ### based on the cell markers and add MOdule score

micro_cells <- rownames(spec_immune_obj@meta.data[grep(paste("^", mic_clus, "$", sep = "", collapse = "|"), spec_immune_obj@meta.data$seurat_clusters), ])
microglia_obj <- subset(spec_immune_obj, cells = micro_cells)

DefaultAssay(microglia_obj) <- "RNA"
microglia_obj <- NormalizeData(microglia_obj)

DefaultAssay(microglia_obj) <- "MAGIC_RNA"
microglia_obj <- AddModuleScore(microglia_obj, features = markers, slot = "data")
featurename <- c("M1_Abhinav", "M2_Abhinav")
colnames(microglia_obj@meta.data)[63:64] <- paste(featurename, "_imputed_seurat", sep = "")

p1 <- VlnPlot(microglia_obj, features = paste(featurename, "_imputed_seurat", sep = ""), group.by = "seurat_clusters")

rm(plot_list)
plot_list <- list()
for (i in 1:length(featurename)) {
    plot_list[[i]] <- FeaturePlot(microglia_obj,
        features = paste(featurename, "_imputed_seurat", sep = "")[i],
        reduction = "sampleharmonyumap"
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/"
dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/Abhinav_M1_M2_Ucell_microglia.pdf", sep = ""), width = 15, height = 11)
print(p1)
dev.off()

pdf(paste(savedir, "featureplot/Abhinav_M1_M2_Ucell_imputed_Addmod_2.pdf", sep = ""))
print(plot_list)
dev.off()

### M0 M1 and M2
### Now for the data score for AddModuleScore making the M1 and M2
microglia_obj@meta.data$microglia_subtype <- "M0"
microglia_obj@meta.data[microglia_obj@meta.data$M1_Abhinav_imputed_seurat > microglia_obj@meta.data$M2_Abhinav_imputed_seurat, "microglia_subtype"] <- "M1"
microglia_obj@meta.data[microglia_obj@meta.data$M1_Abhinav_imputed_seurat < microglia_obj@meta.data$M2_Abhinav_imputed_seurat, "microglia_subtype"] <- "M2"
microglia_obj@meta.data[microglia_obj@meta.data$M1_Abhinav_imputed_seurat < 0 & microglia_obj@meta.data$M2_Abhinav_imputed_seurat < 0, "microglia_subtype"] <- "M0"

pdf(paste(savedir, "UMAP/microglia_subtype.pdf", sep = ""))
DimPlot(microglia_obj, group.by = "microglia_subtype", reduction = "sampleharmonyumap") + scale_color_manual(values = c(
    "#279e68", "#d62728", "orange", "#1f77b4", "#aa40fc",
    "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white"
))
dev.off()

pdf(paste(savedir, "UMAP/microglia_seurat_cluster.pdf", sep = ""))
DimPlot(microglia_obj, reduction = "sampleharmonyumap", label = TRUE)
dev.off()

## Finding the cluster distribution
seurat_macrophage <- table(microglia_obj@meta.data$seurat_clusters, microglia_obj@meta.data$microglia_subtype)
seurat_macrophage_req <- seurat_macrophage[rowSums(seurat_macrophage) != 0, ]

dir.create(paste(savedir, "Table", sep = ""), showWarnings = FALSE)
write.table(seurat_macrophage_req,
    paste(savedir, "Table/seurat_microglia.txt", sep = ""),
    quote = F, col.names = T, row.names = T, sep = "\t"
)

seurat_microglia <- read.table("//diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/remove_GCP_cells/myeloid_cells/Table/seurat_microglia.txt", header = TRUE, row.names = 1)

n_cells <- t(seurat_microglia)
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

pdf(paste(savedir, "Table/microglia_seurat_clusters_2.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

pdf(paste(savedir, "UMAP/dimplot.pdf", sep = ""))
DimPlot(microglia_obj, label = TRUE, reduction = "sampleharmonyumap")
dev.off()

microglia_obj@meta.data$SampleID <- gsub("-", "_", microglia_obj@meta.data$SampleID)
celltypes_ind <- table(microglia_obj@meta.data$SampleID, microglia_obj@meta.data$microglia_subtype)

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

# df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]

df_melted -> df_melted_paired

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    # t_test(percentage ~ timepts, paired = TRUE, alternative = "two.sided") %>%
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

# savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/immune_cells/"
pdf(paste(savedir, "Table/t_Test_microglia_subtypes_paired.pdf", sep = ""), width = 5, height = 4)
bxp4
dev.off()

### Dotplot for the markers
M1_genes <- read.table("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/M1_genes.txt")[, 1]
M2_genes <- read.table("/diazlab/data3/.abhinav/resources/gene_list/immune_cells/Abhinav_M1_M2/M2_genes.txt")[, 1]
genes <- unique(c(M1_genes, M2_genes))

dir.create(paste(savedir, "dotplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "dotplot/M1_M2_seurat_clus_RNA.pdf", sep = ""), width = 9, height = 10)
DotPlot(microglia_obj, features = genes, group.by = "seurat_clusters", assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() + scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

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

pdf(paste(savedir, "Table/celltype_individual_M0_M1_microglia.pdf", sep = ""), width = 10, height = 7)
p
dev.off()
