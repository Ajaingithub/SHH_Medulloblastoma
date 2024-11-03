#### Analysis for TP53 snRNA
library(stringr)
library(ArchR)
library(Seurat)

sampleid_TP53 <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/sampleid_TP53.txt", header = TRUE)

patterns <- sampleid_TP53$Sample
replacements <- sampleid_TP53$TP53_status
snRNA@meta.data$TP53_status <- snRNA@meta.data$SampleID

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    snRNA@meta.data$TP53_status <- str_replace_all(snRNA@meta.data$TP53_status, pattern, replacements[i])
}

#### Based on the seurat_cluster
# objsubset@meta.data$celltypes <- gsub("Neuronal","GN_intermediate",objsubset@meta.data$celltypes)
snRNA@meta.data$broad_CT <- snRNA@meta.data$celltypes

patterns <- c(
    "Oligo", "GCP", "Astro", "Pericyte", "OPC", "Neuron_other", "GN_early",
    "GCP_cycling", "GN_late", "GN_cycling", "Immune", "RL_like", "Endo"
)

replacements <- c(
    "non-malign", "GCP", "non-malign", "non-malign", "non-malign", "non-malign", "GN",
    "GCP", "GN", "GN", "non-malign", "non-malign", "non-malign"
)

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    snRNA@meta.data$broad_CT <- str_replace_all(snRNA@meta.data$broad_CT, pattern, replacements[i])
}


snRNA@meta.data$sample_TP53 <- paste(snRNA@meta.data$sample, snRNA@meta.data$TP53_status, sep = "_")
celltypes_ind <- table(snRNA@meta.data$sample_TP53, snRNA@meta.data$broad_CT)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltypes_sample_TP53.txt"),
    row.names = T, col.names = T, sep = "\t", quote = F
)

n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample", "celltype", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
# c(
#     "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
#     "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
#     "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994",
#     "SF8368", "SF8539", "SF9232"
# ) -> sample_levels

# df_melted$sample <- gsub("-", "_", df_melted$sample)
# df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = sample)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white", "yellow", "green", "red", "grey23", "grey34"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/TP53/Table/celltype_broad_TP53mutated_2.pdf", width = 6, height = 5)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

# long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

# patterns <- long_sample$Sample.ID
# patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
# replacements <- long_sample$Tumor.type
# # names(replacement) <- patterns
# df_melted$timepts <- df_melted$sample
# for (i in seq_along(patterns)) {
#     pattern <- paste0("\\b", patterns[i], "\\b")
#     df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
# }

df_melted$TP53_status <- gsub(".*._", "", df_melted$sample)

df_melted_paired <- df_melted
# df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]
# df_melted_paired <- df_melted

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    t_test(percentage ~ TP53_status, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test
stat.test$p_one_sided <- (stat.test$p) / 2
stat_test_df <- as.data.frame(stat.test)

bxp <- ggboxplot(
    df_melted_paired,
    x = "celltype", y = "percentage",
    color = "TP53_status", palette = c("deepskyblue4", "firebrick1")
)

bxp2 <- bxp + geom_dotplot(
    aes(fill = TP53_status, color = TP53_status),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("deepskyblue4", "firebrick1")) +
    scale_color_manual(values = c("deepskyblue4", "firebrick1"))


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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/TP53/"
dir.create(paste0(savedir, "Table/"), showWarnings = FALSE)
pdf(paste(savedir, "Table/t_Test_one_sided_TP53_broad_CT.pdf", sep = ""), width = 5, height = 4)
bxp4
dev.off()


#### Differential Genes b/w TP53wt vs TP53 mutants
FindMarkers(snRNA, ident.1 = "TP53mut", ident.2 = "TP53wt", group.by = "TP53_status") -> mut_vs_wt

#### Performing Pseudobulk Differential
obj <- snRNA
cluster <- "all"
group1 <- "TP53mut"
group2 <- "TP53wt"
DefaultAssay(obj) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/diazlab/data3/.abhinav/resources/all_scripts/R/pseudobulk_PCA_within_no_batch_correction_AJ.R")
SHH_dds <- pseudobulk_within_cluster_AJ_no_batch_correction(
    obj = obj, savedir = savedir,
    group1 = group1, group2 = group2,
    grouping_by = "TP53_status",
    cluster = cluster,
    cell_freq = 20,
    remove_samples = remove_samples,
    cluster_group = "celltypes",
    sample_col = "orig.ident",
    batch_col = "Run",
    gene_min_counts = 5,
    column_name_split = c("sampleID", "Age", "Sex", "tumor_type", "Run"),
    splitting = "-"
)

savedir2 <- paste(savedir, "pseudobulk/clus_", paste(cluster, collapse = "_"), "_removed_", remove_samples_name, "_", group1, "_vs_", group2, "/", sep = "")
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
source("/diazlab/data3/.abhinav/resources/all_scripts/R/RNAseq_limma_EdgeR.R")
design0 <- model.matrix(~ 0 + TP53_status, data = colData(SHH_dds))
colnames(design0) <- c(group1, group2)
cm <- makeContrasts(TP53mut_VS_TP53wt = TP53mut - TP53wt, levels = design0)
desl_clus <- LimmaEdgeR_differential(
    dds = SHH_dds,
    design0 = design0,
    cm = cm,
    savedir = savedir2,
    logfc = 0.25,
    p_value_adj = 0.05,
    grouping_by = "TP53_status",
    column_name_split = c("sampleID", "Age", "Sex", "tumor_type", "Run")
)

#### Performing the similar analysis for scATAC
scATAC <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_chromvar_impute_rna_linkpeaks.RDS")

sampleid_TP53 <- read.table("/diazlab/data3/.abhinav/projects/SHH/Genotype/sampleid_TP53.txt", header = TRUE)

patterns <- sampleid_TP53$Sample
replacements <- sampleid_TP53$TP53_status
snRNA@meta.data$TP53_status <- snRNA@meta.data$SampleID

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    snRNA@meta.data$TP53_status <- str_replace_all(snRNA@meta.data$TP53_status, pattern, replacements[i])
}

#### Based on the seurat_cluster
# objsubset@meta.data$celltypes <- gsub("Neuronal","GN_intermediate",objsubset@meta.data$celltypes)
snRNA@meta.data$broad_CT <- snRNA@meta.data$celltypes

patterns <- c(
    "Oligo", "GCP", "Astro", "Pericyte", "OPC", "Neuron_other", "GN_early",
    "GCP_cycling", "GN_late", "GN_cycling", "Immune", "RL_like", "Endo"
)

replacements <- c(
    "non-malign", "GCP", "non-malign", "non-malign", "non-malign", "non-malign", "GN",
    "GCP", "GN", "GN", "non-malign", "non-malign", "non-malign"
)

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    snRNA@meta.data$broad_CT <- str_replace_all(snRNA@meta.data$broad_CT, pattern, replacements[i])
}


snRNA@meta.data$sample_TP53 <- paste(snRNA@meta.data$sample, snRNA@meta.data$TP53_status, sep = "_")
celltypes_ind <- table(snRNA@meta.data$sample_TP53, snRNA@meta.data$broad_CT)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltypes_sample_TP53.txt"),
    row.names = T, col.names = T, sep = "\t", quote = F
)

n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
    df[j, ] <- (n_cells[j, ] / n_cells_sum[j]) * 100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample", "celltype", "percentage")
library(ggplot2)


# df_melted$sample <- gsub("-", "_", df_melted$sample)
# df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = sample)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
        "yellow", "black", "blanchedalmond", "blue", "white", "yellow", "green", "red", "grey23", "grey34"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/TP53/Table/celltype_broad_TP53mutated_2.pdf", width = 6, height = 5)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

# long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

# patterns <- long_sample$Sample.ID
# patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
# replacements <- long_sample$Tumor.type
# # names(replacement) <- patterns
# df_melted$timepts <- df_melted$sample
# for (i in seq_along(patterns)) {
#     pattern <- paste0("\\b", patterns[i], "\\b")
#     df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
# }

df_melted$TP53_status <- gsub(".*._", "", df_melted$sample)

df_melted_paired <- df_melted
# df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]
# df_melted_paired <- df_melted

stat.test <- df_melted_paired %>%
    group_by(celltype) %>%
    t_test(percentage ~ TP53_status, paired = FALSE, alternative = "two.sided") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
stat.test
stat.test$p_one_sided <- (stat.test$p) / 2
stat_test_df <- as.data.frame(stat.test)

bxp <- ggboxplot(
    df_melted_paired,
    x = "celltype", y = "percentage",
    color = "TP53_status", palette = c("deepskyblue4", "firebrick1")
)

bxp2 <- bxp + geom_dotplot(
    aes(fill = TP53_status, color = TP53_status),
    trim = FALSE,
    binaxis = "y", stackdir = "center", dotsize = 0.15,
    position = position_dodge(0.8)
) +
    scale_fill_manual(values = c("deepskyblue4", "firebrick1")) +
    scale_color_manual(values = c("deepskyblue4", "firebrick1"))


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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/TP53/"
dir.create(paste0(savedir, "Table/"), showWarnings = FALSE)
pdf(paste(savedir, "Table/t_Test_one_sided_TP53_broad_CT.pdf", sep = ""), width = 5, height = 4)
bxp4
dev.off()


#### Differential Genes b/w TP53wt vs TP53 mutants
FindMarkers(snRNA, ident.1 = "TP53mut", ident.2 = "TP53wt", group.by = "TP53_status") -> mut_vs_wt

#### Performing Pseudobulk Differential
obj <- snRNA
cluster <- "all"
group1 <- "TP53mut"
group2 <- "TP53wt"
DefaultAssay(obj) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/diazlab/data3/.abhinav/resources/all_scripts/R/pseudobulk_PCA_within_no_batch_correction_AJ.R")
SHH_dds <- pseudobulk_within_cluster_AJ_no_batch_correction(
    obj = obj, savedir = savedir,
    group1 = group1, group2 = group2,
    grouping_by = "TP53_status",
    cluster = cluster,
    cell_freq = 20,
    remove_samples = remove_samples,
    cluster_group = "celltypes",
    sample_col = "orig.ident",
    batch_col = "Run",
    gene_min_counts = 5,
    column_name_split = c("sampleID", "Age", "Sex", "tumor_type", "Run"),
    splitting = "-"
)

savedir2 <- paste(savedir, "pseudobulk/clus_", paste(cluster, collapse = "_"), "_removed_", remove_samples_name, "_", group1, "_vs_", group2, "/", sep = "")
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
source("/diazlab/data3/.abhinav/resources/all_scripts/R/RNAseq_limma_EdgeR.R")
design0 <- model.matrix(~ 0 + TP53_status, data = colData(SHH_dds))
colnames(design0) <- c(group1, group2)
cm <- makeContrasts(TP53mut_VS_TP53wt = TP53mut - TP53wt, levels = design0)
desl_clus <- LimmaEdgeR_differential(
    dds = SHH_dds,
    design0 = design0,
    cm = cm,
    savedir = savedir2,
    logfc = 0.25,
    p_value_adj = 0.05,
    grouping_by = "TP53_status",
    column_name_split = c("sampleID", "Age", "Sex", "tumor_type", "Run")
)
