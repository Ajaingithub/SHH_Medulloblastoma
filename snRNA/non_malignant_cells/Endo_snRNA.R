### Extracting endothelial cells
library(Seurat)
library(Signac)

Endocells <- rownames(RNA_integrated@meta.data[grep("Endo", RNA_integrated@meta.data$celltypes), ])
Endothelial <- subset(RNA_integrated, cells = Endocells)
DefaultAssay(Endothelial) <- "RNA"

Endothelial <- NormalizeData(Endothelial, normalization.method = "LogNormalize", scale.factor = 10000)
Endothelial <- ScaleData(Endothelial, features = rownames(Endothelial))

Endothelial <- FindVariableFeatures(Endothelial)
Endothelial <- Endothelial %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/EndoDimPlot_patient.pdf")
DimPlot(Endothelial, reduction = "umap", group.by = "patient_tumortype")
dev.off()

endo_rec_vs_pri <- FindMarkers(Endothelial, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type")

#### Performing the differential between primary and Recurrent
rm(plot_list)
plot_list <- list()
genes <- c("LAMA4", "LAMB1", "COL15A1", "COL4A2", "APLN", "SERPINH1", "NID2", "PGF", "COL15A1", "NES")
genes <- c("APLN", "EGLN3", "ANGPT2", "CD34", "KDR", "ITGA2", "AIMP1", "MED1", "APLNR")
genes <- c("HIF1A", "IGFBP7", "LAMA4", "LAMB1", "COL15A1", "COL4A2", "SERPINH1", "NID2", "PGF")
Endothelial@meta.data$patient_tumortype_SampleID <- paste(Endothelial@meta.data$patient_tumortype, Endothelial@meta.data$SampleID, sep = "_")

for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(Endothelial, genes[i], group.by = "tumor_type", assay = "MAGIC_RNA") + geom_boxplot() + NoLegend()
}

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/Endothelial_Recurrent_high_angiogenesis_other_markers.pdf")
plot_list
dev.off()

rm(plot_list)
plot_list <- list()
genes <- c("LAMA4", "LAMB1", "COL15A1", "COL4A2", "APLN", "SERPINH1", "NID2", "PGF", "COL15A1", "NES")
genes <- c("ANGPT2", "EGLN3", "KCNJ8")

genes <- c(
    "DLL4", "JAG1", "NOTCH1", "FLT1", "KDR", "NRP1", "ROBO4", "YAP1", "ESM1", "RAMP3",
    "PECAM1", "CD31", "VWF", "CLDN5", "PROM1", "ANGPT2", "APLN", "APLNR", "PLVAP", "CA2",
    "LAMA4", "LAMB1", "COL15A1", "COL4A2", "SERPINH1", "NID2", "COL15A1", "NES"
)

genes <- c("FLT4", "RAMP3", "NOTCH1", "VCAM1", "ICAM1", "NFKB", "FOXO1", "ANGPT1", "TIE1")

genes <- grep(paste("^", paste(genes, collapse = "$|^", sep = ""), "$", sep = ""), rownames(Endothelial), value = TRUE)

rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(Endothelial, "APLNR", group.by = "tumor_type", assay = "MAGIC_RNA") + geom_boxplot() + NoLegend()
}

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/Endothelial_Recurrent_high_APLNR.pdf")
plot_list
dev.off()

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/VCAM1_Endothelial_cells.pdf")
p
dev.off()


Endo_snRNA <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/saveRDS_Endothelial/Endothelial_snRNA.RDS")

genes <- unique(c(
    "AIMP1", "MED1", "CLIC4", "FLT1", "ANGPT2", "PDCD10", "APLNR", "CCL2", "TGFBI", "FGF1", "APLN",
    "ITGB1", "FLT1", "ANGPT2", "AKT3", "APLNR", "FGF1", "CD34", "HK2", "KDR"
))

for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(Endothelial, genes[i], group.by = "tumor_type", assay = "MAGIC_RNA") + geom_boxplot() + NoLegend()
}

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/Endothelial_Recurrent_high_tumortype_3.pdf")
plot_list
dev.off()

#### Quantifying the cell cycling cell percentage in the object
# Load cell cycle genes
cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score cell cycle phases
Endothelial <- CellCycleScoring(Endothelial, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Get the cell cycle scores
Endothelial$CT_phase <- paste(Endothelial$celltype_marker, Endothelial$Phase, sep = "_")

p <- DimPlot(Endothelial)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/UMAP/cellcycle.pdf")
p
dev.off()

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/"
pdf(paste(savedir, "S_score_prim_rec.pdf"))
VlnPlot(Endothelial, features = "S.Score", group.by = "tumor_type")
dev.off()

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Endothelial/VlnPlot/"
pdf(paste(savedir, "G2M_score_prim_rec.pdf"))
VlnPlot(Endothelial, features = "G2M.Score", group.by = "tumor_type")
dev.off()

cycling_table <- table(Endothelial$CT_phase) %>% as.data.frame()

CT <- gsub("_.*.", "", cycling_table$Var1) %>% unique()

data2 <- cycling_table[, -1]
# Calculate the percentages
data2 <- data2 %>%
    group_by(CT) %>%
    mutate(
        total = sum(Freq),
        percentage = (Freq / total) * 100
    ) %>%
    select(-total)

# Recast the table
recast_data <- data2 %>%
    select(-Freq) %>%
    spread(key = phase, value = percentage)

write.table(recast_data, "/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/cycling_distribution.txt", sep = "\t", quote = F, col.names = T, row.names = T)
