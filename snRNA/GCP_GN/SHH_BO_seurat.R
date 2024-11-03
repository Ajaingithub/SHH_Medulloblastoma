### Running the SHH project. Using the Bo's Seurat object as currently we donot have the fastq files to generate our own

library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
shh_snRNA <- readRDS(paste(savedir, "saveRDS_obj/SHH.merged_all_samples.hg38.qc.filtered.rds", sep = ""))

counts <- GetAssayData(shh_snRNA, assay = "RNA")
counts <- counts[-(grep("^MT-|^RPL|^RPS", rownames(counts))), ]
shh_snRNA_subset <- subset(shh_snRNA, features = rownames(counts))

dir.create(paste(savedir, "dimplot", sep = ""), showWarnings = FALSE)

## QC and creating the Seurat Object made a minute change in scRNA_QC so it will read table
## We will move forward without running this
source("/diazlab/data3/.abhinav/resources/all_scripts/R/scRNA_QC1.R")
for (i in 1:length(obj_name_2)) {
  try({
    real_name <- paste("obj_", obj_name_2[i], sep = "")
    savedir2 <- paste(savedir, "subset/", real_name, "/", sep = "")
    dir.create(savedir2, showWarnings = FALSE, recursive = TRUE)

    obj <- scRNA_QC(
      GEX_obj = get(real_name),
      Sample = real_name,
      saveDir = savedir2,
      min_cells = 3,
      min_genes = 200,
      max_genes = 5000,
      var_genes = 4000,
      mitopercent = 10,
      regress_features = NULL
    )

    real_name_2 <- paste(real_name, "_2", sep = "")
    assign(real_name_2, obj)
  })
}

### Adding the Run information to it
library(stringr)
sample_batch <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_batch.txt", header = TRUE)

shh_snRNA_subset_subset@meta.data$SampleID <- gsub("_.*.", "", shh_snRNA_subset_subset@meta.data$sample_id)
patterns <- sample_batch$SampleID
replacements <- sample_batch$Run
shh_snRNA_subset@meta.data$Run <- shh_snRNA_subset@meta.data$SampleID
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  shh_snRNA_subset@meta.data$Run <- str_replace_all(shh_snRNA_subset@meta.data$Run, pattern, replacements[i])
}

#### Need to make some changes in the samplename in the orig.ident
all_metadata <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/snRNA_metadata.txt", header = TRUE, sep = "\t")
all_metadata$Sample.ID2 <- gsub("7316-", "", all_metadata$Sample.ID)
all_metadata$Tumor.type <- gsub("Primary", "P", all_metadata$Tumor.type) %>% gsub("Recurrent", "R", .)
all_metadata$Batch..snRNA.seq. <- gsub("Batch", "Run", all_metadata$Batch..snRNA.seq.)
all_metadata$orig.ident <- paste(all_metadata$Project, all_metadata$Sample.ID2, "_",
  all_metadata$Age, "_", all_metadata$Sex, "_",
  all_metadata$Tumor.type, "_", all_metadata$Batch..snRNA.seq.,
  sep = ""
)

patterns <- all_metadata$Sample.ID
replacements <- all_metadata$orig.ident
shh_snRNA_subset@meta.data$orig.ident <- shh_snRNA_subset@meta.data$SampleID
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  shh_snRNA_subset@meta.data$orig.ident <- str_replace_all(shh_snRNA_subset@meta.data$orig.ident, pattern, replacements[i])
}

# We will not be integrating the dataset since the tumor is highly heterogenous and integration does not really work on that
# instead we will perform for each sample reference annotation individually
source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
objname <- "SHH"
Assay <- "RNA"
process <- "sctransform"
shh_sctransformed <- sctransform_V2_integration(
  obj = shh_snRNA_subset,
  saveDir = savedir,
  ngenes = 4000,
  regress = c("nCount_RNA"),
  dims = 30,
  Assay = Assay,
  process = process,
  objname = objname,
  split_by = "Run",
  reference = NULL,
  sample_tree = NULL
)

### Integration by batches is not really looking so I am trying to see do we need the integration or not
# DefaultAssay(shh_sctransformed) <- "RNA"
# shh_sctransformed <- NormalizeData(shh_sctransformed, normalization.method = "LogNormalize", scale.factor = 10000)
# shh_sctransformed <- FindVariableFeatures(shh_sctransformed, selection.method = "vst", nfeatures = 4000)
# shh_sctransformed <- ScaleData(shh_sctransformed, features = rownames(shh_sctransformed))
# shh_sctransformed <- RunPCA(shh_sctransformed, features = VariableFeatures(object = shh_sctransformed))
# shh_sctransformed <- FindNeighbors(shh_sctransformed, dims = 1:30)
# shh_sctransformed <- RunUMAP(shh_sctransformed, dims = 1:30)

# pdf(paste(savedir,"UMAP/no_integration.pdf",sep = ""),width = 10, height = 8)
# print(DimPlot(shh_sctransformed, split.by = "orig.ident", group.by="orig.ident", ncol = 5) + NoLegend())
# dev.off()

# pdf(paste(savedir,"UMAP/no_integration_combine.pdf",sep = ""),width = 10, height = 5)
# print(DimPlot(shh_sctransformed, group.by = "orig.ident"))
# dev.off()

### So we will be adopting three strategy
# 1. Unbiased including all the celltypes from the
# 2. Biased including only Glutamatergic celltypes
# 3. Performing FindAllMarkers to identify the markers for each celltypes

# In the SCTransform workflow, we perform a better version of variance stabilization, so we do not scale in this case.
# We fix the slope parameter of the GLM to ln(10) with log10(total UMI) used as the predictor as proposed by Lause et al.
# We utilize an improved parameter estimation procedure that alleviates uncertainty and bias that result from fitting GLM models for very lowly expressed genes.
# We place a lower bound on gene-level standard deviation when calculating Pearson residuals. This prevents genes with extremely low expression (only 1-2 detected UMIs) from having a high pearson residual.
# https://github.com/satijalab/seurat/issues/3003
split_obj <- SplitObject(shh_snRNA_subset, split.by = "sample")
obj_name <- names(split_obj)
obj_name_2 <- gsub("-", "_", obj_name)

for (i in 1:length(obj_name)) {
  obj <- split_obj[[obj_name[i]]]
  assign(paste("obj", obj_name_2[i], sep = "_"), obj)
}

object_name <- grep("_name", ls(pattern = "obj_"), invert = TRUE, value = TRUE)

for (i in 1:length(object_name)) {
  obj <- get(object_name[i])
  # s_genes <- cc.genes$s.genes
  # g2m_genes <- cc.genes$g2m.genes
  # # Calculate percentage of mitochondrial genes
  # obj[["mito_percent"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # # Calculate percentage of ribosomal genes
  # obj[["ribo_percent"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- CellCycleScoring(obj, g2m.features = g2m_genes, s.features = s_genes)

  obj <- SCTransform(obj,
    vst.flavor = "v2", # using version 2
    vars.to.regress = c("nCount_RNA", "nFeature_RNA"),
    assay = "RNA",
    variable.features.n = 2000,
    ncells = nrow(obj@meta.data)
  )

  obj <- RunPCA(object = obj)

  dir.create(paste(savedir, "subset/", object_name[i], "/elbow_plots/", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, "subset/", object_name[i], "/elbow_plots/", object_name[i], "_sctransform.pdf", sep = ""), width = 8, height = 8)
  print(ElbowPlot(obj, ndims = 50))
  dev.off()

  saveDir <- paste(savedir, "subset/", object_name[i], "/", sep = "")
  source("
  ")
  process <- "SCT"
  Assay <- "RNA"
  obj <- RNA_integration(obj, saveDir,
    dims = 20, RNA_features = c("GLI1", "RBFOX3"),
    Assay = Assay, process = process, objname = paste(object_name[i], sep = ""),
    ncol = 1, ndims = 50, doimputation = FALSE
  )

  source("/diazlab/data3/.abhinav/resources/all_scripts/R/cluster_UMAP_QC.R")
  res <- c(0.2, 0.4, 0.6, 0.8, 1)
  res <- 0.6
  for (j in 1:length(res)) {
    process <- paste("RNA_UMAP_QC", res[j], sep = "_")
    Assay <- "SCT"
    obj <- cluster_UMAP_and_QC(
      obj_path = obj, dims = 20, res = res[j],
      saveDir = saveDir, Assay = Assay,
      QC_features = c(
        "nCount_RNA", "nFeature_RNA",
        "percent.mt", "percent.rb",
        "S.Score", "G2M.Score"
      ),
      objname = paste(object_name[i], sep = ""),
      process = process, col_sel = c("orig.ident")
    )
  }
  assign(paste("obj", object_name[i], sep = "_"), obj)
  dir.create(paste(saveDir, "saveRDS_obj/", sep = ""), showWarnings = FALSE)
  saveRDS(obj, paste(saveDir, "saveRDS_obj/", object_name[i], ".RDS", sep = ""))
}


#### Performing reference mapping
### Unbiased
library(Seurat)
reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/reference_all.RDS")
maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
filepath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", recursive = TRUE, pattern = ".RDS", full.names = TRUE)
objname <- gsub(".RDS", "", basename(filepath))

for (i in 1:length(filepath)) {
  try({
    object <- readRDS(filepath[i])
    DefaultAssay(object) <- "SCT"

    anchors <- FindTransferAnchors(
      reference = reference,
      query = object,
      normalization.method = "SCT",
      dims = 1:50
    )

    predictions.assay <- TransferData(
      anchorset = anchors,
      refdata = reference$figure_clusters,
      prediction.assay = TRUE,
      weight.reduction = object[["pca"]],
      dims = 1:30, k.weight = 50
    )

    object[["predictions"]] <- predictions.assay

    celltypes <- unique(reference@meta.data$figure_clusters)

    DefaultAssay(object) <- "predictions"
    pdf(paste(maindir, objname[i], "/featureplot/featureplot_reference_mapped.pdf", sep = ""), width = 10, height = 9)
    FeaturePlot(object, features = c(celltypes), ncol = 3)
    dev.off()

    prediction <- TransferData(
      anchorset = anchors,
      refdata = reference$figure_clusters,
      prediction.assay = FALSE,
      weight.reduction = object[["pca"]],
      dims = 1:30
    )

    object <- AddMetaData(object, metadata = prediction)

    p <- DimPlot(object,
      group.by = "predicted.id", reduction = "umap",
      cols = c(
        "10-Glia" = "#1f77b4", "02-RL" = "#ff7f0e", "03-GCP" = "#279e68", "01-PC" = "#d62728", "15-Meninges" = "#aa40fc",
        "05-eCN/UBC" = "#8c564b", "08-BG" = "#e377c2", "09-Ast" = "#b5bd61", "19-Ast/Ependymal" = "#b5bd61", "11-OPC" = "#17becf",
        "04-GN" = "#aec7e8", "07-PIP" = "#ffbb78", "14-Microglia" = "cadetblue3", "06-iCN" = "cornsilk4", "20-Choroid" = "plum2",
        "18-MLI" = "yellow", "13-Endothelial" = "hotpink", "16-Pericytes" = "black", "12-Committed OPC" = "cyan"
      )
    )

    dir.create(paste(maindir, objname[i], "/UMAP/", sep = ""), showWarnings = FALSE)
    pdf(paste(maindir, objname[i], "/UMAP/dimplot_celltypes_umap.pdf", sep = ""))
    print(p)
    dev.off()


    assign(paste(objname[i], "_CT", sep = ""), object)
    # saveRDS(object, paste(maindir,objname[i],"saveRDS_obj/",objname[i],"_ref_mapped.rds",sep = ""))
  })
}

### Making a table for celltype distribution
library(Seurat)
library(dplyr)
obj_path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset", pattern = "_refglut_mapped_imputed.rds", full.names = TRUE, recursive = TRUE)
objname <- gsub("_refglut_mapped_imputed.rds", "", basename(obj_path))

for (i in 1:length(obj_path)) {
  obj <- readRDS(obj_path[i])
  celltype_table <- table(obj@meta.data$predicted.id) %>% as.data.frame()
  rownames(celltype_table) <- (celltype_table$Var1)
  assign(paste(objname[i], "_glut_table", sep = ""), celltype_table)
}

all_table <- ls(pattern = "ref_table")

# List of data frames
df_list <- list(
  get(all_table[1]), get(all_table[2]), get(all_table[3]), get(all_table[4]), get(all_table[5]),
  get(all_table[6]), get(all_table[7]), get(all_table[8]), get(all_table[9]), get(all_table[10]),
  get(all_table[11]), get(all_table[12]), get(all_table[13]), get(all_table[14]), get(all_table[15]),
  get(all_table[16]), get(all_table[17]), get(all_table[18]), get(all_table[19])
)

# Perform a full outer join
library(purrr)
merged_df <- reduce(df_list, full_join, by = "Var1")

# Print the merged data frame
print(merged_df)
colnames(merged_df) <- c("celltypes", gsub("_table", "", all_table))
rownames(merged_df) <- merged_df[, 1]
merged_df <- merged_df[, -1]
merged_df[is.na(merged_df)] <- 0 # converting NA to 0
total <- colSums(merged_df)
merged_df_2 <- merged_df

## Calculating the percentage
for (i in 1:ncol(merged_df)) {
  merged_df_2[, i] <- round((merged_df[, i] / total[i]) * 100, 0)
}

write.table(merged_df_2, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/ref_mapped_combined_percentage_integer.txt",
  sep = "\t", quote = F, row.names = T, col.names = T
)

## Did this on the laptop as it is easy
setwd("~/Documents/Industry/UCSF/DIaz_lab/projects/SHH/snRNAseq/subset/table/")
data2 <- read.table("ref_mapped_combined_all_percentage.txt", sep = "\t", row.names = 1, header = TRUE)
data2$celltype <- rownames(data2)
data2_melted <- melt(data2)
colnames(data2_melted) <- c("celltype", "samples", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent
c(
  "obj_7316_4529_ref", "obj_7316_5881_ref", "obj_7316_278_ref", "obj_7316_333_ref", "obj_7316_737_ref", "obj_7316_2118_ref",
  "obj_7316_931_ref", "obj_7316_3023_ref", "obj_7316_2978_ref", "obj_7316_311_ref", "obj_7316_1666_ref",
  "obj_7316_1676_ref", "obj_DOD4182_ref", "obj_SF10961_ref", "obj_SF12930_ref", "obj_SF7994_ref",
  "obj_SF8368_ref", "obj_SF8539_ref", "obj_SF9232_ref"
) -> sample_levels

data2_melted$samples <- factor(data2_melted$samples, levels = sample_levels)


p <- ggplot(data2_melted, aes(fill = celltype, y = percentage, x = samples)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "yellow", "hotpink", "black", "blanchedalmond", "blue", "white"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#### Projecting it on the samplen integrated object
## creating a table
SHH_sample <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/saveRDS_obj/shh_integrated_RNA_NN_cluster.RDS")
objname <- grep("GCP|BO|path", ls(pattern = "obj_"), invert = TRUE, value = TRUE)
for (i in 1:length(objname)) {
  obj <- get(objname[i])
  CT <- as.data.frame(obj@meta.data$predicted.id)
  rownames(CT) <- rownames(obj@meta.data)
  colnames(CT) <- "celltype"
  assign(paste(objname[i], "_CT", sep = ""), CT)
}

rbind(
  obj_7316_1666_CT, obj_7316_1676_CT, obj_7316_2118_CT, obj_7316_278_CT, obj_7316_2978_CT,
  obj_7316_3023_CT, obj_7316_311_CT, obj_7316_333_CT, obj_7316_4529_CT, obj_7316_5881_CT,
  obj_7316_737_CT, obj_7316_931_CT, obj_DOD4182_CT, obj_SF10961_CT, obj_SF12930_CT, obj_SF7994_CT,
  obj_SF8368_CT, obj_SF8539_CT, obj_SF9232_CT
) -> all_sample_CT

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
write.table(all_sample_CT, paste(savedir, "table/all_reference_mapped.txt", sep = ""),
  quote = F, col.names = T, row.names = T, sep = "\t"
)

all_sample_CT <- read.table(paste(savedir, "table/all_reference_mapped.txt", sep = ""), header = TRUE, row.names = 1, sep = "\t")
all_sample_CT$cellname <- rownames(all_sample_CT)
all_sample_CT_order <- all_sample_CT[match(rownames(SHH_sample@meta.data), rownames(all_sample_CT)), ]
stopifnot(all(rownames(all_sample_CT_order) == rownames(SHH_sample@meta.data))) ### Sanity check
SHH_sample@meta.data$celltype_samplewise <- all_sample_CT_order$celltype

p <- DimPlot(SHH_sample,
  reduction = "umap", label = FALSE,
  group.by = "celltype_samplewise", raster = TRUE,
  cols = c(
    "10-Glia" = "#1f77b4", "02-RL" = "#ff7f0e", "03-GCP" = "#279e68", "01-PC" = "#aec7e8", "15-Meninges" = "#aa40fc",
    "05-eCN/UBC" = "#8c564b", "08-BG" = "#e377c2", "09-Ast" = "#b5bd61", "19-Ast/Ependymal" = "#b5bd61", "11-OPC" = "#17becf",
    "04-GN" = "#d62728", "07-PIP" = "#ffbb78", "14-Microglia" = "cadetblue3", "06-iCN" = "cornsilk4", "20-Choroid" = "plum2",
    "17-Brainstem" = "black", "18-MLI" = "yellow", "13-Endothelial" = "hotpink", "16-Pericytes" = "deepskyblue",
    "12-Committed OPC" = "cyan", "21-BS Choroid/Ependymal" = "darkslategrey"
  )
)

pdf(paste(savedir, "plots/sample_integrated_individual_samplewise_dimplot.pdf", sep = ""), width = 11, height = 7)
p
dev.off()

p <- DimPlot(SHH_sample,
  reduction = "umap", label = FALSE,
  group.by = "celltype_samplewise", raster = TRUE
)

pdf(paste(savedir, "plots/sample_integrated_individual_samplewise_dimplot_Seurat.pdf", sep = ""), width = 10, height = 7)
p
dev.off()

#### Plotting each celltype separately
library(Seurat)

# Directory to save the plots
rm(plot_list)
plot_list <- list()
Idents(SHH_sample) <- SHH_sample@meta.data$celltype_samplewise
celltypes <- unique(SHH_sample@meta.data$celltype_samplewise)
for (i in 1:length(celltypes)) {
  # Subset the Seurat object for the current cell type
  subset_obj <- subset(SHH_sample, idents = celltypes[i])

  # Create a DimPlot for the subset
  plot_list[[i]] <- DimPlot(subset_obj, reduction = "umap")

  # Save the plot to a PDF
  dir.create(paste0(savedir, "plots/celltype_UMAP/"), showWarnings = FALSE)
  pdf(paste0(savedir, "plots/celltype_UMAP/", celltypes[i], ".pdf"), width = 10, height = 7)
  print(plot_list[[i]])
  dev.off()
}

require(gridExtra)
pdf(paste0(savedir, "plots/celltype_UMAP/combined_CT.pdf"), width = 22, height = 15)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]],
  plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]],
  plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]], plot_list[[21]],
  ncol = 6, nrow = 4
)
dev.off()


#### Reference mapping
### Biased on considering Glutamatergic celltypes
### Running it on r/4.3.3 (Default)
reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/reference_rawcount.RDS")

### We will be considering the celltypes PC, RL, GCP, GN, eCN/UBC
cellname_glut <- rownames(reference@meta.data[grep(c("01-PC|02-RL|03-GCP|04-GN|05-eCN/UBC"), reference@meta.data$figure_clusters), ])
reference_subset <- subset(reference, cells = cellname_glut)

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
DefaultAssay(reference_subset) <- "RNA"
reference_subset <- SCTransform(reference_subset, ncells = 3000, verbose = TRUE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

saveRDS(reference_subset, paste("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/human_cerebellar_glut_subset.RDS"))

#### Mapping the patient object on the reference data
reference_subset <- readRDS("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/human_cerebellar_glut_subset.RDS")
maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
filepath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", recursive = TRUE, pattern = ".RDS", full.names = TRUE)
objname <- gsub(".RDS", "", basename(filepath))

for (i in 1:length(filepath)) {
  try({
    object <- readRDS(filepath[i])
    DefaultAssay(object) <- "SCT"

    anchors <- FindTransferAnchors(
      reference = reference_subset,
      query = object,
      normalization.method = "SCT",
      dims = 1:50
    )

    predictions.assay <- TransferData(
      anchorset = anchors,
      refdata = reference_subset$figure_clusters,
      prediction.assay = TRUE,
      weight.reduction = object[["pca"]],
      dims = 1:30, k.weight = 50
    )

    object[["predictions"]] <- predictions.assay

    celltypes <- unique(reference_subset@meta.data$figure_clusters)

    DefaultAssay(object) <- "predictions"
    pdf(paste(maindir, objname[i], "/featureplot/featureplot_reference_glut_mapped.pdf", sep = ""), width = 10, height = 9)
    FeaturePlot(object, features = c(celltypes), ncol = 3)
    dev.off()

    prediction <- TransferData(
      anchorset = anchors,
      refdata = reference_subset$figure_clusters,
      prediction.assay = FALSE,
      weight.reduction = object[["pca"]],
      dims = 1:30
    )

    object <- AddMetaData(object, metadata = prediction)

    p <- DimPlot(object,
      group.by = "predicted.id", reduction = "umap",
      cols = c(
        "10-Glia" = "#1f77b4", "02-RL" = "#ff7f0e", "03-GCP" = "#279e68", "01-PC" = "#d62728", "15-Meninges" = "#aa40fc",
        "05-eCN/UBC" = "#8c564b", "08-BG" = "#e377c2", "09-Ast" = "#b5bd61", "19-Ast/Ependymal" = "#b5bd61", "11-OPC" = "#17becf",
        "04-GN" = "#aec7e8", "07-PIP" = "#ffbb78", "14-Microglia" = "cadetblue3", "06-iCN" = "cornsilk4", "20-Choroid" = "plum2",
        "18-MLI" = "yellow", "13-Endothelial" = "hotpink", "16-Pericytes" = "black", "12-Committed OPC" = "cyan"
      )
    )

    dir.create(paste(maindir, objname[i], "/UMAP/", sep = ""), showWarnings = FALSE)
    pdf(paste(maindir, objname[i], "/UMAP/dimplot_celltypes_umap_glut.pdf", sep = ""))
    print(p)
    dev.off()

    ### Imputing the data
    DefaultAssay(object) <- "RNA"
    library(reticulate)
    library(Rmagic)
    use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
    py_discover_config("magic") # to check
    object <- magic(object, npca = 20) ## imputing the RNA data as for RNA PCs are 20
    DefaultAssay(object) <- "MAGIC_RNA"
    object <- ScaleData(object, features = rownames(object)) # Sacling data on all the features as we need it during the feature plot scaled data for gene that are not variable

    saveRDS(object, paste(maindir, objname[i], "saveRDS_obj/", objname[i], "_refglut_mapped_imputed.rds", sep = ""))
  })
}

### Making a table for celltype distribution
obj_path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset", pattern = "_refglut_mapped.rds", full.names = TRUE, recursive = TRUE)
objname <- gsub("refglut_mapped.rds", "", basename(obj_path))

for (i in 1:length(obj_path)) {
  obj <- readRDS(obj_path[i])
  celltype_table <- table(obj@meta.data$predicted.id) %>% as.data.frame()
  rownames(celltype_table) <- (celltype_table$Var1)
  assign(paste(objname[i], "_glut_table", sep = ""), celltype_table)
}
all_table <- ls(pattern = "glut_table")

# List of data frames
df_list <- list(
  get(all_table[1]), get(all_table[2]), get(all_table[3]), get(all_table[4]), get(all_table[5]),
  get(all_table[6]), get(all_table[7]), get(all_table[8]), get(all_table[9]), get(all_table[10]),
  get(all_table[11]), get(all_table[12]), get(all_table[13]), get(all_table[14]), get(all_table[15]),
  get(all_table[16]), get(all_table[17]), get(all_table[18]), get(all_table[19])
)

# Perform a full outer join
library(purrr)
merged_df <- reduce(df_list, full_join, by = "Var1")

# Print the merged data frame
print(merged_df)
colnames(merged_df) <- c("celltypes", gsub("glut_table", "", all_table))
rownames(merged_df) <- merged_df[, 1]
merged_df <- merged_df[, -1]
merged_df[is.na(merged_df)] <- 0 # converting NA to 0
total <- colSums(merged_df)
merged_df_2 <- merged_df

## Calculating the percentage
for (i in 1:ncol(merged_df)) {
  merged_df_2[, i] <- round((merged_df[, i] / total[i]) * 100, 0)
}

write.table(merged_df_2, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Glut_mapped_combined_percentage_integer.txt",
  sep = "\t", quote = F, row.names = T, col.names = T
)

write.table(merged_df, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Glut_mapped_combined_counts.txt",
  sep = "\t", quote = F, row.names = T, col.names = T
)


#### Projecting it on the samplen integrated object
## creating a table
SHH_sample <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/saveRDS_obj/shh_integrated_RNA_NN_cluster.RDS")
### Loading the BO objects
glut_objpath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", pattern = "_refglut_mapped_imputed.rds", recursive = TRUE, full.names = TRUE)
glut_objname <- gsub("_BO_mapped.rds", "", basename(Bo_objpath))

for (i in 1:length(glut_objname)) {
  obj <- readRDS(glut_objpath[i])
  CT <- as.data.frame(obj@meta.data$predicted.id)
  rownames(CT) <- rownames(obj@meta.data)
  colnames(CT) <- "celltype"
  assign(paste(glut_objname[i], "_CT", sep = ""), CT)
  assign(paste(glut_objname[i], "_glut", sep = ""), obj)
}

rbind(
  obj_7316_1666_CT, obj_7316_1676_CT, obj_7316_2118_CT, obj_7316_278_CT, obj_7316_2978_CT,
  obj_7316_3023_CT, obj_7316_311_CT, obj_7316_333_CT, obj_7316_4529_CT, obj_7316_5881_CT,
  obj_7316_737_CT, obj_7316_931_CT, obj_DOD4182_CT, obj_SF10961_CT, obj_SF12930_CT, obj_SF7994_CT,
  obj_SF8368_CT, obj_SF8539_CT, obj_SF9232_CT
) -> all_sample_CT

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
write.table(all_sample_CT, paste(savedir, "table/glutamatergic_cells/glut_celltype_mapped.txt", sep = ""),
  quote = F, col.names = T, row.names = T, sep = "\t"
)

all_sample_CT <- read.table(paste(savedir, "table/glutamatergic_cells/glut_celltype_mapped.txt", sep = ""), header = TRUE, row.names = 1, sep = "\t")
all_sample_CT$cellname <- rownames(all_sample_CT)
all_sample_CT_order <- all_sample_CT[match(rownames(SHH_sample@meta.data), rownames(all_sample_CT)), ]
stopifnot(all(rownames(all_sample_CT_order) == rownames(SHH_sample@meta.data))) ### Sanity check
SHH_sample@meta.data$glut_mapped_CT_samplewise <- all_sample_CT_order$celltype

p <- DimPlot(SHH_sample,
  reduction = "umap", label = FALSE,
  group.by = "glut_mapped_CT_samplewise", raster = TRUE,
  cols = c(
    "10-Glia" = "#1f77b4", "02-RL" = "#ff7f0e", "03-GCP" = "#279e68", "01-PC" = "#d62728", "15-Meninges" = "#aa40fc",
    "05-eCN/UBC" = "#8c564b", "08-BG" = "#e377c2", "09-Ast" = "#b5bd61", "19-Ast/Ependymal" = "#b5bd61", "11-OPC" = "#17becf",
    "04-GN" = "#aec7e8", "07-PIP" = "#ffbb78", "14-Microglia" = "cadetblue3", "06-iCN" = "cornsilk4", "20-Choroid" = "plum2",
    "18-MLI" = "yellow", "13-Endothelial" = "hotpink", "16-Pericytes" = "black", "12-Committed OPC" = "cyan"
  )
)

pdf(paste(savedir, "plots/glutamatergic_cells/glut_sample_integrated_individual_samplewise_dimplot.pdf", sep = ""), width = 9, height = 7)
p
dev.off()

#### Plotting each celltype separately
library(Seurat)

# Directory to save the plots
rm(plot_list)
plot_list <- list()
Idents(SHH_sample) <- SHH_sample@meta.data$glut_mapped_CT_samplewise
celltypes <- unique(SHH_sample@meta.data$glut_mapped_CT_samplewise)
for (i in 1:length(celltypes)) {
  # Subset the Seurat object for the current cell type
  subset_obj <- subset(SHH_sample, idents = celltypes[i])

  # Create a DimPlot for the subset
  plot_list[[i]] <- DimPlot(subset_obj, reduction = "umap")

  # Save the plot to a PDF
  dir.create(paste0(savedir, "plots/glutamatergic_cells/celltype_UMAP/"), showWarnings = FALSE)
  pdf(paste0(savedir, "plots/glutamatergic_cells/celltype_UMAP/", celltypes[i], ".pdf"), width = 10, height = 7)
  print(plot_list[[i]])
  dev.off()
}

require(gridExtra)
pdf(paste0(savedir, "plots/glutamatergic_cells/celltype_UMAP/combined_CT.pdf"), width = 14, height = 12.5)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]],
  ncol = 3, nrow = 2
)
dev.off()

### Biased on considering BO Glutamatergic celltypes
### Running it on r/4.3.3 (Default)
BO_reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/BO_human_cerebellar/human.dev_cerebellum.glutamatergic.lineage_cells.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
DefaultAssay(reference_subset) <- "RNA"
BO_reference <- SCTransform(BO_reference, ncells = 3000, verbose = TRUE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

saveRDS(BO_reference, paste("/diazlab/data3/.abhinav/resources/reference/BO_human_cerebellar/human.dev_cerebellum.glutamatergic.lineage_cells_SCT.rds"))

#### Mapping the patient object on the reference data
# maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
# filepath = list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", recursive = TRUE, pattern = ".RDS", full.names = TRUE)
# objname = gsub(".RDS","",basename(filepath))

objname <- grep("_GCP|_path", ls(pattern = "obj_"), invert = TRUE, value = TRUE)


for (i in 1:length(objname)) {
  try({
    object <- get(objname[i])
    DefaultAssay(object) <- "SCT"

    anchors <- FindTransferAnchors(
      reference = BO_reference,
      query = object,
      normalization.method = "SCT",
      dims = 1:50
    )

    predictions.assay <- TransferData(
      anchorset = anchors,
      refdata = BO_reference$celltype,
      prediction.assay = TRUE,
      weight.reduction = object[["pca"]],
      dims = 1:30, k.weight = 50
    )

    object[["BO_predictions"]] <- predictions.assay

    celltypes <- as.vector(unique(BO_reference@meta.data$celltype))

    DefaultAssay(object) <- "BO_predictions"
    pdf(paste(maindir, objname[i], "/featureplot/featureplot_reference_BO_ref_mapped.pdf", sep = ""), width = 10, height = 9)
    FeaturePlot(object, features = celltypes)
    dev.off()

    prediction <- TransferData(
      anchorset = anchors,
      refdata = BO_reference$celltype,
      prediction.assay = FALSE,
      weight.reduction = object[["pca"]],
      dims = 1:30
    )

    object <- AddMetaData(object, metadata = prediction)

    p <- DimPlot(object, group.by = "predicted.id", reduction = "umap")

    dir.create(paste(maindir, objname[i], "/UMAP/", sep = ""), showWarnings = FALSE)
    pdf(paste(maindir, objname[i], "/UMAP/dimplot_celltypes_umap_BO_ref.pdf", sep = ""))
    print(p)
    dev.off()

    ### Imputing the data
    # DefaultAssay(object) <- "RNA"
    # library(reticulate)
    # library(Rmagic)
    # use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
    # py_discover_config("magic") # to check
    # object <- magic(object, npca=20) ## imputing the RNA data as for RNA PCs are 20
    # DefaultAssay(object) <- "MAGIC_RNA"
    # object <- ScaleData(object,features=rownames(object)) # Sacling data on all the features as we need it during the feature plot scaled data for gene that are not variable

    saveRDS(object, paste(maindir, objname[i], "/saveRDS_obj/", objname[i], "_BO_mapped.rds", sep = ""))
  })
}

### Making a table for celltype distribution
obj_path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset", pattern = "_BO_mapped.rds", full.names = TRUE, recursive = TRUE)
objname <- gsub("_BO_mapped.rds", "", basename(obj_path))

for (i in 1:length(obj_path)) {
  obj <- readRDS(obj_path[i])
  celltype_table <- table(obj@meta.data$predicted.id) %>% as.data.frame()
  rownames(celltype_table) <- (celltype_table$Var1)
  assign(paste(objname[i], "_BO_mapped_table", sep = ""), celltype_table)
}
all_table <- ls(pattern = "BO_mapped_table")

# List of data frames
df_list <- list(
  get(all_table[1]), get(all_table[2]), get(all_table[3]), get(all_table[4]), get(all_table[5]),
  get(all_table[6]), get(all_table[7]), get(all_table[8]), get(all_table[9]), get(all_table[10]),
  get(all_table[11]), get(all_table[12]), get(all_table[13]), get(all_table[14]), get(all_table[15]),
  get(all_table[16]), get(all_table[17]), get(all_table[18]), get(all_table[19])
)

# Perform a full outer join
library(purrr)
merged_df <- reduce(df_list, full_join, by = "Var1")

# Print the merged data frame
colnames(merged_df) <- c("celltypes", gsub("_BO_mapped_table|obj_", "", all_table))
rownames(merged_df) <- merged_df[, 1]
merged_df <- merged_df[, -1]
merged_df[is.na(merged_df)] <- 0 # converting NA to 0
total <- colSums(merged_df)
merged_df_2 <- merged_df

## Calculating the percentage
for (i in 1:ncol(merged_df)) {
  merged_df_2[, i] <- round((merged_df[, i] / total[i]) * 100, 2)
}

write.table(merged_df_2, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/BO_mapped_combined_percentage.txt",
  sep = "\t", quote = F, row.names = T, col.names = T
)

write.table(merged_df, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/BO_mapped_combined_counts.txt",
  sep = "\t", quote = F, row.names = T, col.names = T
)

merged_df_2 <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Bo_mapped/BO_mapped_combined_percentage.txt", header = )
merged_df_2$celltype <- rownames(merged_df_2)
merged_df_2_melted <- melt(merged_df_2)
colnames(merged_df_2_melted) <- c("celltype", "samples", "value")
library(ggplot2)
# create a dataset
# Stacked + percent
c(
  "X7316_4529", "X7316_5881", "X7316_278", "X7316_333", "X7316_737", "X7316_2118",
  "X7316_931", "X7316_3023", "X7316_2978", "X7316_311", "X7316_1666",
  "X7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994",
  "SF8368", "SF8539", "SF9232"
) -> sample_levels

merged_df_2_melted$samples <- factor(merged_df_2_melted$samples, levels = sample_levels)

p <- ggplot(merged_df_2_melted, aes(fill = celltype, y = value, x = samples)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "yellow", "hotpink", "black", "blanchedalmond", "blue", "white"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/BO_mapped_bargraph_2.pdf", width = 10, height = 7)
p
dev.off()

p <- DimPlot(BO_reference,
  reduction = "umap", group.by = "celltype", pt.size = 1,
  cols = c(
    "RL-SVZ" = "#1f77b4", "GN-late" = "#ff7f0e", "GCP" = "#279e68", "GN-early" = "#d62728", "RL-VZ" = "#aa40fc",
    "UBC-early" = "#8c564b", "UBC-mid" = "#e377c2", "UBC-late" = "#b5bd61"
  )
)

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/plots/BO_reference_mapped.pdf", width = 6, height = 5)
p
dev.off()

#### Projecting it on the samplen integrated object
## creating a table
SHH_sample <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/saveRDS_obj/shh_integrated_RNA_NN_cluster.RDS")
### Loading the BO objects
Bo_objpath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", pattern = "_BO_mapped.rds", recursive = TRUE, full.names = TRUE)
Bo_objname <- gsub("_BO_mapped.rds", "", basename(Bo_objpath))

for (i in 1:length(Bo_objname)) {
  obj <- readRDS(Bo_objpath[i])
  CT <- as.data.frame(obj@meta.data$predicted.id)
  rownames(CT) <- rownames(obj@meta.data)
  colnames(CT) <- "celltype"
  assign(paste(Bo_objname[i], "_CT", sep = ""), CT)
  assign(paste(Bo_objname[i], "_Bo", sep = ""), obj)
}

rbind(
  obj_7316_1666_CT, obj_7316_1676_CT, obj_7316_2118_CT, obj_7316_278_CT, obj_7316_2978_CT,
  obj_7316_3023_CT, obj_7316_311_CT, obj_7316_333_CT, obj_7316_4529_CT, obj_7316_5881_CT,
  obj_7316_737_CT, obj_7316_931_CT, obj_DOD4182_CT, obj_SF10961_CT, obj_SF12930_CT, obj_SF7994_CT,
  obj_SF8368_CT, obj_SF8539_CT, obj_SF9232_CT
) -> all_sample_CT

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
write.table(all_sample_CT, paste(savedir, "table/Bo_mapped/Bo_celltype_mapped.txt", sep = ""),
  quote = F, col.names = T, row.names = T, sep = "\t"
)

all_sample_CT <- read.table(paste(savedir, "table/Bo_mapped/Bo_celltype_mapped.txt", sep = ""), header = TRUE, row.names = 1, sep = "\t")
all_sample_CT$cellname <- rownames(all_sample_CT)
all_sample_CT_order <- all_sample_CT[match(rownames(SHH_sample@meta.data), rownames(all_sample_CT)), ]
stopifnot(all(rownames(all_sample_CT_order) == rownames(SHH_sample@meta.data))) ### Sanity check
SHH_sample@meta.data$BO_mapped_CT_samplewise <- all_sample_CT_order$celltype

p <- DimPlot(SHH_sample,
  reduction = "umap", label = FALSE,
  group.by = "BO_mapped_CT_samplewise", raster = TRUE,
  cols = c(
    "RL-SVZ" = "#1f77b4", "GN-late" = "#ff7f0e", "GCP" = "#279e68", "GN-early" = "#d62728", "RL-VZ" = "#aa40fc",
    "UBC-early" = "#8c564b", "UBC-mid" = "#e377c2", "UBC-late" = "#b5bd61"
  )
)

pdf(paste(savedir, "plots/Bo_mapped/sample_integrated_individual_samplewise_dimplot.pdf", sep = ""), width = 9, height = 7)
p
dev.off()

#### Plotting each celltype separately
library(Seurat)

# Directory to save the plots
rm(plot_list)
plot_list <- list()
Idents(SHH_sample) <- SHH_sample@meta.data$BO_mapped_CT_samplewise
celltypes <- unique(SHH_sample@meta.data$BO_mapped_CT_samplewise)
for (i in 1:length(celltypes)) {
  # Subset the Seurat object for the current cell type
  subset_obj <- subset(SHH_sample, idents = celltypes[i])

  # Create a DimPlot for the subset
  plot_list[[i]] <- DimPlot(subset_obj, reduction = "umap")

  # Save the plot to a PDF
  dir.create(paste0(savedir, "plots/Bo_mapped/celltype_UMAP/"), showWarnings = FALSE)
  pdf(paste0(savedir, "plots/Bo_mapped/celltype_UMAP/", celltypes[i], ".pdf"), width = 10, height = 7)
  print(plot_list[[i]])
  dev.off()
}

require(gridExtra)
pdf(paste0(savedir, "plots/Bo_mapped/celltype_UMAP/combined_CT.pdf"), width = 15, height = 12.5)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]],
  plot_list[[8]],
  ncol = 3, nrow = 3
)
dev.off()


#### By looking both the unbiased and Glutamatergic cell type referencing. The unbiased will provide us better resolution.
### It would be worthwhile to look at the genes that are highly expressed in GCP and finding the common gene that are high in GCP cluster
library(Seurat)
path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset",
  pattern = "_refall_mapped.rds",
  full.names = TRUE, recursive = TRUE
)
objname <- gsub("_refall_mapped.rds", "", basename(path))
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"

for (i in 1:length(objname)) {
  obj <- readRDS(obj_path[i])
  DefaultAssay(obj) <- "RNA"
  assign(objname[i], obj)
  GCP_marker <- FindMarkers(obj, ident.1 = "03-GCP", group.by = "predicted.id", assay = "RNA", only.pos = TRUE)
  assign(paste(objname[i], "_GCP_markers", sep = ""), GCP_marker)
  write.table(GCP_marker, paste(savedir, objname[i], "/Table/", objname[i], "_ref_GCP_markers.txt", sep = ""),
    sep = "\t", quote = F, row.names = T, col.names = T
  )
}

#### Making an upset plot for the GCP markers
### Making an Upset plot ####
GCP_markers <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset",
  pattern = "_ref_GCP_markers.txt",
  full.names = TRUE, recursive = TRUE
)

markers <- ls(pattern = "_GCP_markers")
markers <- grep("1666|1676|SF9232", markers, invert = TRUE, value = TRUE) ### these samples are less than 20% in GCP

library(data.table)
mar_7316_2118 <- rownames(get(markers[1]))
mar_7316_278 <- rownames(get(markers[2]))
mar_7316_2978 <- rownames(get(markers[3]))
mar_7316_3023 <- rownames(get(markers[4]))
mar_7316_311 <- rownames(get(markers[5]))
mar_7316_333 <- rownames(get(markers[6]))
mar_7316_4529 <- rownames(get(markers[7]))
mar_7316_5881 <- rownames(get(markers[8]))
mar_7316_737 <- rownames(get(markers[9]))
mar_7316_931 <- rownames(get(markers[10]))
mar_DOD4182 <- rownames(get(markers[11]))
mar_SF10961 <- rownames(get(markers[12]))
mar_SF12930 <- rownames(get(markers[13]))
mar_SF7994 <- rownames(get(markers[14]))
mar_SF8368 <- rownames(get(markers[15]))
mar_SF8539 <- rownames(get(markers[16]))

n <- max(
  length(mar_7316_2118), length(mar_7316_278), length(mar_7316_2978), length(mar_7316_3023), length(mar_7316_311),
  length(mar_7316_333), length(mar_7316_4529), length(mar_7316_5881), length(mar_7316_737), length(mar_7316_931),
  length(mar_DOD4182), length(mar_SF10961), length(mar_SF12930), length(mar_SF7994), length(mar_SF8368), length(mar_SF8539)
)

length(mar_7316_2118) <- n
length(mar_7316_278) <- n
length(mar_7316_2978) <- n
length(mar_7316_3023) <- n
length(mar_7316_311) <- n
length(mar_7316_333) <- n
length(mar_7316_4529) <- n
length(mar_7316_5881) <- n
length(mar_7316_737) <- n
length(mar_7316_931) <- n
length(mar_DOD4182) <- n
length(mar_SF10961) <- n
length(mar_SF12930) <- n
length(mar_SF7994) <- n
length(mar_SF8368) <- n
length(mar_SF8539) <- n

df <- as.data.frame(cbind(
  (mar_7316_2118), (mar_7316_278), (mar_7316_2978), (mar_7316_3023), (mar_7316_311),
  (mar_7316_333), (mar_7316_4529), (mar_7316_5881), (mar_7316_737), (mar_7316_931),
  (mar_DOD4182), (mar_SF10961), (mar_SF12930), (mar_SF7994), (mar_SF8368), (mar_SF8539)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_GCP_markers", "", markers)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
pdf(paste(savedir, "plots/GCP_marker_upset_plot.pdf", sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_278" = mar_7316_278, "mar_7316_2978" = mar_7316_2978,
  "mar_7316_3023" = mar_7316_3023, "mar_7316_311" = mar_7316_311,
  "mar_7316_333" = mar_7316_333, "mar_7316_4529" = mar_7316_4529, "mar_7316_5881" = mar_7316_5881,
  "mar_7316_737" = mar_7316_737, "mar_7316_931" = mar_7316_931,
  "mar_DOD4182" = mar_DOD4182, "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930,
  "mar_SF7994" = mar_SF7994, "mar_SF8368" = mar_SF8368, "mar_SF8539" = mar_SF8539
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
df_int$sample_num <- str_count(df_int$int, "mar_")
GCP_common <- df_int[(df_int$sample_num > 13), 1]

write.table(df_int, paste(savedir, "table/GCP_marker_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "table/GCP_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Putting the filter for log2FC < 0.5 and adjusted pvalue 0.05

markers <- grep("1666|1676|SF9232", markers, invert = TRUE, value = TRUE)

for (i in 1:length(markers)) {
  GCP_marker <- get(markers[i])
  GCP_marker_common <- GCP_marker[match(GCP_common, rownames(GCP_marker), nomatch = 0), ]
  GCP_common_filt <- GCP_marker_common[GCP_marker_common$avg_log2FC > 0.5 & GCP_marker_common$p_val_adj < 0.05, ]
  assign(paste(markers[i], "filter", sep = ""), GCP_common_filt)
}

### Performing upset plot on these samples
#### Making an upset plot for the GCP markers
### Making an Upset plot ####
markers <- ls(pattern = "markersfilter")

library(data.table)
mar_7316_2118 <- rownames(get(markers[1]))
mar_7316_278 <- rownames(get(markers[2]))
mar_7316_2978 <- rownames(get(markers[3]))
mar_7316_3023 <- rownames(get(markers[4]))
mar_7316_311 <- rownames(get(markers[5]))
mar_7316_333 <- rownames(get(markers[6]))
mar_7316_4529 <- rownames(get(markers[7]))
mar_7316_5881 <- rownames(get(markers[8]))
mar_7316_737 <- rownames(get(markers[9]))
mar_7316_931 <- rownames(get(markers[10]))
mar_DOD4182 <- rownames(get(markers[11]))
mar_SF10961 <- rownames(get(markers[12]))
mar_SF12930 <- rownames(get(markers[13]))
mar_SF7994 <- rownames(get(markers[14]))
mar_SF8368 <- rownames(get(markers[15]))
mar_SF8539 <- rownames(get(markers[16]))

n <- max(
  length(mar_7316_2118), length(mar_7316_278), length(mar_7316_2978), length(mar_7316_3023), length(mar_7316_311),
  length(mar_7316_333), length(mar_7316_4529), length(mar_7316_5881), length(mar_7316_737), length(mar_7316_931),
  length(mar_DOD4182), length(mar_SF10961), length(mar_SF12930), length(mar_SF7994), length(mar_SF8368), length(mar_SF8539)
)

length(mar_7316_2118) <- n
length(mar_7316_278) <- n
length(mar_7316_2978) <- n
length(mar_7316_3023) <- n
length(mar_7316_311) <- n
length(mar_7316_333) <- n
length(mar_7316_4529) <- n
length(mar_7316_5881) <- n
length(mar_7316_737) <- n
length(mar_7316_931) <- n
length(mar_DOD4182) <- n
length(mar_SF10961) <- n
length(mar_SF12930) <- n
length(mar_SF7994) <- n
length(mar_SF8368) <- n
length(mar_SF8539) <- n

df <- as.data.frame(cbind(
  (mar_7316_2118), (mar_7316_278), (mar_7316_2978), (mar_7316_3023), (mar_7316_311),
  (mar_7316_333), (mar_7316_4529), (mar_7316_5881), (mar_7316_737), (mar_7316_931),
  (mar_DOD4182), (mar_SF10961), (mar_SF12930), (mar_SF7994), (mar_SF8368), (mar_SF8539)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_GCP_markersfilter", "", markers)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
pdf(paste(savedir, "plots/GCP_marker_upset_plot_filter.pdf", sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_278" = mar_7316_278, "mar_7316_2978" = mar_7316_2978,
  "mar_7316_3023" = mar_7316_3023, "mar_7316_311" = mar_7316_311,
  "mar_7316_333" = mar_7316_333, "mar_7316_4529" = mar_7316_4529, "mar_7316_5881" = mar_7316_5881,
  "mar_7316_737" = mar_7316_737, "mar_7316_931" = mar_7316_931,
  "mar_DOD4182" = mar_DOD4182, "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930,
  "mar_SF7994" = mar_SF7994, "mar_SF8368" = mar_SF8368, "mar_SF8539" = mar_SF8539
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
df_int$sample_num <- str_count(df_int$int, "mar_")
GCP_filter_common <- df_int[(df_int$sample_num > 13), 1]

write.table(df_int, paste(savedir, "table/GCP_marker_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "table/GCP_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(GCP_filter_common, paste(savedir, "table/GCP_markers.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

GCP_filter_common <- read.table(paste(savedir, "table/unbiased_mapped/GCP_markers.txt", sep = ""), header = T, sep = "\t")[, 1]
#### VlnPLot to check the gene selected for GCP makes sense
objname <- grep("GCP|_path", ls(pattern = "obj_"), invert = TRUE, value = TRUE)

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GCP_filter_common[1:20], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GCP_first_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GCP_filter_common[21:40], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GCP_second_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GCP_filter_common[41:57], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GCP_third_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

#### GN markers
#### By looking both the unbiased and Glutamatergic cell type referencing. The unbiased will provide us better resolution.
### It would be worthwhile to look at the genes that are highly expressed in GCP and finding the common gene that are high in GCP cluster
library(Seurat)
obj_path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset",
  pattern = "_refall_mapped.rds",
  full.names = TRUE, recursive = TRUE
)
objname <- gsub("_refall_mapped.rds", "", basename(path))
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"

for (i in 1:length(objname)) {
  obj <- readRDS(obj_path[i])
  DefaultAssay(obj) <- "RNA"
  assign(objname[i], obj)
  GN_marker <- FindMarkers(obj, ident.1 = "04-GN", group.by = "predicted.id", assay = "RNA", only.pos = TRUE)
  assign(paste(objname[i], "_GN_markers", sep = ""), GN_marker)
  write.table(GN_marker, paste(savedir, objname[i], "/Table/", objname[i], "_ref_GN_markers.txt", sep = ""),
    sep = "\t", quote = F, row.names = T, col.names = T
  )
}

#### Making an upset plot for the GCP markers
### Making an Upset plot ####
GN_marker <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset",
  pattern = "_ref_GN_markers.txt",
  full.names = TRUE, recursive = TRUE
)

markers <- ls(pattern = "_GN_markers")
markers <- grep("_1666|_1676|_278|_3023|_333|_5881|_931|DOD4182|SF7994|SF8368|SF8539", markers, invert = TRUE, value = TRUE) ### these samples are less than 20% in GCP

library(data.table)
mar_7316_2118 <- rownames(get(markers[1]))
mar_7316_2978 <- rownames(get(markers[2]))
mar_7316_311 <- rownames(get(markers[3]))
mar_7316_4529 <- rownames(get(markers[4]))
mar_7316_737 <- rownames(get(markers[5]))
mar_SF10961 <- rownames(get(markers[6]))
mar_SF12930 <- rownames(get(markers[7]))
mar_SF9232 <- rownames(get(markers[8]))

n <- max(
  length(mar_7316_2118), length(mar_7316_2978), length(mar_7316_311), length(mar_7316_4529),
  length(mar_7316_737), length(mar_SF10961), length(mar_SF12930), length(mar_SF9232)
)

length(mar_7316_2118) <- n
length(mar_7316_2978) <- n
length(mar_7316_311) <- n
length(mar_7316_4529) <- n
length(mar_7316_737) <- n
length(mar_SF10961) <- n
length(mar_SF12930) <- n
length(mar_SF9232) <- n

df <- as.data.frame(cbind(
  (mar_7316_2118), (mar_7316_2978), (mar_7316_311), (mar_7316_4529),
  (mar_7316_737), (mar_SF10961), (mar_SF12930), (mar_SF9232)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_GN_markers", "", markers)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
pdf(paste(savedir, "plots/GN_marker_upset_plot.pdf", sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_2978" = mar_7316_2978, "mar_7316_311" = mar_7316_311,
  "mar_7316_4529" = mar_7316_4529, "mar_7316_737" = mar_7316_737,
  "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930, "mar_SF9232" = mar_SF9232
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
library(stringr)
df_int$sample_num <- str_count(df_int$int, "mar_")
GN_common <- df_int[(df_int$sample_num > 6), 1]

write.table(df_int, paste(savedir, "table/GN_marker_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "table/GN_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Putting the filter for log2FC < 0.5 and adjusted pvalue 0.05
markers <- grep("_1666|_1676|_278|_3023|_333|_5881|_931|DOD4182|SF7994|SF8368|SF8539", markers, invert = TRUE, value = TRUE)

for (i in 1:length(markers)) {
  GN_marker <- get(markers[i])
  GN_marker_common <- GN_marker[match(GN_common, rownames(GN_marker), nomatch = 0), ]
  GN_common_filt <- GN_marker_common[GN_marker_common$avg_log2FC > 0.5 & GN_marker_common$p_val_adj < 0.05, ]
  assign(paste(markers[i], "filter", sep = ""), GN_common_filt)
}

### Performing upset plot on these samples
#### Making an upset plot for the GCP markers
### Making an Upset plot ####
markers <- ls(pattern = "markersfilter")

library(data.table)
mar_7316_2118 <- rownames(get(markers[1]))
mar_7316_2978 <- rownames(get(markers[2]))
mar_7316_311 <- rownames(get(markers[3]))
mar_7316_4529 <- rownames(get(markers[4]))
mar_7316_737 <- rownames(get(markers[5]))
mar_SF10961 <- rownames(get(markers[6]))
mar_SF12930 <- rownames(get(markers[7]))
mar_SF9232 <- rownames(get(markers[8]))

list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_2978" = mar_7316_2978, "mar_7316_311" = mar_7316_311,
  "mar_7316_4529" = mar_7316_4529, "mar_7316_737" = mar_7316_737,
  "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930, "mar_SF9232" = mar_SF9232
)

length(mar_7316_2118) <- n
length(mar_7316_2978) <- n
length(mar_7316_311) <- n
length(mar_7316_4529) <- n
length(mar_7316_737) <- n
length(mar_SF10961) <- n
length(mar_SF12930) <- n
length(mar_SF9232) <- n

df <- as.data.frame(cbind(
  (mar_7316_2118), (mar_7316_2978), (mar_7316_311), (mar_7316_4529),
  (mar_7316_737), (mar_SF10961), (mar_SF12930), (mar_SF9232)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_GN_markersfilter", "", markers)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
pdf(paste(savedir, "plots/GN_marker_upset_plot_filter.pdf", sep = ""), width = 18, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_2978" = mar_7316_2978, "mar_7316_311" = mar_7316_311,
  "mar_7316_4529" = mar_7316_4529, "mar_7316_737" = mar_7316_737,
  "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930, "mar_SF9232" = mar_SF9232
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
df_int$sample_num <- str_count(df_int$int, "mar_")
GN_filter_common <- df_int[(df_int$sample_num > 6), 1]

write.table(df_int, paste(savedir, "table/GN_marker_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "table/GN_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(GN_filter_common, paste(savedir, "table/GN_markers.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)


#### VlnPLot to check the gene selected for GN makes sense
objname <- grep("GN|_path", ls(pattern = "obj_"), invert = TRUE, value = TRUE)

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GN_filter_common[1:20], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GN_first_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GN_filter_common[21:40], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GN_second_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GN_filter_common[41:60], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GN_third_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GN_filter_common[61:80], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GN_fourth_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GN_filter_common[81:100], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GN_fifth_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

for (i in 1:length(objname)) {
  obj <- get(objname[i])
  p <- VlnPlot(obj, features = GN_filter_common[100:124], assay = "RNA", group.by = "predicted.id")
  dir.create(paste(savedir, objname[i], "/vlnplot", sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, objname[i], "/vlnplot/", objname[i], "_GN_sixth_20.pdf", sep = ""), width = 50, height = 30)
  print(p)
  dev.off()
}

#### Now we are using Bo reference mapped
#### GN markers
#### By looking both the unbiased and Glutamatergic cell type referencing. The unbiased will provide us better resolution.
### It would be worthwhile to look at the genes that are highly expressed in GCP and finding the common gene that are high in GCP cluster
library(Seurat)
obj_path <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset",
  pattern = "_BO_mapped.rds",
  full.names = TRUE, recursive = TRUE
)
objname <- gsub("_BO_mapped.rds", "", basename(obj_path))
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"

for (i in 1:length(objname)) {
  try({
    obj <- readRDS(obj_path[i])
    DefaultAssay(obj) <- "RNA"
    assign(objname[i], obj)
    GN_marker <- FindMarkers(obj, ident.1 = c("GN-early"), group.by = "predicted.id", assay = "RNA", only.pos = TRUE)
    assign(paste(objname[i], "_GN_markers", sep = ""), GN_marker)
    write.table(GN_marker, paste(savedir, objname[i], "/Table/", objname[i], "_Bo_GN_early_markers.txt", sep = ""),
      sep = "\t", quote = F, row.names = T, col.names = T
    )
  })
}

#### Making an upset plot for the GCP markers
### Making an Upset plot ####
GN_marker <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset",
  pattern = "_Bo_GN_early_markers.txt",
  full.names = TRUE, recursive = TRUE
)

markers <- ls(pattern = "_GN_markers")
markers <- grep("_1666|_1676|_278|_3023|_333|_5881|_931|DOD4182|SF7994|SF8368|SF8539", markers, invert = TRUE, value = TRUE) ### these samples are less than 20% in GCP
markers <- grep("_GN_markersfilter", markers, invert = TRUE, value = TRUE)

library(data.table)
mar_7316_2118 <- rownames(get(markers[1]))
mar_7316_2978 <- rownames(get(markers[2]))
mar_7316_311 <- rownames(get(markers[3]))
mar_7316_4529 <- rownames(get(markers[4]))
mar_7316_737 <- rownames(get(markers[5]))
mar_SF10961 <- rownames(get(markers[6]))
mar_SF12930 <- rownames(get(markers[7]))
mar_SF9232 <- rownames(get(markers[8]))

n <- max(
  length(mar_7316_2118), length(mar_7316_2978), length(mar_7316_311), length(mar_7316_4529),
  length(mar_7316_737), length(mar_SF10961), length(mar_SF12930), length(mar_SF9232)
)

length(mar_7316_2118) <- n
length(mar_7316_2978) <- n
length(mar_7316_311) <- n
length(mar_7316_4529) <- n
length(mar_7316_737) <- n
length(mar_SF10961) <- n
length(mar_SF12930) <- n
length(mar_SF9232) <- n

df <- as.data.frame(cbind(
  (mar_7316_2118), (mar_7316_2978), (mar_7316_311), (mar_7316_4529),
  (mar_7316_737), (mar_SF10961), (mar_SF12930), (mar_SF9232)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_GN_markers", "", markers)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
pdf(paste(savedir, "plots/Bo_GN_early_marker_upset_plot.pdf", sep = ""), width = 20, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_2978" = mar_7316_2978, "mar_7316_311" = mar_7316_311,
  "mar_7316_4529" = mar_7316_4529, "mar_7316_737" = mar_7316_737,
  "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930, "mar_SF9232" = mar_SF9232
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
library(stringr)
df_int$sample_num <- str_count(df_int$int, "mar_")
GN_common <- df_int[(df_int$sample_num > 6), 1]

write.table(df_int, paste(savedir, "table/Bo_GN_early_marker_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "table/Bo_GN_early_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Putting the filter for log2FC < 0.5 and adjusted pvalue 0.05
markers <- grep("_1666|_1676|_278|_3023|_333|_5881|_931|DOD4182|SF7994|SF8368|SF8539", markers, invert = TRUE, value = TRUE)

for (i in 1:length(markers)) {
  GN_marker <- get(markers[i])
  GN_marker_common <- GN_marker[match(GN_common, rownames(GN_marker), nomatch = 0), ]
  GN_common_filt <- GN_marker_common[GN_marker_common$avg_log2FC > 0.5 & GN_marker_common$p_val_adj < 0.05, ]
  assign(paste(markers[i], "filter", sep = ""), GN_common_filt)
}

### Performing upset plot on these samples
#### Making an upset plot for the GCP markers
### Making an Upset plot ####
markers <- ls(pattern = "markersfilter")

library(data.table)
mar_7316_2118 <- rownames(get(markers[1]))
mar_7316_2978 <- rownames(get(markers[2]))
mar_7316_311 <- rownames(get(markers[3]))
mar_7316_4529 <- rownames(get(markers[4]))
mar_7316_737 <- rownames(get(markers[5]))
mar_SF10961 <- rownames(get(markers[6]))
mar_SF12930 <- rownames(get(markers[7]))
mar_SF9232 <- rownames(get(markers[8]))

list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_2978" = mar_7316_2978, "mar_7316_311" = mar_7316_311,
  "mar_7316_4529" = mar_7316_4529, "mar_7316_737" = mar_7316_737,
  "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930, "mar_SF9232" = mar_SF9232
)

length(mar_7316_2118) <- n
length(mar_7316_2978) <- n
length(mar_7316_311) <- n
length(mar_7316_4529) <- n
length(mar_7316_737) <- n
length(mar_SF10961) <- n
length(mar_SF12930) <- n
length(mar_SF9232) <- n

df <- as.data.frame(cbind(
  (mar_7316_2118), (mar_7316_2978), (mar_7316_311), (mar_7316_4529),
  (mar_7316_737), (mar_SF10961), (mar_SF12930), (mar_SF9232)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_GN_markersfilter", "", markers)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
pdf(paste(savedir, "plots/Bo_GN_early_marker_upset_plot_filter.pdf", sep = ""), width = 20, height = 14)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
  "mar_7316_2118" = mar_7316_2118, "mar_7316_2978" = mar_7316_2978, "mar_7316_311" = mar_7316_311,
  "mar_7316_4529" = mar_7316_4529, "mar_7316_737" = mar_7316_737,
  "mar_SF10961" = mar_SF10961, "mar_SF12930" = mar_SF12930, "mar_SF9232" = mar_SF9232
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
df_int$sample_num <- str_count(df_int$int, "mar_")
GN_filter_common <- df_int[(df_int$sample_num > 5), 1]

write.table(df_int, paste(savedir, "table/Bo_GN_early_marker_gene_intersect.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "table/Bo_GN_early_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(GN_filter_common, paste(savedir, "table/Bo_GN_early_markers.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

#### Finding the marker in the integrated object
GN_marker <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Bo_mapped/Bo_GN_early_markers.txt")[, 1]
CT_unbias <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/unbiased_mapped/all_reference_mapped.txt", header = TRUE, sep = "\t")
CT_unbias$cellname <- rownames(CT_unbias)
CT_unbias_order <- CT_unbias[match(rownames(shh_integrated_RNA_NN_cluster@meta.data), rownames(CT_unbias)), ]
stopifnot(all(rownames(CT_unbias_order) == rownames(shh_integrated_RNA_NN_cluster@meta.data)))
shh_integrated_RNA_NN_cluster@meta.data$unbiase_CT <- CT_unbias_order$celltype

DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GN_marker, group.by = "Bo_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/Bo_GN_markers.pdf", width = 11, height = 18)
p
dev.off()

GN_lit_marker <- c("PCNA", "CCDN1", "CCND2", "CNTN2", "RBFOX3", "NEUROD1", "ZIC1", "ZIC2", "TUBB3", "MEG3", "GRIN2B", "GRIN2C")
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GN_lit_marker, group.by = "unbiase_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GN_lit_markers.pdf", width = 9, height = 7)
p
dev.off()

DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GN_marker, group.by = "Bo_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GN_markers_Bo_mapped.pdf", width = 10, height = 25)
p
dev.off()

GN_lit_marker <- c("PCNA", "CCDN1", "CCND2", "CNTN2", "RBFOX3", "NEUROD1", "ZIC1", "ZIC2", "TUBB3", "MEG3", "GRIN2B", "GRIN2C")
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GN_lit_marker, group.by = "unbiase_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GN_lit_markers.pdf", width = 9, height = 7)
p
dev.off()

GCP_marker <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/unbiased_mapped/GCP_markers.txt")[, 1]
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GCP_marker, group.by = "unbiase_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GCP_markers.pdf", width = 10, height = 15)
p
dev.off()

#### Bo Mapped adding to the seurat object
CT_Bo <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Bo_mapped/Bo_celltype_mapped.txt", header = TRUE, sep = "\t")
CT_Bo$cellname <- rownames(CT_Bo)
CT_Bo_order <- CT_Bo[match(rownames(shh_integrated_RNA_NN_cluster@meta.data), rownames(CT_Bo)), ]
stopifnot(all(rownames(CT_Bo_order) == rownames(shh_integrated_RNA_NN_cluster@meta.data)))
shh_integrated_RNA_NN_cluster@meta.data$Bo_CT <- CT_Bo_order$celltype

GCP_marker <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/unbiased_mapped/GCP_markers.txt")[, 1]
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GCP_marker, group.by = "Bo_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GCP_markers_Bo_mapped.pdf", width = 10, height = 15)
p
dev.off()

GCP_marker <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/unbiased_mapped/GCP_markers_with_good_visium.txt")[, 1]
GCP_marker <- factor(GCP_marker, levels = rev(GCP_marker))
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GCP_marker, group.by = "unbiase_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GCP_markers_with_good_vis_unbias_mapped.pdf", width = 7, height = 10)
p
dev.off()

GCP_marker <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Bo_mapped/Bo_GN_markers_with_good_visium_quality.txt")[, 1]
GCP_marker <- factor(GCP_marker, levels = rev(GCP_marker))
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GCP_marker, group.by = "Bo_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/GN_markers_with_good_vis_Bo_mapped.pdf", width = 7, height = 10)
p
dev.off()

GCP_marker <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/table/Bo_mapped/Bo_GN_early_marker_good_vis_quality.txt")[, 1]
GCP_marker <- factor(GCP_marker, levels = rev(GCP_marker))
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
p <- DotPlot(shh_integrated_RNA_NN_cluster, feature = GCP_marker, group.by = "Bo_CT") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/Bo_GN_early_markers_with_good_vis_Bo_mapped.pdf", width = 7, height = 10)
p
dev.off()

#### InferCNV
# module load jags/4.3.0
library(infercnv)

#### Generating Gene position file
### Used gene_to_position.py
# python gtf_to_position_file.py /diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf /diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/genes_position.txt
## Changing the Ensemblid with the gene symbol
ensembl_pos <- read.table("/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/genes_position.txt", header = FALSE)
ensembl_gene <- read.table("/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/gene_id_gene_name_3.txt", header = FALSE)
ensembl_gene_ordered <- ensembl_gene[match(ensembl_pos$V1, ensembl_gene$V1), ]
stopifnot(all(ensembl_gene_ordered$V1 == ensembl_pos$V1)) ## Sanity check
ensembl_pos$V1 <- ensembl_gene_ordered$V2
ensembl_pos_remove_dup <- ensembl_pos[!duplicated(ensembl_pos$V1), ]
rownames(ensembl_pos_remove_dup) <- ensembl_pos_remove_dup$V1
ensembl_pos_remove_dup <- ensembl_pos_remove_dup[, -1]
write.table(ensembl_pos_remove_dup,
  "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/final_gene_position.txt",
  sep = "\t", quote = F, row.names = T, col.names = F
)


### Generating the input
## SF9232
write.table(obj_SF9232@assays$RNA@counts, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/SF9232_raw_counts.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
obj_SF9232@meta.data$tum_norm <- obj_SF9232@meta.data$predicted.id
obj_SF9232@meta.data[grep("03-GCP|04-GN", obj_SF9232@meta.data$tum_norm, invert = TRUE), "tum_norm"] <- "normal"
obj_SF9232@meta.data[grep("03-GCP|04-GN", obj_SF9232@meta.data$tum_norm, invert = FALSE), "tum_norm"] <- "tumor"
tum_norm <- as.data.frame(obj_SF9232@meta.data$tum_norm)
rownames(tum_norm) <- rownames(obj_SF9232@meta.data)
write.table(tum_norm, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/tum_norm_annotation.txt", sep = "\t", col.names = FALSE, row.names = TRUE, quote = F)
write.ta

### Creating the infercnv object
SF9232_count <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/SF9232_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/final_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/tum_norm_annotation.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", rownames(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = SF9232_count,
  gene_order_file = gene_order,
  annotations_file = annot_file,
  ref_group_names = c("normal")
)

# perform infercnv operations to reveal cnv signal
infercnv_obj <- infercnv::run(infercnv_obj,
  cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/output", # dir is auto-created for storing outputs
  cluster_by_groups = T, # cluster
  denoise = T,
  HMM = T,
  no_prelim_plot = TRUE,
  png_res = 60
)


### Test run was successful. Now generating Input it for each sample and will run separately
for (i in 2:(length(objname) - 1)) {
  obj <- get(objname[i])
  dir.create(paste(savedir, objname[i], "/CNV/", sep = ""), showWarnings = FALSE)
  write.table(obj@assays$RNA@counts, paste(savedir, objname[i], "/CNV/", objname[i], "_raw_counts.txt", sep = ""),
    sep = "\t", quote = F, row.names = T, col.names = T
  )

  obj@meta.data$tum_norm <- obj@meta.data$predicted.id
  obj@meta.data[grep("03-GCP|04-GN", obj@meta.data$tum_norm, invert = TRUE), "tum_norm"] <- "normal"
  obj@meta.data[grep("03-GCP|04-GN", obj@meta.data$tum_norm, invert = FALSE), "tum_norm"] <- "tumor"
  tum_norm <- as.data.frame(obj@meta.data$tum_norm)
  rownames(tum_norm) <- rownames(obj@meta.data)

  write.table(tum_norm, paste(savedir, objname[i], "/CNV/", objname[i], "_tum_norm_annotation.txt", sep = ""),
    sep = "\t", col.names = FALSE, row.names = TRUE, quote = F
  )
}

### Generating script for each sample
# List of sample names
samples <- c(
  "obj_7316_1666", "obj_7316_1676", "obj_7316_2118", "obj_7316_278",
  "obj_7316_2978", "obj_7316_3023", "obj_7316_311", "obj_7316_333",
  "obj_7316_4529", "obj_7316_5881", "obj_7316_737", "obj_7316_931"
)

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
num_cores <- detectCores()
cat("Number of cores detected:", num_cores, "\n")

# Register the parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sfinal_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- paste("X",gsub("-", ".", row.names(annot_file)),sep="")
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/output",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/output",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             no_prelim_plot=TRUE,
                             png_res=300)

saveRDS(infercnv_obj, "%s%s/CNV/output/%s_infercnv.RDS")

# Stop the cluster
stopCluster(cl)
', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)

  # Write the script content to a file
  script_file <- paste0("script_", sample, ".R")
  writeLines(script_content, script_file)
}

### for other samples
samples <- c(
  "obj_DOD4182", "obj_SF10961", "obj_SF12930",
  "obj_SF7994", "obj_SF8368", "obj_SF8539", "obj_SF9232"
)

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
num_cores <- detectCores()
cat("Number of cores detected:", num_cores, "\n")

# Register the parallel backend
cl <- makeCluster(10)
registerDoParallel(cl)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sfinal_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", row.names(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/output",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/output",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             no_prelim_plot=TRUE,
                             png_res=300)

saveRDS(infercnv_obj, "%s%s/CNV/output/%s_infercnv.RDS")

# Stop the cluster
stopCluster(cl)
', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)

  # Write the script content to a file
  script_file <- paste0("script_", sample, ".R")
  writeLines(script_content, script_file)
}

### Running it on Slurm
# at this location /diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV
# cat submit_template.sh
# !/bin/bash
# #SBATCH --job-name=%s
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=4
# #SBATCH --mem=200G
# #SBATCH --time=1-00:00:00
# #SBATCH --output=%s_%j.out
# #SBATCH --error=%s_%j.err
# #SBATCH --mail-user=abhinav.jain@ucsf.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50

# source ~/.bashrc
# module load jags/4.3.0
# R CMD BATCH %s

#### cat submit_job.sh
# #!/bin/bash

# List of sample names
# samples=("obj_7316_1666" "obj_7316_1676" "obj_7316_2118" "obj_7316_278"
#          "obj_7316_2978" "obj_7316_3023" "obj_7316_311" "obj_7316_333"
#          "obj_7316_4529" "obj_7316_5881" "obj_7316_737" "obj_7316_931"
#          "obj_DOD4182" "obj_SF10961" "obj_SF12930" "obj_SF7994"
#          "obj_SF8368" "obj_SF8539","obj_SF9232")

# # Base directory path for R scripts
# r_script_dir="."

# # Loop through each sample and create and submit the SLURM job script
# for sample in "${samples[@]}"; do
#   r_script="${r_script_dir}/script_${sample}.R"

#   # Create the SLURM job script from the template
#   slurm_script="submit_${sample}.sh"
#   cp submit_template.sh $slurm_script
#   sed -i "s/%s/${sample}/g" $slurm_script
#   sed -i "s|%s|${r_script}|g" $slurm_script

#   # Submit the SLURM job
#   sbatch $slurm_script
# done

### First making a final object for each sample
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
glut_filepath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",
  pattern = "_refglut_mapped_imputed.rds",
  recursive = TRUE, full.names = TRUE
)
objname <- gsub("_refglut_mapped_imputed.rds", "", basename(glut_filepath))
BO_filepath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",
  pattern = "_BO_mapped.rds",
  recursive = TRUE, full.names = TRUE
)

refall_filepath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",
  pattern = "_refall_mapped.rds",
  recursive = TRUE, full.names = TRUE
)

infercnv_objpath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/",
  pattern = "run.final.infercnv_obj", recursive = TRUE, full.names = TRUE
)

for (i in 1:length(objname)) {
  library(dplyr)
  library(Seurat)
  library(infercnv)

  glut_obj <- readRDS(glut_filepath[i])

  ### Adding the BO predicted celltypes
  BO_obj <- readRDS(BO_filepath[i])
  stopifnot(all(rownames(glut_obj@meta.data) == rownames(BO_obj@meta.data)))
  glut_obj@meta.data$Bo_mapped_predicted <- BO_obj@meta.data$predicted.id

  ### Adding the reference all predicted celltypes
  refall_obj <- readRDS(refall_filepath[i])
  stopifnot(all(rownames(glut_obj@meta.data) == rownames(refall_obj@meta.data)))
  glut_obj@meta.data$refall_mapped_predicted <- refall_obj@meta.data$predicted.id

  #### Adding the infercnv output to the Seurat object
  infercnv_obj <- readRDS(infercnv_objpath[i])
  infercnv_output <- paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", objname[i], "/CNV/output/", sep = "") # dir is auto-created for storing outputs

  infercnv_obj <- infercnv::run(infercnv_obj,
    cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
    out_dir = infercnv_output,
    cluster_by_groups = T, # cluster
    denoise = T,
    HMM = T,
    no_prelim_plot = TRUE,
    BayesMaxPNormal = 0.2,
    png_res = 300
  )

  glut_obj_2 <- glut_obj ### just to generate the files
  rownames(glut_obj_2@meta.data) <- gsub("-", ".", rownames(glut_obj_2@meta.data))

  if (grepl("7316", objname[i])) {
    rownames(glut_obj_2@meta.data) <- paste("X", rownames(glut_obj_2@meta.data), sep = "")
    ### Addding the inference object to seurat object
    glut_obj_2 <- add_to_seurat(
      seurat_obj = glut_obj_2,
      infercnv_output_path = infercnv_output
    )
  } else {
    glut_obj_2 <- add_to_seurat(
      seurat_obj = glut_obj_2,
      infercnv_output_path = infercnv_output
    )
  }

  ### Since addtoseurat is just generating the NA in all the column we have to manually add it up
  infercnv_meta <- read.delim(paste(infercnv_output, "map_metadata_from_infercnv.txt", sep = ""), header = TRUE, sep = "\t")
  ### We have to make . to - in the top_loss, top_duplis and infercnv_meta
  rownames(infercnv_meta) <- gsub("\\.", "-", rownames(infercnv_meta)) %>% gsub("^X", "", .)
  stopifnot(all(rownames(glut_obj@meta.data) == rownames(infercnv_meta)))

  merged_data <- merge(glut_obj@meta.data, infercnv_meta,
    by = "row.names", all = TRUE
  )
  rownames(merged_data) <- merged_data$Row.names
  merged_data <- merged_data[, -1]
  stopifnot(all(rownames(glut_obj@meta.data) == rownames(merged_data)))
  glut_obj@meta.data <- merged_data
  saveRDS(glut_obj, paste(savedir, objname[i], "/saveRDS_obj/", objname[i], "_final.rds", sep = ""))
}

#### Identifying the cell type
WNT <- c("CTNNB1", "DDX3X", "SMARC4", "TP53", "CSNK2B", "KMT2D", "PIK3CA", "BAI3", "EPHA7", "ARID1A", "ARID2", "SYNCRIP", "ATM")
SHH <- c("PTCH1", "TERT", "TP53", "SUFU", "ELP1", "U1snRNA", "MYCN", "GLI1", "GLI2", "SMO", "DDX3X", "KMT2D", "CREBBP", "TCF4", "PTEN", "KMT2C", "FBXW7", "GSE1", "BCOR", "PRKAR1A", "IDH1", "NF1", "PIK3CA", "Nesting", "SOX2", "PDLIM3", "EYA1", "HHIP", "ATOH1", "SFRP1", "LGALS1")
ped_SHH <- c("COL1A1", "COL3A1", "COL4A1", "LAMA1", "LAMA2", "LAMA4", "LAMB1", "CADM2", "CDH11", "PECAM1", "SMOC2", "SPARC", "FN1", "LUM", "FREM2")
G3 <- c("MYC", "GFI1", "GFI1B", "SMARC4", "OTX2", "KBTBD4", "VTN")
G4 <- c("KDM6A", "MYCN", "CDK6", "PRDM6", "CBFA2T2", "CBFA2T3", "UTX", "OTX2", "EOMES", "KCNA1", "RBM24", "UNC5D", "OAS1", "KHDRBS2")
DMG <- c("H3F3A", "HIST1H3B/C", "MYB", "MYBL1", "MAPK", "CDKN2A/CDKN2B", "MICA", "ATRX", "PPM1D", "TERT", "MCL1", "IDH1", "IDH2")
GABA <- c("Pax3", "Pax2", "Ascl1", "Sox2", "Gad1", "PTF1A", "ASCL1", "LHX5", "LBX1")
GCP <- c("Gli1", "Atoh1", "Mki67", "MYCN", "JAG1", "BMP4", "WNT3", "APC", "PAX6", "MATH1", "PTF1A")
GN <- c("Pcna", "Ccdn1", "Ccnd2", "barhl1", "Cntn2", "Rbfox3", "Grin2b", "Calb2", "Eomes", "NEUROD1", "ZIC2")
PN <- c("Calb1", "Slc1a3", "Ppp1r17", "Car8", "Pcp4", "SHH", "PTF1A")
SC <- c("Sox2", "Gfap", "Olig1", "Olig2", "Blbp", "Fabp7", "Pdgfra", "Nestin")
WM <- c("Mbp", "Fth1", "Plp1", "Mobp", "Cts3")
IGL <- c("Snap25", "Garbra6", "Cbln3", "Diras2", "Gm2694")
ML <- c("Kit", "Gm42418", "Nd4")
neurogenesis <- c("Gfap", "APOE", "ALDOC", "AQP4")
Pre_OPC <- c("PDGFRA", "OLIG1", "EGR1", "ASCL1", "SOX10", "ASCL1", "HES6", "BTG2", "DLL1", "EGFR")
OPC <- c("BCAS1", "MBP", "TF", "MOG", "TFEB")
MSC <- c("GAP43", "DDIT3", "TIMP1", "SPP1")
TAM <- c("CD163", "CSF1R")
PDOX <- c("CD276", "GD2", "IL13R2", "EphA2", "HER2")
Other <- c("MEG3", "MKI67", "SLC17A6")
Astrocytes <- c("FABP7", "SLC1A3", "ALDH1L1", "AQP4")

celltypes <- c(
  "WNT", "SHH", "ped_SHH", "DMG", "GABA", "GCP", "GN", "PN", "SC", "WM", "IGL", "ML",
  "neurogenesis", "Pre_OPC", "OPC", "MSC", "TAM",
  "PDOX", "Other", "Astrocytes"
)

### We will make the heatmap, bubbleplot and featureplot for all the celltypes
library(ArchR)
library(Seurat)
library(reticulate)
library(Rmagic)
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"

DefaultAssay(shh_integrated_RNA_NN_cluster) <- "RNA"
shh_integrated_RNA_NN_cluster <- NormalizeData(shh_integrated_RNA_NN_cluster)
shh_integrated_RNA_NN_cluster <- ScaleData(shh_integrated_RNA_NN_cluster)

library(ArchR)
library(Seurat)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

DefaultAssay(shh_integrated_RNA_NN_cluster) <- "RNA"
shh_integrated_RNA_NN_cluster <- magic(shh_integrated_RNA_NN_cluster, npca = 30) ## imputing the RNA data as for RNA PCs are 20
DefaultAssay(shh_integrated_RNA_NN_cluster) <- "MAGIC_RNA"
shh_integrated_RNA_NN_cluster <- ScaleData(shh_integrated_RNA_NN_cluster)

for (i in 1:length(celltypes)) {
  DefaultAssay(shh_integrated_RNA_NN_cluster) <- "RNA"
  gene_ct <- get(celltypes[i])
  gene_ct <- toupper(gene_ct)
  gene_ct <- grep(paste("^", gene_ct, "$", collapse = "|", sep = ""), rownames(shh_integrated_RNA_NN_cluster), value = TRUE, ignore.case = TRUE)

  ### featureplot
  rm(plot_list)
  plot_list <- list()
  source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")
  for (j in 1:length(gene_ct)) {
    p <- FeaturePlot(shh_integrated_RNA_NN_cluster, gene_ct[j], reduction = "umap", raster = TRUE) + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
    plot_list[[j]] <- p
  }


  dir.create(paste(savedir, "/featureplot/RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
  pdf(paste(savedir, "/featureplot/RNA/", celltypes[i], "_RNA_featureplot.pdf", sep = ""), width = 6.5, height = 5)
  print(plot_list)
  dev.off()

  ### BubblePlot
  p <- DotPlot(shh_integrated_RNA_NN_cluster, gene_ct) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) + coord_flip()

  dir.create(paste(savedir, "/dotplot/RNA", sep = ""), showWarnings = FALSE, recursive = TRUE)
  pdf(paste(savedir, "/dotplot/RNA/", celltypes[i], "_RNA_dotplot.pdf", sep = ""))
  print(p)
  dev.off()

  p <- DoHeatmap(obj, gene_ct)

  dir.create(paste(savedir, "/heatmap/RNA", sep = ""), showWarnings = FALSE, recursive = TRUE)
  pdf(paste(savedir, "/heatmap/RNA/", celltypes[i], "_RNA_heatmap.pdf", sep = ""))
  print(p)
  dev.off()

  # DefaultAssay(obj) <- "MAGIC_SCT"
  # gene_ct <- get(celltypes[i])
  # gene_ct <- grep(paste("^",gene_ct,"$",collapse = "|",sep = ""),rownames(obj),value=TRUE, ignore.case = TRUE)

  # ### featureplot
  # rm(plot_list)
  # plot_list <- list()
  # source("/diazlab/data3/.abhinav/resources/all_scripts/R/express_cell_front1.R")

  # for (j in 1:length(gene_ct)) {
  #   p <- featureplot_front(obj, gene_ct[j],
  #                    reduction = "umap", x="umap_1",
  #                    y="umap_2",size=0.8) +
  #                    scale_color_gradientn(colours = ArchRPalettes$solarExtra)

  #   plot_list[[j]] <- p
  # }

  # dir.create(paste(savedir,"featureplot/MAGIC_SCT/",sep =""),showWarnings = FALSE, recursive = TRUE)
  # pdf(paste(savedir,"featureplot/MAGIC_SCT/",celltypes[i],"_MAGICSCT_featureplot.pdf",sep = ""))
  # print(plot_list)
  # dev.off()

  # ### BubblePlot
  # p <- DotPlot(obj, gene_ct) +scale_color_gradientn(colours = ArchRPalettes$solarExtra) + coord_flip()

  # dir.create(paste(savedir,"dotplot/MAGIC_SCT",sep = ""), showWarnings=FALSE)
  # pdf(paste(savedir,"dotplot/MAGIC_SCT/",celltypes[i],"_MAGICSCT_dotplot.pdf",sep = ""))
  # print(p)
  # dev.off()

  # p <- DoHeatmap(obj, gene_ct)
  # dir.create(paste(savedir,"heatmap/MAGIC_SCT",sep = ""), showWarnings=FALSE)
  # pdf(paste(savedir,"heatmap/MAGIC_SCT/",celltypes[i],"_MAGICSCT_heatmap.pdf",sep = ""))
  # print(p)
  # dev.off()
}

saveRDS(shh_integrated_RNA_NN_cluster, "/diazlab/data3/.abhinav/projects/SHH/snRNA/saveRDS_obj/shh_integrated_RNA_NN_cluster_2.RDS")
g
#### Identifying the non-malignant cells gene markers
### Taking the cellnames generated from Bohyeon
all_samples@meta.data$cellnames <- rownames(all_samples@meta.data)
Bo_non_malignant_CT <- all_samples@meta.data[, c("cellnames", "celltype", "celltype2")]
write.table(Bo_non_malignant_CT, "/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/Bo_non_malignant_CT.txt",
  quote = F, row.names = F, col.names = T, sep = "\t"
)

Bo_non_malignant_CT <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/Bo_non_malignant_CT.txt", header = TRUE, sep = "\t")
Idents(all_samples) <- all_samples@meta.data$celltype

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_markers/"
non_mal_cellmarker <- FindAllMarkers(all_samples, group.by = "celltype")
write.table(non_mal_cellmarker, paste(savedir, "seurat/non_mal_findallmarker.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

for (i in 1:length(celltype)) {
  marker_CT <- non_mal_cellmarker[grep(celltype[i], non_mal_cellmarker$cluster), ]
  marker_CT_pos <- marker_CT[marker_CT$avg_log2FC > 0, ]
  marker_CT_pos_top100 <- head(marker_CT_pos, 100)
  write.table(marker_CT_pos_top100, paste0("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_markers/seurat/", celltype[i], "_top100.txt"),
    quote = F, col.names = T, row.names = T, sep = "\t"
  )
}

DefaultAssay(all_samples) <- "RNA"
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check
all_samples <- magic(all_samples, npca = 20) ## imputing the RNA data as for RNA PCs are 20
DefaultAssay(all_samples) <- "MAGIC_RNA"
all_samples <- ScaleData(all_samples, features = rownames(all_samples)) # Sacling data on all the features as we need it during the feature plot scaled data for gene that are not variable

#### Selected genes based on the good quality of Visium
filepath <- list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_markers/seurat/vis_qual",
  pattern = "_selected_genes.txt", full.names = TRUE
)
celltype <- gsub("_selected_genes.txt", "", basename(filepath))


for (i in 1:length(celltype)) {
  genes <- read.table(filepath[i])[, 1]
  p <- DotPlot(all_samples, feature = genes, group.by = "celltype") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  pdf(paste0("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_markers/seurat/", celltype[i], "_markers.pdf"), width = 7, height = 8)
  print(p)
  dev.off()
}

genes <- c("PGF", "IL16", "LAMB1", "PECAM1", "SPP1")
p <- DotPlot(snRNA, feature = genes, group.by = "celltypes", assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/dotplot/"
pdf(paste0(savedir, "interesting_markers_RNA.pdf"), width = 7, height = 4)
print(p)
dev.off()


library(scran)
out <- pairwiseBinom(all_samples@assays$RNA@counts, all_samples$celltype, direction = "up")
markers <- getTopMarkers(out$statistics, out$pairs, n = 100)
label.markers <- lapply(markers, unlist)
label.markers <- lapply(label.markers, unique)
str(label.markers)

### Saving the files for each celltypes
celltypes <- names(label.markers)


for (i in 1:length(celltypes)) {
  CT_marker <- label.markers[[celltypes[i]]]
  write.table(CT_marker, paste(savedir, "singleR/", celltypes[i], "_marker.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
}


p <- DotPlot(all_samples, feature = label.markers$Astrocytes[1:40], group.by = "celltype") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_markers/singleR/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_markers/singleR/Astrocyte_markers.pdf", width = 11, height = 18)
p
dev.off()


#### Identifying the CNV from each samples
### Performing the celltype annotation for each samples
library(Seurat)
library(ArchR)
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
obj_path <- list.files(savedir, pattern = "_final.rds", recursive = TRUE, full.names = TRUE)
objname <- gsub("_final.rds", "", basename(obj_path))

markers_path <- list.files("/diazlab/data3/.abhinav/resources/gene_list/marker_list/", pattern = ".txt", recursive = FALSE, full.names = TRUE)
marker_name <- gsub(".txt", "", basename(markers_path))

for (i in 1:length(obj_path)) {
  obj <- readRDS(obj_path[i])
  for (j in 1:length(marker_name)) {
    DefaultAssay(obj) <- "MAGIC_RNA"
    genes <- read.table(markers_path[j], header = FALSE)[, 1]
    p <- DotPlot(obj, genes) + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste(savedir, objname[i], "/dotplot/MAGIC_RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
    pdf(paste(savedir, objname[i], "/dotplot/MAGIC_RNA/", objname[i], "_", marker_name[j], "_impute_dotplot.pdf", sep = ""), width = 11, height = 18)
    print(p)
    dev.off()

    ### featureplot
    rm(plot_list)
    plot_list <- list()
    for (k in 1:length(genes)) {
      try({
        p <- FeaturePlot(obj, genes[k], reduction = "umap", pt.size = 2) + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
        plot_list[[k]] <- p
      })
    }

    dir.create(paste(savedir, objname[i], "/featureplot/MAGIC_RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
    pdf(paste(savedir, objname[i], "/featureplot/MAGIC_RNA/", objname[i], "_", marker_name[j], "_impute_featureplot.pdf", sep = ""))
    print(plot_list)
    dev.off()

    DefaultAssay(obj) <- "RNA"
    p <- DotPlot(obj, genes) + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste(savedir, objname[i], "/dotplot/RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
    pdf(paste(savedir, objname[i], "/dotplot/RNA/", objname[i], "_", marker_name[j], "_RNA_dotplot.pdf", sep = ""), width = 11, height = 18)
    print(p)
    dev.off()

    ### featureplot
    rm(plot_list)
    plot_list <- list()
    for (k in 1:length(genes)) {
      try({
        p <- FeaturePlot(obj, genes[k], reduction = "umap", pt.size = 2) + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
        plot_list[[k]] <- p
      })
    }

    dir.create(paste(savedir, objname[i], "/featureplot/RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
    pdf(paste(savedir, objname[i], "/featureplot/RNA/", objname[i], "_", marker_name[j], "_RNA_featureplot.pdf", sep = ""), width = 6.5, height = 5)
    print(plot_list)
    dev.off()
  }
}


### Performing Findallmarkers for each cluster in each samples
library(Seurat)
library(ArchR)
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
obj_path <- list.files(savedir, pattern = "_final.RDS", recursive = TRUE, full.names = TRUE)
objname <- gsub("_final.RDS", "", basename(obj_path))

for (i in 1:length(obj_path)) {
  obj <- readRDS(obj_path[i])
  DefaultAssay(obj) <- "RNA"
  # clus_marker <- FindAllMarkers(obj, only.pos = TRUE)
  # write.table(clus_marker, paste(savedir,objname[i],"/Table/",objname[i],"_all_markers.txt",sep = ""), quote = FALSE,
  #             sep = "\t", row.names = TRUE, col.names = TRUE)
  assign(objname[i], obj)
}


### Performed manually annotation now adding it to the seurat_clusters
CT <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/celltype_annotation/cell_markers/CT_annotation_no_proliferating.txt", sep = "\t", header = TRUE)
j <- 1
i <- 1
objname <- grep("_path", ls(pattern = "obj_"), invert = TRUE, value = TRUE)

for (i in 13:(dim(CT)[2] / 2)) {
  subset_CT <- CT[, c(j, j + 1)]
  subset_CT_req <- subset_CT[!(is.na(subset_CT[, 2]) | subset_CT[, 2] == ""), ] # Removing the empty lines

  ### Matching the first column header with the object to make any changes in object metadata
  obj <- get(objname[i])
  seurat_clus <- levels(Idents(obj))

  annot_clus <- strsplit(paste(subset_CT_req[, 2], collapse = ","), ",")[[1]]
  annot_clus_order <- sort(as.numeric(annot_clus))

  ### Matching that the object name and table object is same as well as the cluster numnber
  stopifnot(objname[i] == colnames(subset_CT_req)[1])
  stopifnot(annot_clus_order == seurat_clus)

  obj@meta.data$celltype_marker <- obj@meta.data$seurat_clusters
  for (k in 1:nrow(subset_CT_req)) {
    clus_pattern <- paste("^", gsub(",", "$|^", paste(subset_CT_req[k, ][2])), "$", sep = "")
    obj@meta.data$celltype_marker <- gsub(clus_pattern, subset_CT_req[k, ][1], obj@meta.data$celltype_marker)
  }

  assign(objname[i], obj)
  j <- (i * 2) + 1
}

#### saving the seurat_object
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
objname <- grep("_path", ls(pattern = "obj_"), invert = TRUE, value = TRUE)
for (i in 1:length(objname)) {
  obj <- get(objname[i])
  saveRDS(obj, paste0(savedir, objname[i], "/saveRDS_obj/", objname[i], "_final.RDS"))
}

#### Annotation for all the samples combined
library(Seurat)
library(ArchR)
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"

markers_path <- list.files("/diazlab/data3/.abhinav/resources/gene_list/marker_list/", pattern = ".txt", recursive = FALSE, full.names = TRUE)
marker_name <- gsub(".txt", "", basename(markers_path))

for (j in 1:length(marker_name)) {
  DefaultAssay(obj) <- "MAGIC_RNA"
  genes <- read.table(markers_path[j], header = FALSE)[, 1]
  p <- DotPlot(obj, genes) + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  dir.create(paste(savedir, "/dotplot/MAGIC_RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
  pdf(paste(savedir, "/dotplot/MAGIC_RNA/", marker_name[j], "_impute_dotplot.pdf", sep = ""), width = 11, height = 18)
  print(p)
  dev.off()

  DefaultAssay(obj) <- "RNA"
  p <- DotPlot(obj, genes) + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  dir.create(paste(savedir, "/dotplot/RNA/", sep = ""), showWarnings = FALSE, recursive = TRUE)
  pdf(paste(savedir, "/dotplot/RNA/", marker_name[j], "_RNA_dotplot.pdf", sep = ""), width = 11, height = 18)
  print(p)
  dev.off()
}

celltype <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_annotation.txt", header = TRUE, sep = "\t")
annot_clus <- strsplit(paste(celltype[, 2], collapse = ","), ",")[[1]]
annot_clus_order <- sort(as.numeric(annot_clus))


### Matching that the object name and table object is same as well as the cluster numnber
seurat_clus <- levels(Idents(obj))
stopifnot(annot_clus_order == seurat_clus)

obj@meta.data$celltype_marker <- obj@meta.data$seurat_clusters
for (k in 1:nrow(celltype)) {
  clus_pattern <- paste("^", gsub(",", "$|^", paste(celltype[k, ][2])), "$", sep = "")
  obj@meta.data$celltype_marker <- gsub(clus_pattern, celltype[k, ][1], obj@meta.data$celltype_marker)
}

saveRDS(obj, "/diazlab/data3/.abhinav/projects/SHH/snRNA/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.6_CT.RDS")

p <- DimPlot(obj,
  group.by = "celltype_marker", reduction = "umap",
  cols = c(
    "astrocyte" = "#1f77b4", "endothelial" = "#ff7f0e", "GCP" = "#279e68", "GN early" = "#d62728", "GN late" = "#aa40fc",
    "Monocytes" = "#8c564b", "Oligodendrocytes" = "#e377c2", "OPC" = "#b5bd61", "Pericytes" = "gray"
  )
)

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/UMAP/celltype_marker.pdf", width = 8, height = 6)
p
dev.off()

### Separately collecting the celltype projecting in the integrated UMAP
for (i in 1:length(objname)) {
  obj <- get(objname[i])
  CT <- as.data.frame(obj@meta.data$celltype_marker)
  rownames(CT) <- rownames(obj@meta.data)
  colnames(CT) <- "celltype"
  assign(paste(objname[i], "_CT", sep = ""), CT)
}

rbind(
  obj_7316_1666_CT, obj_7316_1676_CT, obj_7316_2118_CT, obj_7316_278_CT, obj_7316_2978_CT,
  obj_7316_3023_CT, obj_7316_311_CT, obj_7316_333_CT, obj_7316_4529_CT, obj_7316_5881_CT,
  obj_7316_737_CT, obj_7316_931_CT, obj_DOD4182_CT, obj_SF10961_CT, obj_SF12930_CT, obj_SF7994_CT,
  obj_SF8368_CT, obj_SF8539_CT, obj_SF9232_CT
) -> all_sample_CT

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
write.table(all_sample_CT, paste(savedir, "Table/marker_celltype_mapped.txt", sep = ""),
  quote = F, col.names = T, row.names = T, sep = "\t"
)

all_sample_CT <- read.table(paste(savedir, "Table/marker_celltype_mapped.txt", sep = ""), header = TRUE, row.names = 1, sep = "\t")
all_sample_CT$cellname <- rownames(all_sample_CT)
all_sample_CT_order <- all_sample_CT[match(rownames(obj@meta.data), rownames(all_sample_CT)), ]
stopifnot(all(rownames(all_sample_CT_order) == rownames(obj@meta.data))) ### Sanity check
obj@meta.data$celltype_marker_samplewise <- all_sample_CT_order$celltype

p <- DimPlot(obj,
  reduction = "umap", label = FALSE,
  group.by = "celltype_marker_samplewise", raster = TRUE,
  cols = c(
    "Ast/Ependemyal" = "#1f77b4", "astrocyte" = "#ff7f0e", "Endo_mono" = "#279e68", "Endo_Peri" = "#d62728", "Endo_Purkinje" = "#aa40fc",
    "endothelial" = "#8c564b", "GCP" = "#e377c2", "GN_early" = "blue", "GN_late" = "#b5bd61", "Monocyte_T" = "#17becf",
    "Monocytes" = "#aec7e8", "Multiple_Normal" = "#ffbb78", "MYCN_NRG1_high" = "cadetblue3", "Oligodendrocytes" = "cornsilk4", "OPC" = "plum2",
    "Pericytes" = "yellow", "Purkinje" = "hotpink", "Tcells" = "black"
  )
)

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/UMAP/celltype_marker_samplewise.pdf", width = 9, height = 7)
p
dev.off()

#### InferCNV
### Generating the input
## SF9232
obj_SF9232@meta.data$tum_norm <- obj_SF9232@meta.data$celltype_marker
obj_SF9232@meta.data[grep("GCP|GN_early|GN_late", obj_SF9232@meta.data$tum_norm, invert = TRUE), "tum_norm"] <- "normal"
obj_SF9232@meta.data[grep("GCP|GN_early|GN_late", obj_SF9232@meta.data$tum_norm, invert = FALSE), "tum_norm"] <- "tumor"
tum_norm <- as.data.frame(obj_SF9232@meta.data$tum_norm)
rownames(tum_norm) <- rownames(obj_SF9232@meta.data)
write.table(tum_norm, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/obj_SF9232_tum_norm_annotation_marker.txt", sep = "\t", col.names = FALSE, row.names = TRUE, quote = F)

### Creating the infercnv object
SF9232_count <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/obj_SF9232_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/final_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/obj_SF9232_tum_norm_annotation_marker.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", rownames(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = SF9232_count,
  gene_order_file = gene_order,
  annotations_file = annot_file,
  ref_group_names = c("normal")
)

# perform infercnv operations to reveal cnv signal
infercnv_obj <- infercnv::run(infercnv_obj,
  cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/output", # dir is auto-created for storing outputs
  cluster_by_groups = T, # cluster
  denoise = T,
  HMM = T,
  no_prelim_plot = TRUE,
  png_res = 60
)


### Test run was successful. Now generating Input it for each sample and will run separately
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
for (i in 1:(length(objname))) {
  obj <- get(objname[i])
  dir.create(paste(savedir, objname[i], "/CNV/", sep = ""), showWarnings = FALSE)

  obj@meta.data$tum_norm <- obj@meta.data$celltype_marker
  obj@meta.data[grep("GCP|GN_early|GN_late", obj@meta.data$tum_norm, invert = TRUE), "tum_norm"] <- "normal"
  obj@meta.data[grep("GCP|GN_early|GN_late", obj@meta.data$tum_norm, invert = FALSE), "tum_norm"] <- "tumor"
  tum_norm <- as.data.frame(obj@meta.data$tum_norm)
  rownames(tum_norm) <- rownames(obj@meta.data)

  write.table(tum_norm, paste(savedir, objname[i], "/CNV/", objname[i], "_tum_norm_annotation_marker.txt", sep = ""),
    sep = "\t", col.names = FALSE, row.names = TRUE, quote = F
  )
}

### Generating script for each sample
# List of sample names
setwd("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/")
samples <- c(
  "obj_7316_1666", "obj_7316_1676", "obj_7316_2118", "obj_7316_278",
  "obj_7316_2978", "obj_7316_3023", "obj_7316_311", "obj_7316_333",
  "obj_7316_4529", "obj_7316_5881", "obj_7316_737", "obj_7316_931"
)

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
# num_cores <- detectCores()
# cat("Number of cores detected:", num_cores, "\n")

# Register the parallel backend
# cl <- makeCluster(10)
# registerDoParallel(cl)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sfinal_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation_marker.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- paste("X",gsub("-", ".", row.names(annot_file)),sep="")
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/output_marker",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/output_marker",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             BayesMaxPNormal=0.4,
                             HMM=T,
                             no_prelim_plot=TRUE,
                             png_res=300)

saveRDS(infercnv_obj, "%s%s/CNV/output_marker/%s_infercnv.RDS")

# Stop the cluster
# stopCluster(cl)
', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)

  # Write the script content to a file
  script_file <- paste0("marker_script_", sample, ".R")
  writeLines(script_content, script_file)
}

### for other samples
samples <- c(
  "obj_DOD4182", "obj_SF12930",
  "obj_SF7994", "obj_SF8368", "obj_SF8539"
)

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
num_cores <- detectCores()
cat("Number of cores detected:", num_cores, "\n")

# Register the parallel backend
# cl <- makeCluster(10)
# registerDoParallel(cl)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sfinal_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation_marker.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", row.names(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/output_marker",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/output_marker",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             no_prelim_plot=TRUE,
                             png_res=300)

saveRDS(infercnv_obj, "%s%s/CNV/output_marker/%s_infercnv.RDS")

# Stop the cluster
# stopCluster(cl)
', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)

  # Write the script content to a file
  script_file <- paste0("marker_script_", sample, ".R")
  writeLines(script_content, script_file)
}


### Running it on Slurm
# at this location /diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV
# #!/bin/bash
# #SBATCH --job-name=%s
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=50G
# #SBATCH --time=1-00:00:00
# #SBATCH --output=%s_%j.out
# #SBATCH --error=%s_%j.err
# #SBATCH --mail-user=abhinav.jain@ucsf.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50

# source ~/.bashrc
# conda activate
# R CMD BATCH marker_script_%s.R

#### cat submit_job.sh
# #!/bin/bash

# # List of sample names
# samples=("obj_7316_2978" "obj_7316_3023" "obj_7316_5881"
#         "obj_DOD4182" "obj_SF10961"
#          "obj_SF12930" "obj_SF7994"
#          "obj_SF8368" "obj_SF8539" "obj_SF9232")

# # Base directory path for R scripts
# r_script_dir="."

# # Loop through each sample and create and submit the SLURM job script
# for sample in "${samples[@]}"; do
#   r_script="${r_script_dir}/script_${sample}.R"

#   # Create the SLURM job script from the template
#   slurm_script="submit_${sample}.sh"
#   cp marker_submit_template.sh $slurm_script
#   sed -i "s/%s/${sample}/g" $slurm_script
#   sed -i "s|%s|${r_script}|g" $slurm_script

#   # Submit the SLURM job
#   sbatch $slurm_script
# done


### CNV with Q and P arm
#### InferCNV
### Generating the input
## SF9232
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
obj_path <- list.files(savedir, pattern = "_final.RDS", recursive = TRUE, full.names = TRUE)
objname <- gsub("_final.RDS", "", basename(obj_path))

for (i in 1:length(obj_path)) {
  print(i)
  obj <- readRDS(obj_path[i])
  DefaultAssay(obj) <- "RNA"
  assign(objname[i], obj)
}

obj_SF9232@meta.data$tum_norm <- obj_SF9232@meta.data$celltype_marker
obj_SF9232@meta.data[grep("GCP|GN_early|GN_late", obj_SF9232@meta.data$tum_norm, invert = TRUE), "tum_norm"] <- "normal"
obj_SF9232@meta.data[grep("GCP|GN_early|GN_late", obj_SF9232@meta.data$tum_norm, invert = FALSE), "tum_norm"] <- "tumor"
tum_norm <- as.data.frame(obj_SF9232@meta.data$tum_norm)
rownames(tum_norm) <- rownames(obj_SF9232@meta.data)
write.table(tum_norm, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/obj_SF9232_tum_norm_annotation_marker.txt", sep = "\t", col.names = FALSE, row.names = TRUE, quote = F)

### Creating the infercnv object
SF9232_count <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/obj_SF9232_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("/diazlab/data3/.abhinav/resources/reference/gencode/gene_position_add_cytoband.txt", header = FALSE, row.names = 1)
annot_file <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/obj_SF9232_tum_norm_annotation_marker.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", rownames(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = SF9232_count,
  gene_order_file = gene_order,
  annotations_file = annot_file,
  ref_group_names = c("normal")
)

# perform infercnv operations to reveal cnv signal
infercnv_obj <- infercnv::run(infercnv_obj,
  cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_SF9232/CNV/markerannot_P_Q_arm", # dir is auto-created for storing outputs
  cluster_by_groups = T, # cluster
  denoise = T,
  HMM = T,
  no_prelim_plot = TRUE,
  png_res = 300
)


### Test run was successful. Now generating Input it for each sample and will run separately
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
for (i in 1:(length(objname))) {
  obj <- get(objname[i])
  dir.create(paste(savedir, objname[i], "/CNV/", sep = ""), showWarnings = FALSE)

  obj@meta.data$tum_norm <- obj@meta.data$celltype_marker
  obj@meta.data[grep("GCP|GN_early|GN_late", obj@meta.data$tum_norm, invert = TRUE), "tum_norm"] <- "normal"
  obj@meta.data[grep("GCP|GN_early|GN_late", obj@meta.data$tum_norm, invert = FALSE), "tum_norm"] <- "tumor"
  tum_norm <- as.data.frame(obj@meta.data$tum_norm)
  rownames(tum_norm) <- rownames(obj@meta.data)

  write.table(tum_norm, paste(savedir, objname[i], "/CNV/", objname[i], "_tum_norm_annotation_marker_2.txt", sep = ""),
    sep = "\t", col.names = FALSE, row.names = TRUE, quote = F
  )
}

### Generating script for each sample
# List of sample names
dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/")
setwd("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV/chrq_p_CNV_2/")
samples <- c(
  "obj_7316_1666", "obj_7316_1676", "obj_7316_2118", "obj_7316_278",
  "obj_7316_2978", "obj_7316_3023", "obj_7316_311", "obj_7316_333",
  "obj_7316_4529", "obj_7316_5881", "obj_7316_737", "obj_7316_931"
)

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sgene_position_add_cytoband.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation_marker_2.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- paste("X",gsub("-", ".", row.names(annot_file)),sep="")
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/markerannot_P_Q_arm",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/markerannot_P_Q_arm_2",  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             analysis_mode = "subclusters",
                             denoise=T,
                             tumor_subcluster_partition_method = "leiden",
                             HMM=T,
                             HMM_type = "i6",
                             plot_steps = T,
                             write_expr_matrix = T,
                             write_phylo = T,
                             no_prelim_plot=TRUE,
                             output_format = "pdf",
                             num_threads = 32)

saveRDS(infercnv_obj, "%s%s/CNV/markerannot_P_Q_arm/%s_infercnv.RDS")

', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)

  # Write the script content to a file
  script_file <- paste0("marker_script_", sample, ".R")
  writeLines(script_content, script_file)
}

### for other samples
samples <- c("obj_DOD4182", "obj_SF7994", "obj_SF8368", "obj_SF8539")

# Base directory paths
base_data_dir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
base_resource_dir <- "/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/"

# Loop through each sample and create the scripts
for (sample in samples) {
  script_content <- sprintf('
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
num_cores <- detectCores()
cat("Number of cores detected:", num_cores, "\n")

# Register the parallel backend
# cl <- makeCluster(10)
# registerDoParallel(cl)

count <- read.table("%s%s/CNV/%s_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("%sgene_position_add_cytoband.txt", header = FALSE, row.names = 1)
annot_file <- read.table("%s%s/CNV/%s_tum_norm_annotation_marker_2.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- gsub("-", ".", row.names(annot_file))
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("%s%s/CNV/markerannot_P_Q_arm",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="%s%s/CNV/markerannot_P_Q_arm_2",  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             analysis_mode = "subclusters",
                             denoise=T,
                             tumor_subcluster_partition_method = "leiden",
                             HMM=T,
                             HMM_type = "i6",
                             plot_steps = T,
                             write_expr_matrix = T,
                             write_phylo = T,
                             no_prelim_plot=TRUE,
                             output_format = "pdf",
                             num_threads = 32)

saveRDS(infercnv_obj, "%s%s/CNV/markerannot_P_Q_arm/%s_infercnv.RDS")

# Stop the cluster
# stopCluster(cl)
', base_data_dir, sample, sample, base_resource_dir, base_data_dir, sample, sample, base_data_dir, sample, base_data_dir, sample, base_data_dir, sample, sample)

  # Write the script content to a file
  script_file <- paste0("marker_script_", sample, ".R")
  writeLines(script_content, script_file)
}

### Running it on Slurm
# at this location /diazlab/data3/.abhinav/projects/SHH/snRNA/subset/CNV
# marker_submit_template.sh
# #!/bin/bash
# #SBATCH --job-name=%s
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=50G
# #SBATCH --time=1-00:00:00
# #SBATCH --output=%s_%j.out
# #SBATCH --error=%s_%j.err
# #SBATCH --mail-user=abhinav.jain@ucsf.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50

# source ~/.bashrc
# conda activate
# R CMD BATCH marker_script_%s.R

#### cat submit_job.sh
# !/bin/bash

# List of sample names
# samples=("obj_7316_2978" "obj_7316_3023" "obj_7316_5881" "obj_DOD4182" "obj_SF7994" "obj_SF8368" "obj_SF8539" "obj_7316_1666" "obj_7316_1676" "obj_7316_2118" "obj_7316_278" "obj_7316_2978" "obj_7316_3023" "obj_7316_311" "obj_7316_333" "obj_7316_4529" "obj_7316_5881" "obj_7316_737" "obj_7316_931")

# # Base directory path for R scripts
# r_script_dir="."

# # Loop through each sample and create and submit the SLURM job script
# for sample in "${samples[@]}"; do
#   r_script="${r_script_dir}/script_${sample}.R"

#   # Create the SLURM job script from the template
#   slurm_script="submit_${sample}.sh"
#   cp marker_submit_template.sh $slurm_script
#   sed -i "s/%s/${sample}/g" $slurm_script
#   sed -i "s|%s|${r_script}|g" $slurm_script

#   # Submit the SLURM job
#   sbatch $slurm_script
# done

#### Longitudinal analysis
celltypes_ind <- table(objsubset@meta.data$sample, objsubset@meta.data$celltype_marker2)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_individual2_noSF7994.txt"),
  row.names = T, col.names = T, sep = "\t", quote = F
)

celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_individual.txt"), header = TRUE, sep = "\t")
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
c(
  "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
  "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
  "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994",
  "SF8368", "SF8539", "SF9232"
) -> sample_levels

df_melted$sample <- gsub("-", "_", df_melted$sample)
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

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/marker_celltype2.pdf", width = 10, height = 7)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE)

patterns <- long_sample$Sample
patterns <- gsub("SF7994R", "SF7994", patterns)
replacements <- long_sample$Longitudinal
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$sample), ]

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = TRUE, alternative = "greater") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_all_samples_bxplot_one_sided_greater.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()

### barplot
gbr <- ggbarplot(df_melted,
  x = "Celltype", y = "Cell_percentage",
  add = c("mean_se", "jitter"),
  color = "group", palette = c("#00AFBB", "#E7B800"),
  position = position_dodge(0.8)
)

stat.test <- stat.test %>%
  add_xy_position(x = "Celltype", dodge = 0.8)
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

pdf(paste(savedir, "Table/HC_P_barplot_t_test.pdf", sep = ""), width = 8, height = 6)
gbr2
dev.off()

#### Removed sample Recurrent and Primary
objsubset@meta.data$celltypes <- gsub("Neuronal", "GN_intermediate", objsubset@meta.data$celltypes)
celltypes_ind <- table(objsubset@meta.data$sample, objsubset@meta.data$celltypes)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltypes_individual2_noSF7994.txt"),
  row.names = T, col.names = T, sep = "\t", quote = F
)

# celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_individual.txt"), header = TRUE, sep = "\t")
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
c(
  "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
  "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
  "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994",
  "SF8368", "SF8539", "SF9232"
) -> sample_levels

df_melted$sample <- gsub("-", "_", df_melted$sample)
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

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/marker_celltype2_noSF7994.pdf", width = 10, height = 7)
p
dev.off()


#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns) %>% gsub("-", "_", .)
replacements <- long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

# df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",df_melted$sample),]
# df_melted_paired <- df_melted
stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir, "Table/t_test_all_primary_recurrent_paired_samples_bxplot_two_sided_greater.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()

#### Based on the seurat_cluster
# objsubset@meta.data$celltypes <- gsub("Neuronal","GN_intermediate",objsubset@meta.data$celltypes)
celltypes_ind <- table(objsubset@meta.data$sample, objsubset@meta.data$seurat_clusters2)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/seurat_clusters2_individual2_noSF7994.txt"),
  row.names = T, col.names = T, sep = "\t", quote = F
)

# celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_individual.txt"), header = TRUE, sep = "\t")
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
c(
  "7316_4529", "7316_5881", "7316_278", "7316_333", "7316_737", "7316_2118",
  "7316_931", "7316_3023", "7316_2978", "7316_311", "7316_1666",
  "7316_1676", "DOD4182", "SF10961", "SF12930", "SF7994",
  "SF8368", "SF8539", "SF9232"
) -> sample_levels

df_melted$sample <- gsub("-", "_", df_melted$sample)
df_melted$sample <- factor(df_melted$sample, levels = sample_levels)

p <- ggplot(df_melted, aes(fill = celltype, y = percentage, x = sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white", "yellow", "green", "red", "grey23", "grey34"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/seurat_cluster_celltype2_noSF7994.pdf", width = 10, height = 7)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

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
# df_melted_paired <- df_melted

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
stat.test$p_one_sided <- (stat.test$p) / 2
stat_test_df <- as.data.frame(stat.test)

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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"
pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_all_samples_bxplot_two_sided_cluster.pdf", sep = ""), width = 12, height = 6)
bxp4
dev.off()

marker_genes <- c(
  "GFAP", "AQP4", "ALDH1L1", "SLC6A11", "SLC1A2",
  "VWF", "CD34", "CD93", "CDH5", "TIE",
  "GLI2", "SMO", "BCOR", "ATOH1", "PTCH2",
  "NHLH2", "NHLH1", "PRDM8", "STMN2", "TUBB3",
  "NEUROD1", "NEUROD2", "OXR1", "NRXN2", "LMTK3",
  "CD14", "CD163", "CD86", "CSF1R", "CSF3R",
  "MOBP", "PLP1", "MOG", "MAG", "MBP",
  "PDGFRA", "PTPRZ1", "MYT1", "CSPG4", "VCAN",
  "PDGFRB", "ACTA2", "RGS5", "CDH6", "RBPMS"
)

p <- DotPlot(obj, feature = marker_genes, group.by = "celltype_marker") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/markers_genes.pdf", width = 10, height = 15)
p
dev.off()

c(
  "NPAS3", "TNC", "ID4", "FMN2", "CTNND2",
  "ERG", "KDR", "CLDN5", "ENG", "NOTCH4",
  "ITGAX", "MSR1", "CD74", "CLEC7A", "LY86", "CD4", "CD8A", "FOXP3", "GZMK", "CD25", "CD69", "PDCD1", "KLRG1", "NKG7", "PTPRC", "IL2RA",
  "CNP", "CLDN11", "ENPP2", "SOX10", "APLP1",
  "NTNG1", "CNTN1", "OLIG2", "OLIG1", "CD82",
  "NOTCH3", "CD248", "KCNJ8", "ABCC9"
) -> genes

pericyte_marker <- c("PDGFRB", "ACTA2", "RGS5", "CDH6", "RBPMS", "NOTCH3", "CD248", "KCNJ8", "ABCC9")
p <- DotPlot(obj, feature = genes, group.by = "celltype_marker") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/non_maginant_markers_genes.pdf", width = 10, height = 15)
p
dev.off()

p <- DotPlot(objsubset, feature = unique(c(G3, G4)), group.by = "celltypes") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/G3_G4_markers_genes.pdf", width = 10, height = 15)
p
dev.off()

p <- DotPlot(objsubset, feature = unique(c(G3, G4)), group.by = "seurat_clusters") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/", showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/dotplot/MAGIC_RNA/G3_G4_markers_genes_cluster.pdf", width = 10, height = 15)
p
dev.off()

#### Performning the differential gene expression for the paired primary and recurrent
### Bulk Pseudobulk ####
celltype <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_annotation.txt", header = TRUE, sep = "\t")
annot_clus <- strsplit(paste(celltype[, 2], collapse = ","), ",")[[1]]
annot_clus_order <- sort(as.numeric(annot_clus))

### Matching that the object name and table object is same as well as the cluster numnber
seurat_clus <- levels(Idents(obj))
stopifnot(annot_clus_order == seurat_clus)

obj@meta.data$celltype_marker <- obj@meta.data$seurat_clusters
for (k in 1:nrow(celltype)) {
  clus_pattern <- paste("^", gsub(",", "$|^", paste(celltype[k, ][2])), "$", sep = "")
  obj@meta.data$celltype_marker <- gsub(clus_pattern, celltype[k, ][1], obj@meta.data$celltype_marker)
}

#### Splitting the clusters
clusters <- levels(obj@meta.data$seurat_clusters)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
  cellnames <- obj@meta.data[grep(paste("^", clusters[i], "$", sep = ""), obj@meta.data$seurat_clusters), ] %>% rownames()
  p <- DimPlot(obj,
    cells.highlight = cellnames,
    reduction = "umap",
    label = FALSE, cols.highlight = "deeppink2",
    sizes.highlight = 0.5,
    cols = "gray92"
  ) +
    ggtitle(paste(clusters[i])) +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.2))
  plot_list[[i]] <- p
}

pdf(paste(savedir, "UMAP/cluster_splitted.pdf", sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

### Combining the GN early and GN late to GN
obj@meta.data$celltype_marker2 <- gsub("GNearly|GNlate", "GN", obj@meta.data$celltype_marker)

GCP <- c(0, 4, 5, 7, 8, 13)
# GNearly <- c(2,3,6,9,11)
# GNlate <-	c(1,10,17)
GN <- c(2, 3, 6, 9, 11, 1, 10, 17)
astrocyte <- c(16)
endothelial <- c(19)
Monocytes <- c(12)
OPC <- c(18)
Oligodendrocytes <- c(20)
Pericytes <- c(14, 15)

# celltypes = unique(obj@meta.data$celltype_marker2)
celltypes <- "GN"
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
obj@meta.data$tumor_type <- gsub("_Run.*.", "", obj@meta.data$orig.ident) %>% gsub(".*._", "", .)

snRNA_combine <- obj

for (i in 1:length(celltypes)) {
  try({
    cluster <- get(celltypes[i])
    # cluster = "all"
    group1 <- "P"
    group2 <- "R"
    remove_samples <- c("DOD4182", "SF7994R", "SF8368", "SF8539", "SF10961", "SF9232", "SF12930", "CBTN2978", "CBTN311", "CBTN1666", "CBTN1676")
    remove_samples_name <- paste(remove_samples, collapse = "_and_")
    DefaultAssay(obj) <- "RNA"
    # savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
    source("/diazlab/data3/.abhinav/resources/all_scripts/R/pseudobulk_PCA_within_no_batch_correction_AJ.R")
    SHH_dds <- pseudobulk_within_cluster_AJ_no_batch_correction(
      obj = obj, savedir = savedir, group1 = group1, group2 = group2,
      grouping_by = "tumor_type", cluster = cluster, cell_freq = 20, remove_samples = remove_samples,
      cluster_group = "seurat_clusters", sample_col = "orig.ident", batch_col = "Run",
      gene_min_counts = 5,
      column_name_split = c("sampleID", "samplename", "Age", "Sex", "tumor_type", "Run"),
      splitting = "-"
    )


    savedir2 <- paste(savedir, "pseudobulk/clus_", paste(cluster, collapse = "_"), "_removed_", remove_samples_name, "_", group1, "_vs_", group2, "/", sep = "")
    # savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
    source("/diazlab/data3/.abhinav/resources/all_scripts/R/RNAseq_limma_EdgeR.R")
    design0 <- model.matrix(~ 0 + tumor_type, data = colData(SHH_dds))
    colnames(design0) <- c(group1, group2)
    cm <- makeContrasts(P_VS_R = P - R, levels = design0)
    desl_clus <- LimmaEdgeR_differential(
      dds = SHH_dds,
      design0 = design0,
      cm = cm,
      savedir = savedir2,
      logfc = 0.25,
      p_value_adj = 0.05,
      column_name_split = c("sampleID", "samplename", "Age", "Sex", "tumor_type", "Run")
    )
  })
}

### Performing the different of the recurrent vs primary using DESeq2
GCP <- c(0, 4, 5, 7, 8, 13)
# GNearly <- c(2,3,6,9,11)
# GNlate <-	c(1,10,17)
GN <- c(2, 3, 6, 9, 11, 1, 10, 17)
astrocyte <- c(16)
endothelial <- c(19)
Monocytes <- c(12)
OPC <- c(18)
Oligodendrocytes <- c(20)
Pericytes <- c(14, 15)

celltypes <- c("GCP", "GNearly", "GNlate", "astrocyte", "endothelial", "Monocytes", "OPC", "Oligodendrocytes", "Pericytes")
celltypes <- c("GN")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
obj@meta.data$tumor_type <- gsub("_Run.*.", "", obj@meta.data$orig.ident) %>% gsub(".*._", "", .)

snRNA_combine <- obj

for (i in 1:length(celltypes)) {
  try({
    cluster <- get(celltypes[i])
    # cluster = "all"
    group1 <- "P"
    group2 <- "R"
    remove_samples <- c("DOD4182", "SF7994R", "SF8368", "SF8539", "SF10961", "SF9232", "SF12930", "CBTN2978", "CBTN311", "CBTN1666", "CBTN1676")
    remove_samples_name <- paste(remove_samples, collapse = "_and_")
    DefaultAssay(obj) <- "RNA"
    # savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
    source("/diazlab/data3/.abhinav/resources/all_scripts/R/Differential_paired.R")
    savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/"
    result <- Pseudobulk_Differential_paired(obj, savedir, group1, group2,
      grouping_by = "tumor_type", cluster = cluster,
      cluster_group = "seurat_clusters", cell_freq = 10,
      remove_samples = remove_samples, gene_min_counts = 5, sample_col = "orig.ident",
      batch_col = "Run", column_name_split = c("sampleID", "samplename", "Age", "Sex", "tumor_type", "Run"), splitting = "-"
    )
  })
}

diff <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/table/GEX_diff.txt")
celltypes <- c("GCP", "GCP", "GNlate", "GNlate", "Monocytes", "Monocytes", "Pericytes", "Pericytes", "Astrocytes", "Astrocytes", "GNearly", "GNearly", "all", "all")
Regulated <- c("down", "up", "down", "up", "down", "up", "down", "up", "down", "up", "down", "up", "up", "down")
diff <- diff[1:14, ]
colnames(diff) <- c("DEG", "filepath")
diff[, "celltypes"] <- celltypes
diff[, "Regulated"] <- Regulated

diff$celltypes <- factor(diff$celltypes, levels = c("all", "GCP", "GNearly", "GNlate", "Monocytes", "Astrocytes", "Pericytes"))

library(ggpubr)
gbr <- ggbarplot(diff,
  x = "celltypes", y = "DEG",
  add = c("mean_se", "jitter"),
  color = "Regulated", fill = "Regulated", palette = c("#00AFBB", "#E7B800"),
  position = position_dodge(0.8)
)

pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/plots/GEX_number.pdf")
gbr
dev.off()

#### Performing GSEA enrichment for differentiated and undifferentiated
### GSEA for R vs P ####
library(fgsea)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(dplyr)

### For the cluster 7 genes we have to take the Cluster 5, 2, 3 and 2+3 for innate and adaptive
### combining cluster 2 and 3
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/"
dir.create(savedir, showWarnings = FALSE)
maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/"
files <- grep("P_vs_R", list.files(maindir, pattern = "DESeq_result.txt", full.names = TRUE, recursive = TRUE), invert = FALSE, value = TRUE)
Ag_name <- gsub(maindir, "", files) %>%
  gsub("_P_vs_R/DESeq_result.txt", "", .) %>%
  gsub("/", "", .)
Ag_name <- Ag_name[1:7]

# j=6
# add j= x number corresponding to number in Ag_name when only running a subset of comparisons, in that case skip the next line "for (....)" and keep running with "Ag_diff..."


# Pathways to be considered
# Load necessary libraries
library(fgsea)
library(msigdbr)
library(dplyr)
library(clusterProfiler)

# Download C5 Biological Process (BP) gene sets
msigdb_c5_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# Define the specific pathways of interest
pathways_of_interest <- c(
  "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION",
  "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION",
  "GOBP_CEREBRAL_CORTEX_GABAERGIC_INTERNEURON_DIFFERENTIATION",
  "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION",
  "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION",
  "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION",
  "GOBP_GABAERGIC_NEURON_DIFFERENTIATION",
  "GOBP_GLUTAMATERGIC_NEURON_DIFFERENTIATION",
  "GOBP_HYPOTHALAMUS_GONADOTROPHIN_RELEASING_HORMONE_NEURON_DIFFERENTIATION",
  "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION",
  "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION",
  "GOBP_NORADRENERGIC_NEURON_DIFFERENTIATION",
  "GOBP_OLFACTORY_BULB_INTERNEURON_DIFFERENTIATION",
  "GOBP_PERIPHERAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_DOPAMINERGIC_NEURON_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION",
  "GOBP_PYRAMIDAL_NEURON_DIFFERENTIATION",
  "GOBP_REGULATION_OF_DOPAMINERGIC_NEURON_DIFFERENTIATION",
  "GOBP_REGULATION_OF_FOREBRAIN_NEURON_DIFFERENTIATION",
  "GOBP_REGULATION_OF_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION",
  "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION",
  "GOBP_RETINAL_BIPOLAR_NEURON_DIFFERENTIATION",
  "GOBP_SOMATIC_MOTOR_NEURON_DIFFERENTIATION",
  "GOBP_SPINAL_CORD_ASSOCIATION_NEURON_DIFFERENTIATION",
  "GOBP_SPINAL_CORD_MOTOR_NEURON_DIFFERENTIATION",
  "GOBP_VENTRAL_SPINAL_CORD_INTERNEURON_DIFFERENTIATION",
  "GOBP_WNT_SIGNALING_PATHWAY_INVOLVED_IN_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"
)

# Filter the gene sets for pathways of interest
filtered_gene_sets <- msigdb_c5_bp %>% filter(gs_name %in% pathways_of_interest)

# Convert the filtered gene sets to list format
gene_sets <- filtered_gene_sets %>% split(x = .$entrez_gene, f = .$gs_name)

# Run fgsea
fgsea_results <- fgsea(
  pathways = gene_sets, stats = geneList,
  minSize = 8,
  eps = 0.0,
  maxSize = 1500
)

# View results
fgsea_results <- fgsea_results %>% arrange(padj)
print(fgsea_results)

for (j in 1:length(Ag_name)) {
  Ag_diff <- read.table(files[j], header = TRUE, sep = "\t")

  ## extracting the genes required
  Ag_diff_req <- Ag_diff
  Ag_diff_req_Entrez <- bitr(rownames(Ag_diff_req), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  Ag_diff_req_Entrez[duplicated(Ag_diff_req_Entrez$SYMBOL) # Identify the duplicates
  | duplicated(Ag_diff_req_Entrez$SYMBOL, fromLast = TRUE), ]
  # Ag_diff_req_Entrez_unique <- tops_all_Entrez[Ag_diff_req_Entrez$ENTREZID != 7006 & Ag_diff_req_Entrez$ENTREZID != 100124696,]
  ## Adding the fold Change
  Ag_diff_req_Entrez$logFC <- Ag_diff_req[match(Ag_diff_req_Entrez$SYMBOL, rownames(Ag_diff_req)), "log2FoldChange"]

  ### Saving the file

  # write.csv(Ag_diff_req_Entrez, file = paste(savedir,Ag_name[j],"_all_O_vs_Y_diff_req_Entrez.csv",sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE)

  library(clusterProfiler)
  geneList <- Ag_diff_req_Entrez$logFC ## logfc
  names(geneList) <- as.character(Ag_diff_req_Entrez$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)

  library(fgsea)
  library(data.table)
  library(ggplot2)

  # Filter the gene sets for pathways of interest
  # filtered_gene_sets <- msigdb_c5_bp %>% filter(gs_name %in% pathways_of_interest)

  # # Convert the filtered gene sets to list format
  # gene_sets <- filtered_gene_sets %>% split(x = .$entrez_gene, f = .$gs_name)

  # # Run fgsea
  # fgsea_results <- fgsea(pathways = gene_sets, stats = geneList,
  #                        minSize  = 8,
  #                       eps = 0.0,
  #                       maxSize  = 1500)

  # write.table(as.matrix(fgsea_results), "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/clus_0_4_5_7_8_13_P_vs_R/neuron_differentiation_GCP.txt",
  # row.names = TRUE, col.names = TRUE, sep = "\t", quote=F)

  # pathways_combined <- list()
  rm(pathways)
  pathways <- list()
  T_cells <- c("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/GN_top_markers.txt", "/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/GCP_top_markers.txt")
  cell_name <- gsub("_top_markers.txt", "", basename(T_cells))
  for (k in 1:length(T_cells)) {
    tfh_genes <- read.table(T_cells[k])[, 1]
    Tfh_entrez <- bitr(tfh_genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
    Tfh_entrez[, 2]
    #
    pathways[[cell_name[k]]] <- Tfh_entrez[, 2]
    geneList

    fgseaRes <- fgsea(
      pathways = pathways,
      stats = geneList,
      minSize = 8,
      eps = 0.0,
      maxSize = 1500
    )

    ## If you get the warning it is due to the  genes in your preranked list having the same numerical/ranking value. For example, if you rank genes by log2 fold change,
    # genes with the same log2 fold change value will be ordered randomly once processed. You can use rank() function and a ties.method parameter to control how genes
    # with the same value are ranked:

    ### for this it
    pdf(paste(savedir, Ag_name[j], "_P_vs_R/", cell_name[k], "_GSEA.pdf", sep = ""))
    print(plotEnrichment(pathways[[cell_name[k]]], geneList, gseaParam = 1) + labs(title = cell_name[k]))
    dev.off()
  }

  pdf(paste(savedir, Ag_name[j], "_P_vs_R/GCP_GN_GSEA_table.pdf", sep = ""))
  # pdf(paste(savedir,Ag_name[j],"_",filename_pathway,"_all_O_vs_Y_table.pdf",sep = ""))
  print(plotGseaTable(pathways[cell_name], geneList, fgseaRes, gseaParam = 0.5))
  dev.off()

  fgsea_df <- as.matrix(fgseaRes)
  write.table(fgsea_df, paste(savedir, Ag_name[j], "_P_vs_R/GCP_GN_table_2.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

write.table(fgsea_df, "/diazlab/data3/.abhinav/projects/SHH/snRNA/pseudobulk/P_vs_R/paired/clus_all_P_vs_R/neuron_differentiation.txt", quote = F, row.names = T, col.names = T, sep = "\t")


### Multiple testing
T_cells <- c("GSEA_Exhausted-vs-Effector_CD8-LCMV_GSE9650.txt", "GSEA_Exhausted-vs-Memory_CD8-LCMV_GSE41867.txt", "GSEA_StemLikeExhausted-vs-TerminalExhausted_CD8-LCMV_GSE84105.txt", "GSEA_TerminalExhausted-vs-StemLikeExhausted_CD8-LCMV_GSE84105.txt", "GSEA_Senescence_GOBP.txt", "GSEA_Senescence_Reactome.txt", "GSEA_Senescence_SenMayo.txt", "GSEA_Senescence_SenSig_UP.txt")
pathways <- list()
for (i in 1:length(T_cells)) {
  innate <- read.table(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/Resource/", T_cells[i], sep = ""))[, 1]
  innate_entrez <- bitr(innate, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  pathways[[T_cells[i]]] <- innate_entrez[, 2]
}

fgseaRes <- fgsea(
  pathways = pathways,
  stats = geneList,
  minSize = 8,
  eps = 0.0, maxSize = 1500
)

### for this it
pdf(paste(savedir, Ag_name[j], "_", T_cells[k], "_all_O_vs_Y_2.pdf", sep = ""))
print(plotEnrichment(pathways[["adaptiveness"]], geneList) + labs(title = T_cells[k]))
dev.off()

pdf(paste(savedir, Ag_name[j], "_", T_cells[k], "_all_O_vs_Y_table_2.pdf", sep = ""))
print(plotGseaTable(pathways[T_cells[k]], geneList, fgseaRes, gseaParam = 0.5))
dev.off()

### Find all markers
DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj@meta.data$celltype_marker
celltype_markers <- FindAllMarkers(obj)
write.table(celltype_markers, "/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltype_markers_annotated.txt", quote = F, col.names = T, row.names = T, sep = "\t")

### Performing the differential for P vs R for each celltype
obj@meta.data$tumor_sample_celltype <- paste(obj@meta.data$tumor_type, obj@meta.data$patient, obj@meta.data$celltype_marker, sep = "_")
obj@meta.data$tumor_patient_celltype <- paste(obj@meta.data$tumor_type, obj@meta.data$patient, obj@meta.data$celltype_marker, sep = "_")

remove_samples <- c("DOD4182", "SF7994", "SF8368", "SF8539", "SF10961", "SF9232", "SF12930", "CBTN2978", "CBTN311", "CBTN1666", "CBTN1676")
paired_samples <- grep(paste(remove_samples, collapse = "|"), unique(obj@meta.data$sample_id), value = TRUE, invert = TRUE)
paired_cells <- rownames(obj@meta.data[grep(paste(paired_samples, collapse = "|"), obj@meta.data$sample_id), ])

obj_subset <- subset(obj, cells = paired_cells)
obj_subset <- NormalizeData(obj_subset)
R_vs_P_all <- FindMarkers(obj_subset, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type")
write.table(R_vs_P_all, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_all_lognorm.txt", sep = "\t", row.names = T, col.names = T, quote = F)

obj_subset@meta.data$tumor_celltype <- paste(obj_subset@meta.data$tumor_type, obj_subset@meta.data$celltype_marker, sep = "_")
R_vs_P_GCP <- FindMarkers(obj_subset, ident.1 = "Recurrent_GCP", ident.2 = "Primary_GCP", group.by = "tumor_celltype")
write.table(R_vs_P_GCP, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_GCP_lognorm.txt", sep = "\t", row.names = T, col.names = T, quote = F)

R_vs_P_GNearly <- FindMarkers(obj_subset, ident.1 = "Recurrent_GNearly", ident.2 = "Primary_GNearly", group.by = "tumor_celltype")
write.table(R_vs_P_GNearly, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_GNearly_RNA_lognorm.txt", sep = "\t", row.names = T, col.names = T, quote = F)

R_vs_P_GNlate <- FindMarkers(obj_subset, ident.1 = "Recurrent_GNlate", ident.2 = "Primary_GNlate", group.by = "tumor_celltype")
write.table(R_vs_P_GNlate, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_GNlate_RNA_lognorm.txt", sep = "\t", row.names = T, col.names = T, quote = F)

#### USING MAST
R_vs_P_all <- FindMarkers(obj_subset, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type", test.use = "MAST")
write.table(R_vs_P_all, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_all_MAST.txt", sep = "\t", row.names = T, col.names = T, quote = F)

R_vs_P_GCP <- FindMarkers(obj_subset, ident.1 = "Recurrent_GCP", ident.2 = "Primary_GCP", group.by = "tumor_celltype", test.use = "MAST")
write.table(R_vs_P_GCP, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_GCP_MAST.txt", sep = "\t", row.names = T, col.names = T, quote = F)

R_vs_P_GNearly <- FindMarkers(obj_subset, ident.1 = "Recurrent_GNearly", ident.2 = "Primary_GNearly", group.by = "tumor_celltype", test.use = "MAST")
write.table(R_vs_P_GNearly, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_GNearly_RNA_MAST.txt", sep = "\t", row.names = T, col.names = T, quote = F)

R_vs_P_GNlate <- FindMarkers(obj_subset, ident.1 = "Recurrent_GNlate", ident.2 = "Primary_GNlate", group.by = "tumor_celltype", test.use = "MAST")
write.table(R_vs_P_GNlate, "/diazlab/data3/.abhinav/projects/SHH/snRNA/scdifferential/R_vs_P_GNlate_RNA_MAST.txt", sep = "\t", row.names = T, col.names = T, quote = F)


#### Performing the deconvolution of public available bulk RNAseq data using our snRNA
DefaultAssay(obj) <- "RNA"

obj_cbs <- NormalizeData(
  obj,
  normalization.method = "RC",
  scale.factor = 1e6,
  margin = 1,
  block.size = NULL,
  verbose = TRUE
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
pdf(paste(savedir, "UMAP/obj.pdf", sep = ""))
DimPlot(obj_cbs, reduction = "umap", label.size = 6, label = TRUE)
dev.off()

# obj_cbs@meta.data$seurat_clusters_2 <- paste("clus",obj_cbs@meta.data$seurat_clusters,sep="_")
obj_cbs_df <- obj_cbs@assays$RNA@data

# Set the seed for reproducibility (optional)
set.seed(123)

# Generate 10,000 random numbers from 1 to 140,815
random_numbers <- sample(1:140815, 10000, replace = FALSE)
obj_cbs_df_small <- obj_cbs_df[, random_numbers]
metadata_small <- obj_cbs@meta.data[random_numbers, ]

all(colnames(obj_cbs_df_small) == rownames(metadata_small))
celltypes <- as.data.frame(metadata_small$celltype_marker) %>% t()
rownames(celltypes) <- "GeneSymbol"

write.table(obj_cbs_df_small,
  paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/obj_reference_cpm_small.txt"),
  sep = "\t", quote = F, col.names = F, row.names = T
)

write.table(celltypes, paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/celltypes.txt"),
  sep = "\t", quote = F, col.names = F, row.names = T
)

### Converting it to RPKM
## We have the gene length at this location
obj_cbs_df_small <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/obj_reference_cpm_small.txt", sep = "\t", header = TRUE)
gene_pos <- read.table("/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/final_gene_position_only_chr.txt", header = FALSE)
gene_pos_require <- gene_pos[match(rownames(obj_cbs_df_small), gene_pos$V1, nomatch = 0), ]
obj_cbs_df_small_require <- obj_cbs_df_small[match(gene_pos_require$V1, rownames(obj_cbs_df_small), nomatch = 0), ]
stopifnot(all(gene_pos_require$V1 == rownames(obj_cbs_df_small_require)))
gene_pos_require$gene_length <- (gene_pos_require$V4 - gene_pos_require$V3) / 1000 # Since we need in KB

obj_cbs_df_small_require_rpkm <- obj_cbs_df_small_require

for (i in 1:nrow(gene_pos_require)) {
  obj_cbs_df_small_require_rpkm[i, ] <- obj_cbs_df_small_require[i, ] / gene_pos_require$gene_length[i]
}

write.table(obj_cbs_df_small_require_rpkm,
  paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/obj_reference_cpm_small_rpkm.txt"),
  sep = "\t", quote = F, col.names = F, row.names = T
)

#  Since the size is very big we have to decrease the size
# obj_cbs_df_2 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scRNASeq/Yuki/scRNA_download_Immunity/analysis/Table/scRNA_mps_reference_cpm.txt", header = FALSE)
# obj_cbs_df_small <- sample(obj_cbs_df, size=10000)
# obj_cbs_df_small[,1] <- obj_cbs_df$V1

# write.table(obj_cbs_df_small, paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scRNASeq/Yuki/scRNA_download_Immunity/analysis/Table/scRNA_mps_reference_cpm_small.txt"),
#             sep = "\t", quote = F, col.names = F, row.names = F)

## bulk mixture reference
# library(session)
# restore.session("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Yuki/RNA_aortitis/Change_severe_AO_samples/42_samples.Rda")
# desl_cpm <- cpm(desl, log=FALSE)
# Counts <- desl_cpm

Counts <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/Normalized_gene_counts.txt", header = TRUE, sep = "\t")
rownames(Counts) <- Counts$probeset
Counts <- Counts[, c(-1, -2)]

write.table(Counts, paste("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/Counts_cibersoft.txt"),
  sep = "\t", quote = F, col.names = T, row.names = T
)

write.table(colnames(Counts), paste("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/header_cibersoft.txt"),
  sep = "\t", quote = F, col.names = F, row.names = F
)

### Ran the deconvolution in Cibersortx

public_data <- read.csv("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/decon/all_celltypes/CIBERSORTx_Job2_Results.csv", header = TRUE)
public_data2 <- public_data[, 1:10]

public_meta <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/metadata", header = TRUE)
public_metashh <- public_meta[grep("shh", public_meta$group), ]
public_data2_shh <- public_data2[match(public_data2$Mixture, public_metashh$probeset, nomatch = 0), ]
stopifnot(all(public_metashh$probeset == public_data2_shh$Mixture))

public_data2_shh$status <- public_metashh$status
rownames(public_data2_shh) <- public_data2_shh$Mixture

library(reshape2)
public_data2_shh_melted <- melt(public_data2_shh)

colnames(public_data2_shh_melted) <- c("sample", "status", "celltype", "percentage")
library(ggplot2)

public_data2_shh_melted$sample <- factor(public_data2_shh_melted$sample, levels = unique(public_data2_shh_melted$sample))

p <- ggplot(public_data2_shh_melted, aes(fill = celltype, y = percentage, x = sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/decon/all_celltypes/all_celltype_decon.pdf", width = 10, height = 7)
p
dev.off()

#### Running Deconvolution with RPKM
public_data <- read.csv("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/decon/all_celltypes/CIBERSORTx_Job4_Results_RPKM_Ref.csv", header = TRUE)
public_data2 <- public_data[, 1:10]

public_meta <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/metadata", header = TRUE)
public_metashh <- public_meta[grep("shh", public_meta$group), ]
public_data2_shh <- public_data2[match(public_data2$Mixture, public_metashh$probeset, nomatch = 0), ]
stopifnot(all(public_metashh$probeset == public_data2_shh$Mixture))

public_data2_shh$status <- public_metashh$status
rownames(public_data2_shh) <- public_data2_shh$Mixture

library(reshape2)
public_data2_shh_melted <- melt(public_data2_shh)

colnames(public_data2_shh_melted) <- c("sample", "status", "celltype", "percentage")
library(ggplot2)

public_data2_shh_melted$sample <- factor(public_data2_shh_melted$sample, levels = unique(public_data2_shh_melted$sample))

p <- ggplot(public_data2_shh_melted, aes(fill = celltype, y = percentage, x = sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/decon/all_celltypes/all_celltype_decon_RPKM.pdf", width = 10, height = 7)
p
dev.off()

#####
library(rstatix)
library(ggpubr)
stat.test <- public_data2_shh_melted %>%
  group_by(celltype) %>%
  t_test(percentage ~ status, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

bxp <- ggboxplot(
  public_data2_shh_melted,
  x = "celltype", y = "percentage",
  color = "status", palette = c("#00AFBB", "#E7B800")
)

bxp2 <- bxp + geom_dotplot(
  aes(fill = status, color = status),
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


pdf(paste("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/decon/all_celltypes/t_test_paired_primary_recurrent_public_bxplot_rpkm.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()


### Checking the gene expression in the public data
genes <- c("RPL39", "RPL26", "RPL29", "RPS25", "RPL34", "RPL35A", "RSL24D1", "RPL35A", "RPS18", "RSL24D1", "EIF3E", "EIF2AK4", "GCN2", "EIF2AK3")
genes <- c("AKT1", "FOXO1", "RICTOR", "CDKN1B", "AKT3", "VEGFA", "HSPB1", "CCND2", "TP53", "RPTOR", "PTEN")
genes <- c("EIF3E", "EIF2AK4", "GCN2", "EIF2AK3", "EIF4EBP1", "MAPK8", "FOXO3", "EIF2A", "EIF3E", "ATF4")
genes <- c("HSPA4", "HSPB8", "HSP90B1", "DNAJC10", "DNAJC19", "DNAJC28", "DNAJC19P1", "DNAJC30", "DNAJB11")
genes <- c("NES", "CD34")
genes <- c("FGFR1", "EGFR", "IGF1R", "FLT1", "KDR", "FLT4", "ITGA5", "ITGB1", "ITGA6", "ITGB3", "ITGA2")
genes <- c("MTOR", "RICTOR", "PDK1", "PIK3CA", "AKT2", "AKT3", "FOXO1", "CCND2", "RAPTOR", "ILK", "AKT1", "CDK2", "CDKN2")

counts <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/Counts_cibersoft.txt", header = TRUE, row.names = 1)
genes <- row.names(counts)[match(genes, row.names(counts), nomatch = 0)]
SHH_counts <- log2(counts[genes, 1:50]) ### Only SHH
# SHH_melted <- melt(SHH_counts)
SHH_counts$gene <- rownames(SHH_counts)

# Load necessary libraries
library(tidyverse)
library(tidyr)
library(ggpubr)

# Reshape the data to long format
data_long <- SHH_counts %>%
  tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  mutate(type = ifelse(as.numeric(sub(".*_(\\d+)$", "\\1", sample)) %% 2 == 1, "Primary", "Recurrent"))

# Perform t-tests
t_test_results <- data_long %>%
  group_by(gene) %>%
  summarize(t_test = list(t.test(expression ~ type, data = cur_data())), .groups = "drop") %>%
  mutate(p_value = map_dbl(t_test, "p.value"))

# Merge t-test results with the long data
data_long <- left_join(data_long, t_test_results, by = "gene")

# Create boxplots with t-test results
plot <- ggplot(data_long, aes(x = type, y = expression, color = type)) +
  geom_boxplot() +
  geom_point(aes(color = type)) +
  facet_wrap(~gene, scales = "free") +
  stat_compare_means(aes(label = ..p.format..),
    method = "t.test",
    label.y = min(data_long$expression) + 0.1
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Comparison of Gene Expression between Primary and Recurrent Samples",
    x = "Sample Type", y = "Expression Level"
  )

pdf("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/downstream_receptor_genes_unpaired.pdf")
plot
dev.off()

counts <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/Counts_cibersoft.txt", header = TRUE, row.names = 1)
endo <- counts[c("NES", "CD34", "FLT1"), ]
endo_trans <- t(endo)
cor(endo_trans[, "NES"], endo_trans[, "CD34"])



#### Performing the paired t-test
# genes <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/primary_high_genes_log_0.5.txt", header=FALSE)[,1]
genes <- unique(c("NES", "CD34", "FLT1"))
counts <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/Counts_cibersoft.txt", header = TRUE, row.names = 1)
genes <- row.names(counts)[match(genes, row.names(counts), nomatch = 0)]
SHH_counts <- log2(counts[genes, 1:50]) ### Only SHH
# SHH_melted <- melt(SHH_counts)
SHH_counts$gene <- rownames(SHH_counts)

# Load necessary libraries
library(tidyverse)
library(tidyr)
library(ggpubr)

# Reshape the data to long format with pairs
data_long <- SHH_counts %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  mutate(pair = rep(1:(ncol(SHH_counts) - 1) / 2, each = nrow(SHH_counts)))

# Separate primary and recurrent measurements
data_long <- data_long %>%
  mutate(type = ifelse(as.numeric(sub(".*_(\\d+)$", "\\1", sample)) %% 2 == 1, "Primary", "Recurrent"))

# Perform paired t-tests
t_test_results <- data_long %>%
  group_by(gene) %>%
  summarize(t_test = list(t.test(expression[type == "Primary"], expression[type == "Recurrent"], paired = TRUE)), .groups = "drop") %>%
  mutate(p_value = map_dbl(t_test, "p.value"))

# Merge t-test results with the long data
data_long <- left_join(data_long, t_test_results, by = "gene")

# Calculate y-position for p-value labels
y_max <- data_long %>%
  group_by(gene) %>%
  summarize(y_pos = max(expression))

data_long <- left_join(data_long, y_max, by = "gene")

# Create boxplots with t-test results
plot <- ggplot(data_long, aes(x = type, y = expression, color = type)) +
  geom_boxplot() +
  geom_point(aes(color = type)) +
  facet_wrap(~gene, scales = "free") +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black", position = position_dodge(0.75)) +
  geom_text(aes(x = 1.5, y = (y_pos - 0.5), label = paste0("p = ", format.pval(p_value, 3))), vjust = -1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Comparison of Gene Expression between Primary and Recurrent Samples",
    x = "Sample Type", y = "Expression Level"
  )

pdf("/diazlab/data3/.abhinav/projects/SHH/public_data/primary_vs_recurrent_36635768/HSP_unpaired.pdf", width = 8, height = 12)
plot
dev.off()


#### Quantifying the cell cycling cell percentage in the object
# Load cell cycle genes
cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score cell cycle phases
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Get the cell cycle scores
obj$CT_phase <- paste(obj$celltype_marker, obj$Phase, sep = "_")

p <- DimPlot(obj)
pdf("/diazlab/data3/.abhinav/projects/SHH/snRNA/UMAP/cellcycle.pdf")
p
dev.off()

cycling_table <- table(obj$CT_phase) %>% as.data.frame()

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

#### Identified the high expression of HOX gene in the recurrent cases
## Lets check it
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
dir.create(paste(savedir, "vlnplot/", sep = ""), showWarnings = FALSE)
setwd(savedir)
# file_list = list.files("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/all_run_analysis/removed_unknown/Vlnplot/",
#                        pattern = ".txt")
genes <- c(
  "HOXA9", "HOXA3", "HOXA10", "HOXA10-AS", "HOXA-AS3", "HOXA5", "HOXA-AS2", "HOXA6",
  "PPM1D", "APPBP2", "RPS6KB1", "VMP1", "PTRH2", "TUBD1", "HEATR6", "USP32"
)

for (j in 1:length(genes)) {
  p <- VlnPlot(obj,
    features = genes[j], group.by = "orig.ident", split.by = "celltype_marker",
    idents = c("GCP", "GNearly", "GNlate"), assay = "RNA"
  )
  p$data$ident <- factor(p$data$ident, levels = c(
    "CBTN4529_1_F_P_Run8", "CBTN278_15_F_P_Run10", "CBTN737_0_M_P_Run8", "CBTN931_9_M_P_Run11",
    "CBTN1666_3_F_P_Run11", "CBTN1676_3_F_P_Run11", "CBTN2978_NA_F_P_Run8",
    "CBTN311_NA_F_P_Run10", "UCSFDOD4182_NA_F_P_Run1", "UCSFSF10961_14_F_P_Run2",
    "UCSFSF12930_NA_M_P_Run7", "UCSFSF8368_10_M_P_Run3", "UCSFSF8539_28_F_P_Run6", "UCSFSF9232_53_M_P_Run2",
    "CBTN5881_1_F_R_Run8", "CBTN333_15_F_R_Run8", "CBTN2118_0_M_R_Run8", "CBTN3023_9_M_R_Run9", "UCSFSF7994R_11_F_R_Run2"
  ))
  sample_color <- c(rep("chartreuse2", 4), rep("chartreuse4", 10), rep("dodgerblue1", 4), rep("dodgerblue4", 1))
  colnames(p$data) <- c("genes", "ident", "split")
  q <- ggplot(p$data, aes_string(x = "ident", y = "genes", fill = "ident")) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    facet_grid(split ~ ., scales = "free") +
    NoLegend() +
    theme_cowplot(font_size = 12) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, size = 0.001) +
    scale_fill_manual(values = sample_color) +
    theme(
      legend.position = "none", panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  pdf(paste("vlnplot/", genes[j], "_GCP_GN_split.pdf", sep = ""), width = 10, height = 8)
  print(q)
  dev.off()
}

for (j in 1:length(genes)) {
  p <- VlnPlot(obj, features = genes[j], group.by = "orig.ident", split.by = "celltype_marker", assay = "RNA")
  p$data$ident <- factor(p$data$ident, levels = c(
    "CBTN4529_1_F_P_Run8", "CBTN278_15_F_P_Run10", "CBTN737_0_M_P_Run8", "CBTN931_9_M_P_Run11",
    "CBTN1666_3_F_P_Run11", "CBTN1676_3_F_P_Run11", "CBTN2978_NA_F_P_Run8",
    "CBTN311_NA_F_P_Run10", "UCSFDOD4182_NA_F_P_Run1", "UCSFSF10961_14_F_P_Run2",
    "UCSFSF12930_NA_M_P_Run7", "UCSFSF8368_10_M_P_Run3", "UCSFSF8539_28_F_P_Run6", "UCSFSF9232_53_M_P_Run2",
    "CBTN5881_1_F_R_Run8", "CBTN333_15_F_R_Run8", "CBTN2118_0_M_R_Run8", "CBTN3023_9_M_R_Run9", "UCSFSF7994R_11_F_R_Run2"
  ))
  sample_color <- c(rep("chartreuse2", 4), rep("chartreuse4", 10), rep("dodgerblue1", 4), rep("dodgerblue4", 1))
  colnames(p$data) <- c("genes", "ident", "split")
  q <- ggplot(p$data, aes_string(x = "ident", y = "genes", fill = "ident")) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    facet_grid(split ~ ., scales = "free") +
    NoLegend() +
    theme_cowplot(font_size = 12) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, size = 0.001) +
    scale_fill_manual(values = sample_color) +
    theme(
      legend.position = "none", panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  pdf(paste("vlnplot/", genes[j], "_split_all.pdf", sep = ""), width = 10, height = 8)
  print(q)
  dev.off()
}

for (j in 1:length(genes)) {
  p <- VlnPlot(obj, features = genes[j], group.by = "orig.ident", assay = "RNA")
  p$data$ident <- factor(p$data$ident, levels = c(
    "CBTN4529_1_F_P_Run8", "CBTN278_15_F_P_Run10", "CBTN737_0_M_P_Run8", "CBTN931_9_M_P_Run11",
    "CBTN1666_3_F_P_Run11", "CBTN1676_3_F_P_Run11", "CBTN2978_NA_F_P_Run8",
    "CBTN311_NA_F_P_Run10", "UCSFDOD4182_NA_F_P_Run1", "UCSFSF10961_14_F_P_Run2",
    "UCSFSF12930_NA_M_P_Run7", "UCSFSF8368_10_M_P_Run3", "UCSFSF8539_28_F_P_Run6", "UCSFSF9232_53_M_P_Run2",
    "CBTN5881_1_F_R_Run8", "CBTN333_15_F_R_Run8", "CBTN2118_0_M_R_Run8", "CBTN3023_9_M_R_Run9", "UCSFSF7994R_11_F_R_Run2"
  ))
  sample_color <- c(rep("chartreuse2", 4), rep("chartreuse4", 10), rep("dodgerblue1", 4), rep("dodgerblue4", 1))
  colnames(p$data) <- c("genes", "ident")
  q <- ggplot(p$data, aes_string(x = "ident", y = "genes", fill = "ident")) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    # facet_grid(split~., scales = "free") +
    NoLegend() +
    theme_cowplot(font_size = 12) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, size = 0.001) +
    scale_fill_manual(values = sample_color) +
    theme(
      legend.position = "none", panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  pdf(paste("vlnplot/", genes[j], "_split_combined.pdf", sep = ""), width = 7, height = 5)
  print(q)
  dev.off()
}

#### MAGIC_RNA
savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/"
dir.create(paste(savedir, "vlnplot/", sep = ""), showWarnings = FALSE)
setwd(savedir)
# file_list = list.files("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/all_run_analysis/removed_unknown/Vlnplot/",
#                        pattern = ".txt")
genes <- c(
  "HOXA9", "HOXA3", "HOXA10", "HOXA10-AS", "HOXA-AS3", "HOXA5", "HOXA-AS2", "HOXA6",
  "PPM1D", "APPBP2", "RPS6KB1", "VMP1", "PTRH2", "TUBD1", "HEATR6", "USP32"
)

dir.create(paste("vlnplot/MAGIC_RNA/GCP_GN/", sep = ""), showWarning = FALSE, recursive = TRUE)
dir.create(paste("vlnplot/MAGIC_RNA/celltypes_splitted/", sep = ""), showWarning = FALSE, recursive = TRUE)
dir.create(paste("vlnplot/MAGIC_RNA/combined/", sep = ""), showWarning = FALSE, recursive = TRUE)

DefaultAssay(obj) <- "MAGIC_RNA"
for (j in 1:length(genes)) {
  colors <- c("GCP" = "blue", "GNearly" = "green", "GNlate" = "red") # Customize as needed
  p <- VlnPlot(obj,
    features = genes[j],
    group.by = "orig.ident",
    split.by = "celltype_marker",
    idents = c("GCP", "GNearly", "GNlate"),
    assay = "MAGIC_RNA", cols = colors
  )

  p$data$ident <- factor(p$data$ident, levels = c(
    "CBTN4529_1_F_P_Run8", "CBTN278_15_F_P_Run10", "CBTN737_0_M_P_Run8", "CBTN931_9_M_P_Run11",
    "CBTN1666_3_F_P_Run11", "CBTN1676_3_F_P_Run11", "CBTN2978_NA_F_P_Run8",
    "CBTN311_NA_F_P_Run10", "UCSFDOD4182_NA_F_P_Run1", "UCSFSF10961_14_F_P_Run2",
    "UCSFSF12930_NA_M_P_Run7", "UCSFSF8368_10_M_P_Run3", "UCSFSF8539_28_F_P_Run6", "UCSFSF9232_53_M_P_Run2",
    "CBTN5881_1_F_R_Run8", "CBTN333_15_F_R_Run8", "CBTN2118_0_M_R_Run8", "CBTN3023_9_M_R_Run9", "UCSFSF7994R_11_F_R_Run2"
  ))
  sample_color <- c(rep("chartreuse2", 4), rep("chartreuse4", 10), rep("dodgerblue1", 4), rep("dodgerblue4", 1))
  colnames(p$data) <- c("genes", "ident", "split")
  q <- ggplot(p$data, aes_string(x = "ident", y = "genes", fill = "ident")) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    facet_grid(split ~ ., scales = "free") +
    NoLegend() +
    theme_cowplot(font_size = 12) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, size = 0.001) +
    scale_fill_manual(values = sample_color) +
    theme(
      legend.position = "none", panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  pdf(paste("vlnplot/MAGIC_RNA/GCP_GN/", genes[j], "_GCP_GN_split.pdf", sep = ""), width = 10, height = 8)
  print(q)
  dev.off()
}

for (j in 1:length(genes)) {
  p <- VlnPlot(obj, features = genes[j], group.by = "orig.ident", split.by = "celltype_marker", assay = "RNA")
  p$data$ident <- factor(p$data$ident, levels = c(
    "CBTN4529_1_F_P_Run8", "CBTN278_15_F_P_Run10", "CBTN737_0_M_P_Run8", "CBTN931_9_M_P_Run11",
    "CBTN1666_3_F_P_Run11", "CBTN1676_3_F_P_Run11", "CBTN2978_NA_F_P_Run8",
    "CBTN311_NA_F_P_Run10", "UCSFDOD4182_NA_F_P_Run1", "UCSFSF10961_14_F_P_Run2",
    "UCSFSF12930_NA_M_P_Run7", "UCSFSF8368_10_M_P_Run3", "UCSFSF8539_28_F_P_Run6", "UCSFSF9232_53_M_P_Run2",
    "CBTN5881_1_F_R_Run8", "CBTN333_15_F_R_Run8", "CBTN2118_0_M_R_Run8", "CBTN3023_9_M_R_Run9", "UCSFSF7994R_11_F_R_Run2"
  ))
  sample_color <- c(rep("chartreuse2", 4), rep("chartreuse4", 10), rep("dodgerblue1", 4), rep("dodgerblue4", 1))
  colnames(p$data) <- c("genes", "ident", "split")
  q <- ggplot(p$data, aes_string(x = "ident", y = "genes", fill = "ident")) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    facet_grid(split ~ ., scales = "free") +
    NoLegend() +
    theme_cowplot(font_size = 12) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, size = 0.001) +
    scale_fill_manual(values = sample_color) +
    theme(
      legend.position = "none", panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  pdf(paste("vlnplot/MAGIC_RNA/celltypes_splitted/", genes[j], "_split_all.pdf", sep = ""), width = 10, height = 8)
  print(q)
  dev.off()
}
genes <- c("XIST", "RPS4Y1", "DDX3Y", "EIF1AX", "KDM5D", "ZFY")
for (j in 1:length(genes)) {
  p <- VlnPlot(obj, features = genes[j], group.by = "orig.ident", assay = "MAGIC_RNA")
  p$data$ident <- factor(p$data$ident, levels = c(
    "CBTN4529_1_F_P_Run8", "CBTN278_15_F_P_Run10", "CBTN737_0_M_P_Run8", "CBTN931_9_M_P_Run11",
    "CBTN1666_3_F_P_Run11", "CBTN1676_3_F_P_Run11", "CBTN2978_NA_F_P_Run8",
    "CBTN311_NA_F_P_Run10", "UCSFDOD4182_NA_F_P_Run1", "UCSFSF10961_14_F_P_Run2",
    "UCSFSF12930_NA_M_P_Run7", "UCSFSF8368_10_M_P_Run3", "UCSFSF8539_28_F_P_Run6", "UCSFSF9232_53_M_P_Run2",
    "CBTN5881_1_F_R_Run8", "CBTN333_15_F_R_Run8", "CBTN2118_0_M_R_Run8", "CBTN3023_9_M_R_Run9", "UCSFSF7994R_11_F_R_Run2"
  ))
  sample_color <- c(rep("chartreuse2", 4), rep("chartreuse4", 10), rep("dodgerblue1", 4), rep("dodgerblue4", 1))
  colnames(p$data) <- c("genes", "ident")
  q <- ggplot(p$data, aes_string(x = "ident", y = "genes", fill = "ident")) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    # facet_grid(split~., scales = "free") +
    NoLegend() +
    theme_cowplot(font_size = 12) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, size = 0.001) +
    scale_fill_manual(values = sample_color) +
    theme(
      legend.position = "none", panel.spacing = unit(0, "lines"),
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  pdf(paste("vlnplot/MAGIC_RNA/combined/", genes[j], "_split_combined.pdf", sep = ""), width = 7, height = 5)
  print(q)
  dev.off()
}

genes <- c("OTX2", "OTX2-AS1", "UNC5D", "ROBO3", "EOMES", "LMX1A", "ERBB4")
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(CD4_modality_integrate_cluster, genes[i],
    reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2", size = 0.2
  ) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir, "featureplot/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "featureplot/CD4_Goronzy.pdf", sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()


### Making nice Figures
GCP <- c("PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L", "ATOH1", "GLI2", "GLI1", "PTCH1")
GNearly <- c("NHLH2", "NHLH1", "STMN2", "TUBB3", "NRG1", "GRIN2B", "KIF21B", "DISP3")
GNlate <- c("NEUROD1", "NEUROD2", "OXR1", "NRXN2", "RBFOX3", "BRINP1", "GRIN2C")
Neuronal <- c("TRPM3", "OTX2", "ROBO3", "UNC5D", "ITGA3", "ADRA1A", "OTX2-AS1", "PRDM6")
Stem_cells <- c("NES", "SOX2", "MKI67", "TOP2A")
HSP <- c("HSPH1", "HSPA1B", "HSPA1A", "DNAJB1", "HSPB1", "VEGFA", "HIF1A")
RB <- c("RPL35", "RPS5", "RPL4", "RPL26", "RPL29", "RPS17")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[, 1]
genes <- c(GCP, Stem_cells, HSP, RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

p <- DoHeatmap(RNA_integrated, genes, group.by = "sub_celltypes", assay = "MAGIC_RNA")

# dir.create(paste(savedir, "/heatmap/MAGIC_RNA", sep = ""), showWarnings = FALSE, recursive = TRUE)
pdf(paste(savedir, "/heatmap/MAGIC_RNA/all_heatmap_subcelltypes.pdf", sep = ""), width = 12, height = 24)
print(p)
dev.off()

p <- DoHeatmap(RNA_integrated, genes, group.by = "sub_celltypes", assay = "RNA")

dir.create(paste(savedir, "/heatmap/RNA", sep = ""), showWarnings = FALSE, recursive = TRUE)
pdf(paste(savedir, "/heatmap/RNA/all_subcelltypes_RNA.pdf", sep = ""), width = 12, height = 24)
print(p)
dev.off()


#### Making a dotplot
p <- DotPlot(snRNA, genes, group.by = "sub_celltypes", assay = "RNA") +
  scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
  coord_flip() +
  scale_size(range = c(1, 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/SHH_marker_genes_subcelltypes.pdf"), height = 20, width = 10)
print(p)
dev.off()


SHH_malignant_cells <- rownames(snRNA@meta.data[grep("GCP|GN", snRNA@meta.data$celltypes), ])
write.table(SHH_malignant_cells,
  paste0("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/SHH_malignant_cells.txt"),
  quote = F, row.names = F, col.names = F, sep = "\t"
)

### Making a small figure
GCP <- c("GLI1", "GLI2", "PTCH1")
GCP_cycling <- c("MKI67", "TOP2A")
GCP_ribosomal <- c("RPL35", "RPS5", "RPL4")
GCP_stress <- c("HSPH1", "HSPA1B", "DNAJB1")
GN_pre <- c("NHLH1", "NHLH2", "STMN2")
GN_mig <- c("NEUROD1", "NEUROD2", "BRINP1")
GN_post <- c("NRXN2", "RBFOX3", "GRIN2C")
Astro <- c("GFAP", "AQP4", "SLC6A11")
Endo <- c("VWF", "CD34", "CD93")
Immune <- c("CD14", "CD163", "CD4")
Oligo <- c("MOBP", "PLP1", "MOG")
OPC <- c("MYT1", "VCAN", "CSPG4")
Pericytes <- c("PDGFRB", "ACTA2", "RGS5")

genes <- c(GCP, GCP_cycling, GCP_ribosomal, GCP_stress, GN_pre, GN_mig, GN_post, Astro, Endo, Immune, Oligo, OPC, Pericytes)

# snRNA@meta.data$cell_types <- snRNA@meta.data$sub_celltypes
# snRNA@meta.data$cell_types <- gsub("GCP_cycling_1","GCP_cycling",snRNA@meta.data$cell_types) %>%
# gsub("GCP_cycling_2","GCP_ribo",.) %>%
# gsub("GCP_HSP","GCP_stress",.) %>%
# gsub("Neuron_other","Neuron",.) %>%
# gsub("RL_like","GCP_GN",.) %>%
# gsub("GN_Premigratory","GN_premigratory",.)

snRNA@meta.data$cell_types <- factor(snRNA@meta.data$cell_types, levels = c("GCP", "GCP_cycling", "GCP_ribo", "GCP_stress", "GCP_GN", "GN_cycling", "GN_premigratory", "GN_migratory", "GN_postmigratory", "Neuron", "Astro", "Endo", "Immune", "Oligo", "OPC", "Pericyte"))

p <- DotPlot(snRNA, genes, group.by = "cell_types", assay = "RNA") +
  scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
  coord_flip() +
  scale_size(range = c(1, 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/"

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/SHH_marker_genes_cell_types_small.pdf"), height = 12, width = 9)
print(p)
dev.off()
