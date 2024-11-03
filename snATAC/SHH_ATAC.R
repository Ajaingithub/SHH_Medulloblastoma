library(dplyr)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86) # hg38
library(GenomicRanges)
library(future)

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

ATAC_loc <- list.files("/diazlab/data3/SHH/snATACseq/cellranger/", pattern = "filtered_peak_bc_matrix.h5", recursive = TRUE, full.names = TRUE)
datadir <- gsub("outs/filtered_peak_bc_matrix.h5", "", ATAC_loc)
samplename <- gsub(".*.batch[1-3]/", "", datadir) %>%
  gsub("_ATAC", "", .) %>%
  gsub("/", "", .) %>%
  gsub("-", "_", .)
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/first/"

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) ## Adding Annotation since it is same for all the objects

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "hg38"

for (i in 1:length(samplename)) {
  source("/diazlab/data3/.abhinav/resources/functions/R/scATACseq/QC_plots.R")
  obj <- ATAC_QCplots(datadir[i], savedir, samplename[i], annotations)
  assign(paste0("obj_", samplename[i]), obj)
}

## Subsetting the value keeping the minimum cutoff
object <- ls(pattern = "^obj_") %>% grep("_subset$", ., invert = TRUE, value = TRUE, )
rm(ATAC_df)
ATAC_df <- as.data.frame(matrix(ncol = 4, nrow = length(object)))
colnames(ATAC_df) <- c("samplename", "before", "after_minimal", "after_default")
ATAC_df$samplename <- object

for (i in 1:length(object)) {
  try({
    obj <- get(object[i])
    ATAC_df[i, "before"] <- ncol(obj)

    obj_subset <- subset(
      x = obj,
      subset = nCount_peaks > 1000 &
        nCount_peaks < 40000 &
        pct_reads_in_peaks > 10 &
        blacklist_ratio < 0.1 &
        nucleosome_signal < 6 &
        TSS.enrichment > 2
    )

    ATAC_df[i, "after_minimal"] <- ncol(obj_subset)

    # obj_subset2 <- subset(
    #   x = obj,
    #   subset = nCount_peaks > 3000 &
    #     nCount_peaks < 30000 &
    #     pct_reads_in_peaks > 15 &
    #     blacklist_ratio < 0.05 &
    #     nucleosome_signal < 4 &
    #     TSS.enrichment > 3
    # )

    # ATAC_df[i, "after_default"] <- ncol(obj_subset2)

    assign(paste0(object[i], "_subset"), obj_subset)
  })
}

dir.create(paste0(savedir, "Table"), showWarnings = FALSE)
write.table(ATAC_df, paste0(savedir, "Table/QC_cutoff.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Normalization and linear dimensional reduction

# Normalization: Signac performs term frequency-inverse document frequency (TF-IDF) normalization. This is a two-step normalization procedure, that both
# normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.

# Feature selection: The low dynamic range of scATAC-seq data makes it challenging to perform variable feature selection, as we do for scRNA-seq. Instead, we
# can choose to use only the top n% of features (peaks) for dimensional reduction, or remove features present in less than n cells with the FindTopFeatures()
# function. Here we will use all features, though we have seen very similar results when using only a subset of features (try setting min.cutoff to ‘q75’ to
# use the top 25% all peaks), with faster runtimes. Features used for dimensional reduction are automatically set as VariableFeatures() for the Seurat object
# by this function.

# Dimension reduction: We next run singular value decomposition (SVD) on the TD-IDF matrix, using the features (peaks) selected above. This returns a reduced
#  dimension representation of the object (for users who are more familiar with scRNA-seq, you can think of this as analogous to the output of PCA).

# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI), and were first introduced for the analysis of scATAC-seq data by Cusanovich et al. 2015.

objname <- grep("7316", ls(pattern = "_subset"), value = TRUE)
dir_name <- gsub("obj_", "", objname) %>% gsub("_subset", "", .)
for (i in 2:length(objname)) {
  pbmc <- get(objname[i])
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
  pbmc <- RunSVD(pbmc)

  dir.create(paste0(savedir, dir_name[i], "/QC/"), showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(savedir, dir_name[i], "/QC/", objname[i], "_Depth_Corr.pdf"))
  print(DepthCor(pbmc))
  dev.off()

  pbmc <- RunUMAP(object = pbmc, reduction = "lsi", dims = 2:30)
  pbmc <- FindNeighbors(object = pbmc, reduction = "lsi", dims = 2:30)
  pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

  dir.create(paste0(savedir, dir_name[i], "/UMAP/"), showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(savedir, dir_name[i], "/UMAP/", objname[i], "_LSI_default.pdf"))
  print(DimPlot(object = pbmc, label = TRUE) + NoLegend())
  dev.off()

  assign(objname[i], pbmc)
}

####################################
# Create the Gene Activity Score
####################################

#  We can try to quantify the activity of each gene in the genome by assessing the chromatin accessibility associated with the gene, and create a new gene
# activity assay derived from the scATAC-seq data. Here we will use a simple approach of summing the fragments intersecting the gene body and promoter region
# To create a gene activity matrix, we extract gene coordinates and extend them to include the 2 kb upstream region (as promoter accessibility is often correlated
# with gene expression). We then count the number of fragments for each cell that map to each of these regions, using the using the FeatureMatrix() function.
# These steps are automatically performed by the GeneActivity() function:

# Now we can visualize the activities of canonical marker genes to help interpret our ATAC-seq clusters. Note that the activities will be much noisier than scRNA-seq
# measurements. This is because they represent measurements from sparse chromatin data, and because they assume a general correspondence between gene body/promoter
# accessibility and gene expression which may not always be the case. Nonetheless, we can begin to discern populations of monocytes, B, T, and NK cells based on these
# gene activity profiles. However, further subdivision of these cell types is challenging based on supervised analysis alone.

objname <- grep("7316", ls(pattern = "_subset"), value = TRUE)
dir_name <- gsub("obj_", "", objname) %>% gsub("_subset", "", .)

### Since we will be performing the cell type annotation as well
GCP <- c("PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L", "ATOH1", "GLI2", "GLI1", "PTCH1")
GNearly <- c("NHLH2", "NHLH1", "STMN2", "TUBB3", "NRG1", "GRIN2B", "KIF21B", "DISP3")
GNlate <- c("NEUROD1", "NEUROD2", "OXR1", "NRXN2", "RBFOX3", "BRINP1", "GRIN2C")
Neuronal <- c("TRPM3", "OTX2", "ROBO3", "UNC5D", "ITGA3", "ADRA1A", "OTX2-AS1", "PRDM6")
Stem_cells <- c("NES", "SOX2", "MKI67", "TOP2A")
HSP <- c("HSPH1", "HSPA1B", "HSPA1A", "DNAJB1", "HSPB1")
RB <- c("RPL35", "RPS5", "RPL4", "RPL26", "RPL29", "RPS17")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[, 1]
genes <- c(GCP, Stem_cells, HSP, RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

for (i in 1:length(objname)) {
  try({
    pbmc <- get(objname[i])
    gene.activities <- GeneActivity(pbmc)

    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)
    pbmc <- NormalizeData(
      object = pbmc,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = median(pbmc$nCount_RNA)
    )

    DefaultAssay(pbmc) <- "RNA"

    p <- FeaturePlot(
      object = pbmc,
      features = c("STMN2", "GLI2", "NEUROD1", "NES", "HSPH1", "RPL35"),
      pt.size = 0.1,
      max.cutoff = "q95",
      ncol = 3
    )

    dir.create(paste0(savedir, dir_name[i], "/featureplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, dir_name[i], "/featureplot/", objname[i], "_marker_genes.pdf"), width = 10, height = 8)
    print(p)
    dev.off()


    p <- DotPlot(pbmc, genes, assay = "RNA") + scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
      coord_flip() +
      # scale_size(range = c(1, 10)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, dir_name[i], "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, dir_name[i], "/dotplot/", objname[i], "_marker_genes_noscale_range.pdf"), height = 26, width = 9)
    print(p)
    dev.off()

    assign(objname[i], pbmc)

    ## Saving the Object
    dir.create(paste0(savedir, dir_name[i], "/saveRDS_obj/"))
    saveRDS(pbmc, paste0(savedir, dir_name[i], "/saveRDS_obj/", objname[i], "_QC.RDS"))
  })
}

# Integrating with scRNA-seq data
# To help interpret the scATAC-seq data, we can classify cells based on an scRNA-seq experiment from the same biological system (human PBMC). We utilize methods for
# cross-modality integration and label transfer, described https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub, with a more in-depth
# tutorial here. We aim to identify shared correlation patterns in the gene activity matrix and scRNA-seq dataset to identify matched biological states across
# the two modalities. This procedure returns a classification score for each cell for each scRNA-seq-defined cluster label.

# Here we load a pre-processed scRNA-seq dataset for human PBMCs, also provided by 10x Genomics. You can download the raw data for this experiment from the 10x website,
# and view the code used to construct this object on GitHub. Alternatively, you can download the pre-processed Seurat object here.

### Loading the snRNA seq data
snRNA_obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed.RDS")
DefaultAssay(snRNA_obj) <- "RNA"
snRNA_obj <- FindVariableFeatures(snRNA_obj, nfeatures = 3000)

# Reference mapping
objname <- grep("7316", ls(pattern = "_subset"), value = TRUE)
dir_name <- gsub("obj_", "", objname) %>% gsub("_subset", "", .)

for (i in 3:length(objname)) {
  pbmc <- get(objname[i])
  DefaultAssay(pbmc) <- "RNA"
  pbmc <- FindVariableFeatures(pbmc, nfeatures = 4000)

  transfer.anchors <- FindTransferAnchors(
    reference = snRNA_obj,
    query = pbmc,
    reduction = "cca"
  )

  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = snRNA_obj$sub_celltypes,
    weight.reduction = pbmc[["lsi"]],
    dims = 2:30
  )

  predicted.labels2 <- TransferData(
    anchorset = transfer.anchors,
    refdata = snRNA_obj$celltypes,
    weight.reduction = pbmc[["lsi"]],
    dims = 2:30
  )

  pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

  pbmc@meta.data$predicted_celltypes <- predicted.labels2$predicted.id

  plot1 <- DimPlot(
    object = pbmc,
    group.by = "predicted_celltypes",
    label = TRUE,
    repel = TRUE
  ) + ggtitle("celltypes")

  plot2 <- DimPlot(
    object = pbmc,
    group.by = "predicted.id",
    label = TRUE,
    repel = TRUE
  ) + ggtitle("subcelltypes")

  pdf(paste0(savedir, dir_name[i], "/UMAP/", objname[i], "_projection_CCA.pdf"), width = 12, height = 5.5)
  print(plot1 + plot2)
  dev.off()

  assign(objname[i], pbmc)

  saveRDS(pbmc, paste0(savedir, dir_name[i], "/saveRDS_obj/", objname[i], "_QC_snRNA_ref_mapped.RDS"))
}

#### Also performing without CCA
# snRNA_obj <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed.RDS")
DefaultAssay(snRNA_obj) <- "RNA"
snRNA_obj <- FindVariableFeatures(snRNA_obj, nfeatures = 3000)

# Reference mapping
RDS_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/first/", pattern = "_subset_QC.RDS", recursive = TRUE, full.names = TRUE)

objname <- gsub(".*./obj", "obj", RDS_files) %>% gsub("_QC.RDS", "", .)
dir_name <- gsub("obj_", "", objname) %>% gsub("_subset", "", .)

for (i in 9:length(RDS_files)) {
  try({
    pbmc <- readRDS(RDS_files[i])
    DefaultAssay(pbmc) <- "RNA"
    pbmc <- FindVariableFeatures(pbmc, nfeatures = 3000)

    transfer.anchors <- FindTransferAnchors(
      reference = snRNA_obj,
      query = pbmc,
      # reduction = "cca"
    )

    predicted.labels <- TransferData(
      anchorset = transfer.anchors,
      refdata = snRNA_obj$sub_celltypes,
      weight.reduction = pbmc[["lsi"]],
      dims = 2:30
    )

    pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

    plot2 <- DimPlot(
      object = pbmc,
      group.by = "predicted.id",
      label = TRUE,
      repel = TRUE
    ) + ggtitle("scATAC-seq")

    pdf(paste0(savedir, dir_name[i], "/UMAP/", objname[i], "_projection_PCA.pdf"), width = 7, height = 5.5)
    print(plot2)
    dev.off()

    assign(objname[i], pbmc)

    saveRDS(pbmc, paste0(savedir, dir_name[i], "/saveRDS_obj/", objname[i], "_QC_snRNA_ref_mapped_PCA.RDS"))
  })
}

objname <- ls(pattern = "_subset") %>% grep("obj_7316", ., value = TRUE)
samplename <- gsub("obj_|_subset", "", objname)
for (i in 1:length(objname)) {
  pbmc <- get(objname[i])
  pdf(paste(savedir, samplename[i], "/QC/", samplename[i], "_Vlnplot_allQC_clustering_default.pdf", sep = ""), width = 20, height = 10)
  print(VlnPlot(
    object = pbmc,
    features = c(
      "nCount_peaks", "nFeature_peaks", "TSS.enrichment", "blacklist_ratio",
      "nucleosome_signal", "pct_reads_in_peaks", "nCount_RNA",
      "nFeature_RNA"
    ),
    pt.size = 0.1,
    ncol = 4,
  ))
  dev.off()
}

for (i in 1:length(objname)) {
  pbmc <- get(objname[i])
  df1 <- as.data.frame((table(pbmc@meta.data$predicted.id)))
  assign(paste0(objname[i], "_df_CCA"), df1)

  df2 <- as.data.frame(table(pbmc@meta.data$predicted_celltypes))
  assign(paste0(objname[i], "_df_PCA"), df2)
}

### Merging and Integrating
library(purrr)
PCA_df <- ls(pattern = "_df_PCA")

list_of_dfs <- list(
  get(PCA_df[1]), get(PCA_df[3]), get(PCA_df[4]),
  get(PCA_df[5]), get(PCA_df[6]), get(PCA_df[7]), get(PCA_df[8]),
  get(PCA_df[9]), get(PCA_df[10]), get(PCA_df[11]), get(PCA_df[12]), get(PCA_df[13])
)

# Merge all data frames together
merged_df <- reduce(list_of_dfs, full_join, by = "Var1")

required_col <- (gsub("obj_", "", PCA_df) %>% gsub("_subset_df_PCA", "", .))[c(1, 3:13)]
colnames(merged_df) <- c("celltypes", required_col)

# ColumnSum <- colSums(merged_df)
# df_CC <- (t(cellcycling)/ColumnSum)*100
df_melted <- melt(merged_df)
colnames(df_melted) <- c("Celltypes", "Sample", "cellnumber")
p <- ggplot(df_melted, aes(fill = Celltypes, y = cellnumber, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "yellow", "black", "blanchedalmond", "blue", "white"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/first/projection/PCA/celltype_cellnumber.pdf", width = 6, height = 5)
p
dev.off()

library(purrr)
PCA_df <- ls(pattern = "_df_CCA")

list_of_dfs <- list(
  get(PCA_df[1]), get(PCA_df[3]), get(PCA_df[4]),
  get(PCA_df[5]), get(PCA_df[6]), get(PCA_df[7]), get(PCA_df[8]),
  get(PCA_df[9]), get(PCA_df[10]), get(PCA_df[11]), get(PCA_df[12]), get(PCA_df[13])
)

# Merge all data frames together
merged_df <- reduce(list_of_dfs, full_join, by = "Var1")

required_col <- (gsub("obj_", "", PCA_df) %>% gsub("_subset_df_CCA", "", .))[c(1, 3:13)]
colnames(merged_df) <- c("celltypes", required_col)

# ColumnSum <- colSums(merged_df)
# df_CC <- (t(cellcycling)/ColumnSum)*100
df_melted <- melt(merged_df)
colnames(df_melted) <- c("Celltypes", "Sample", "cellnumber")
p <- ggplot(df_melted, aes(fill = Celltypes, y = cellnumber, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(
    "#279e68", "#d62728", "#ff7f0e", "#1f77b4", "#aa40fc",
    "#8c564b", "#e377c2", "#b5bd61", "#17becf",
    "#aec7e8", "#ffbb78", "darkseagreen3", "cornsilk4", "plum2",
    "green3", "black", "blanchedalmond", "blue", "white"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/first/projection/PCA/subcelltype_cellnumber.pdf", width = 8, height = 6)
p
dev.off()

#####
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

merged_df[is.na(merged_df)] <- 0
ColumnSum <- colSums(merged_df[, 2:ncol(merged_df)])
merged_df_percent <- (t(merged_df[, 2:ncol(merged_df)]) / ColumnSum) * 100
colnames(merged_df_percent) <- merged_df$celltypes
df_melted <- melt(merged_df_percent)

patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns)
patterns <- gsub("-", "_", patterns)
replacements <- long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$Sample <- gsub("_12.*.", "", df_melted$Var1)
df_melted$timepts <- df_melted$Sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

### Patient ID
patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns)
patterns <- gsub("-", "_", patterns)
replacements <- long_sample$Patient.ID
# names(replacement) <- patterns
df_melted$Patient.ID <- df_melted$Sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$Patient.ID <- str_replace_all(df_melted$Patient.ID, pattern, replacements[i])
}

# df_melted$celltype2 <- gsub("GCP_cycling", "GCP", df_melted$celltype) %>% gsub("GNearly|GNlate", "GN", .)

colnames(df_melted)[1:3] <- c("sample.id", "celltype", "percentage")
df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931", df_melted$Sample), ]
df_melted_paired$Patient_Timepts <- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts, sep = "_")

df_wide <- dcast(df_melted_paired, celltype ~ Patient_Timepts, value.var = "percentage", fun.aggregate = mean)

df_melted_paired <- df_melted
stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test$p <- (stat.test$p) / 2

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

pdf(paste(savedir, "Table/t_test_paired_primary_recurrent_samples_bxplot_two_sided_sub_celltype_all_samples.pdf", sep = ""), width = 8, height = 6)
bxp4
dev.off()

# When merging multiple single-cell chromatin datasets, it’s important to be aware that if peak calling was performed on each dataset independently, the peaks are unlikely
#  to be exactly the same. We therefore need to create a common set of peaks across all the datasets to be merged.

# To create a unified set of peaks we can use functions from the GenomicRanges package. The reduce function from GenomicRanges will merge all intersecting peaks. Another
# option is to use the disjoin function, that will create distinct non-overlapping sets of peaks.

# Creating a common peak set
# If the peaks were identified independently in each experiment then they will likely not overlap perfectly. We can merge peaks from all the datasets to create a common
# peak set, and quantify this peak set in each experiment prior to merging the objects.

# First we’ll load the peak coordinates for each experiment and convert them to genomic ranges, the use the GenomicRanges::reduce function to create a common set of peaks to quantify in each dataset.

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

### Creating a common peak set
