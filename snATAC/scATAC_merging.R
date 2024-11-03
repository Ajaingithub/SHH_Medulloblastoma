### Merging Object
# Creating a common peak set
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(dplyr)
library(EnsDb.Hsapiens.v86) # hg38

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
## Getting a filtered bed files
### Removing sample 5452 since it has very low quality cells do not want to use this sample for the common as we will remove in the downstream
peak_loc <- list.files("/diazlab/data3/SHH/snATACseq/cellranger/", pattern = "peaks.bed$", recursive = TRUE, full.names = TRUE)
peak_loc2 <- grep("filtered_peak_bc_matrix", peak_loc, value = TRUE) %>%
    grep("7316-5452", ., invert = TRUE, value = TRUE)
samplename <- gsub("_ATAC/outs/filtered_peak_bc_matrix/peaks.bed", "", peak_loc2) %>%
    gsub(".*.batch[1-3]/", "", .) %>%
    gsub("-", "_", .) %>%
    paste0("gr_", .)

for (i in 1:length(samplename)) {
    peaks <- read.table(
        file = peak_loc2[i],
        col.names = c("chr", "start", "end")
    )

    # convert to genomic ranges
    gr <- makeGRangesFromDataFrame(peaks)

    ## Assiging to each samples
    assign(samplename[i], gr)
}

# Create a unified set of peaks to quantify in each dataset
grnames <- ls(pattern = "gr_")
combined.peaks <- Signac::reduce(x = c(
    get(gr_names[1]), get(gr_names[2]), get(gr_names[3]), get(gr_names[4]),
    get(gr_names[5]), get(gr_names[6]), get(gr_names[7]), get(gr_names[8]),
    get(gr_names[9]), get(gr_names[10]), get(gr_names[11]), get(gr_names[12]),
    get(gr_names[13])
))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks2 <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks2

## Create Fragment objects
# To quantify our combined set of peaks we’ll need to create a Fragment object for each experiment. The Fragment class is a specialized class defined in Signac to hold
# all the information related to a single fragment file.

# First we’ll load the cell metadata for each experiment so that we know what cell barcodes are contained in each file, then we can create Fragment objects using the
# CreateFragmentObject function. The CreateFragmentObject function performs some checks to ensure that the file is present on disk and that it is compressed and indexed,
# computes the MD5 sum for the file and the tabix index so that we can tell if the file is modified at any point, and checks that the expected cells are present in the file.

md_loc <- list.files("/diazlab/data3/SHH/snATACseq/cellranger", pattern = "singlecell.csv", recursive = TRUE, full.names = TRUE)
md_loc2 <- grep("outs", md_loc, value = TRUE) %>%
    grep("7316-5452", ., invert = TRUE, value = TRUE)

mdname <- gsub("_ATAC/outs/singlecell.csv", "", md_loc2) %>%
    gsub(".*.batch[1-3]/", "", .) %>%
    gsub("-", "_", .) %>%
    paste0("md_", .)

fragment_loc <- list.files("/diazlab/data3/SHH/snATACseq/cellranger", pattern = "fragments.tsv.gz$", recursive = TRUE, full.names = TRUE)
fragment_loc2 <- grep("outs", fragment_loc, value = TRUE) %>%
    grep("7316-5452", ., invert = TRUE, value = TRUE)

objname <- gsub("_ATAC/outs/fragments.tsv.gz", "", fragment_loc2) %>%
    gsub(".*.batch[1-3]/", "", .) %>%
    gsub("-", "_", .) %>%
    paste0("obj_", .)

for (i in 1:length(fragmentname)) {
    mdfile <- read.table(
        file = md_loc2[i],
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ]

    # perform an initial filtering of low count cells
    mdfile <- mdfile[mdfile$passed_filters > 500, ]

    # create fragment objects
    fragment <- CreateFragmentObject(
        path = fragment_loc2[i],
        cells = rownames(mdfile)
    )

    # Quantify peaks in each dataset
    # We can now create a matrix of peaks x cell for each sample using the FeatureMatrix function. This function is parallelized using the
    # future package.
    SHH_counts <- FeatureMatrix(
        fragments = fragment,
        features = combined.peaks,
        cells = rownames(mdfile)
    )

    # Create the objects
    # We will now use the quantified matrices to create a Seurat object for each dataset, storing the Fragment object for each dataset in the assay.
    SHH_assay <- CreateChromatinAssay(SHH_counts, fragments = fragment)
    obj <- CreateSeuratObject(SHH_assay, assay = "ATAC", meta.data = mdfile)

    # assigning the objname
    assign(paste0(objname[i]), obj)
}

### Cleaning the object before merging it
# extract gene annotations from EnsDb
objname <- ls(pattern = "obj_")
samplenames <- gsub("obj_", "", objname)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) ## Adding Annotation since it is same for all the objects

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "hg38"

for (i in 1:length(samplenames)) {
    source("/diazlab/data3/.abhinav/resources/functions/R/scATACseq/QC_seu_obj.r")
    obj <- ATAC_QC_wid_obj(get(objname[i]), savedir, samplenames[i], annotations)
    assign(paste0("obj_", samplenames[i]), obj)
}

### Subsetting
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
            subset = nCount_ATAC > 1000 &
                nCount_ATAC < 40000 &
                pct_reads_in_peaks > 10 &
                blacklist_ratio < 0.1 &
                nucleosome_signal < 6 &
                TSS.enrichment > 2
        )

        ATAC_df[i, "after_minimal"] <- ncol(obj_subset)
        assign(paste0(object[i], "_subset_min"), obj_subset)

        obj_subset2 <- subset(
            x = obj,
            subset = nCount_ATAC > 3000 &
                nCount_ATAC < 30000 &
                pct_reads_in_peaks > 15 &
                blacklist_ratio < 0.05 &
                nucleosome_signal < 4 &
                TSS.enrichment > 3
        )

        ATAC_df[i, "after_default"] <- ncol(obj_subset2)
        assign(paste0(object[i], "_subset_def"), obj_subset)
    })
}

objname <- ls(pattern = "_subset_min")
for (i in 1:length(samplenames)) {
    obj <- get(objname[i])
    obj$dataset <- paste0("sample_", samplenames[i])
    assign(objname[i], obj)
}

# compute LSI for each of the object
for (i in 2:length(samplenames)) {
    obj <- get(objname[i])
    obj <- FindTopFeatures(obj, min.cutoff = 10)
    obj <- RunTFIDF(obj)
    obj <- RunSVD(obj)
    assign(objname[i], obj)
}

# Need to rename the cells to avoid the name conflict
for (i in 1:length(samplenames)) {
    obj <- get(objname[i])
    new_cells <- paste(samplenames[i], rownames(obj@meta.data), sep = "_")
    obj <- RenameCells(obj, new.names = new_cells)
    assign(paste0(objname[i], "_2"), obj)
}

objname <- ls(pattern = "min_2")
# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
    x = get(objname[1]),
    y = list(
        get(objname[2]), get(objname[3]), get(objname[4]),
        get(objname[5]), get(objname[6]), get(objname[7]),
        get(objname[8]), get(objname[9]), get(objname[10]),
        get(objname[11]), get(objname[12]), get(objname[13])
    )
    # add.cell.ids = samplenames
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = "lsi")

dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/combined_sample.pdf"))
DimPlot(combined, group.by = "dataset")
dev.off()

pdf(paste0(savedir, "UMAP/combined_sample_splitted.pdf"), width = 15, height = 12)
DimPlot(combined, group.by = "dataset", split.by = "dataset", ncol = 4) + NoLegend()
dev.off()

dir.create(paste0(savedir, "saveRDS_obj"), showWarnings = FALSE, recursive = TRUE)
saveRDS(combined, paste0(savedir, "saveRDS_obj/combined.RDS"))

for (i in 2:length(objname)) {
    dir.create(paste0(savedir, samplenames[i], "/saveRDS_obj"), showWarnings = FALSE, recursive = TRUE)
    saveRDS(get(objname[1]), paste0(savedir, samplenames[i], "/saveRDS_obj/", objname[i], ".RDS"))
}

# Integration

# To find integration anchors between the two datasets, we need to project them into a shared low-dimensional space. To do this, we’ll use reciprocal LSI
# projection (projecting each dataset into the others LSI space) by setting reduction="rlsi". For more information about the data integration methods in Seurat,
# see our recent paper and the Seurat website.

# Rather than integrating the normalized data matrix, as is typically done for scRNA-seq data, we’ll integrate the low-dimensional cell embeddings (the LSI
# coordinates) across the datasets using the IntegrateEmbeddings() function. This is much better suited to scATAC-seq data, as we typically have a very sparse
# matrix with a large number of features. Note that this requires that we first compute an uncorrected LSI embedding using the merged dataset (as we did above).
# find integration anchors

integration.anchors <- FindIntegrationAnchors(
    object.list = list(
        get(objname[1]),
        get(objname[2]), get(objname[3]), get(objname[4]),
        get(objname[5]), get(objname[6]), get(objname[7]),
        get(objname[8]), get(objname[9]), get(objname[10]),
        get(objname[11]), get(objname[12]), get(objname[13])
    ),
    anchor.features = rownames(get(objname[1])),
    reduction = "rlsi",
    dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = combined[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindNeighbors(object = integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)

pdf(paste0(savedir, "UMAP/integrated_sample_2.pdf"), width = 6.5, height = 5)
DimPlot(integrated, group.by = "dataset")
dev.off()

pdf(paste0(savedir, "UMAP/integrated_cluster.pdf"), width = 6.5, height = 5)
DimPlot(integrated, label = TRUE)
dev.off()

pdf(paste0(savedir, "UMAP/integrated_sample_splitted.pdf"), width = 15, height = 12)
DimPlot(integrated, group.by = "dataset", split.by = "dataset", ncol = 4) + NoLegend()
dev.off()

dir.create(paste0(savedir, "saveRDS_obj"), showWarnings = FALSE, recursive = TRUE)
saveRDS(combined, paste0(savedir, "saveRDS_obj/combined.RDS"))

dir.create(paste0(savedir, "saveRDS_obj"), showWarnings = FALSE, recursive = TRUE)
saveRDS(integrated, paste0(savedir, "saveRDS_obj/integrated.RDS"))

### Identifying the gene activity score
gene.activities <- GeneActivity(integrated)

integrated[["RNA"]] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
    object = integrated,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = median(integrated$nCount_RNA)
)

integrated <- ScaleData(integrated, assay = "RNA")
DefaultAssay(integrated) <- "RNA"

p <- FeaturePlot(
    object = integrated,
    features = c("STMN2", "GLI2", "NEUROD1", "NES", "HSPH1", "RPL35"),
    pt.size = 0.1,
    max.cutoff = "q95",
    ncol = 3
)

DefaultAssay(integrated) <- "MAGIC_RNA"
p1 <- FeaturePlot(object = integrated, features = c("VEGFA"), raster = FALSE, reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
p2 <- FeaturePlot(object = integrated, features = c("VEGFA"), raster = FALSE, reduction = "chromharmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
p3 <- FeaturePlot(object = integrated, features = c("VEGFA"), raster = FALSE, reduction = "ATACharmonyumap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)


savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "featureplot/VEGFA_umap_chromumap_ATAC_noraster.pdf", sep = ""), width = 15, height = 5.5)
p1 | p2 | p3
dev.off()

DefaultAssay(integrated) <- "MAGIC_RNA"

rm(plot_list)
plot_list <- list()
genes <- c("HSPH1", "DNAJB1", "HSPA1B", "DNAJB1")
for (i in 1:length(genes)) {
    p1 <- FeaturePlot(object = integrated, features = genes[i], raster = FALSE, reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p1
}

pdf(paste(savedir, "featureplot/HSP_umap_chromumap_ATAC_noraster.pdf", sep = ""), width = 6, height = 5.5)
plot_list
dev.off()

rm(plot_list)
plot_list <- list()
genes <- c("EGR1", "EGR3", "KLF2", "SP1", "SP3", "HIF1A")
for (i in 1:length(genes)) {
    p1 <- FeaturePlot(object = integrated, features = genes[i], raster = FALSE, reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p1
}

pdf(paste(savedir, "featureplot/Recurrent_TF_umap_RNA.pdf", sep = ""), width = 6, height = 5.5)
plot_list
dev.off()


rm(plot_list)
plot_list <- list()
genes <- c("EGR1", "EGR3", "KLF2", "SP1", "SP3", "HIF1A")
for (i in 1:length(genes)) {
    p1 <- FeaturePlot(object = integrated, features = genes[i], raster = FALSE, reduction = "umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p1
}

pdf(paste(savedir, "featureplot/Recurrent_TF_umap_RNA.pdf", sep = ""), width = 6, height = 5.5)
plot_list
dev.off()

genes <- unique(c("VEGFA", "HSPH1", "DNAJB1", "HSPA1B", "EGR1", "EGR3", "KLF2", "SP1", "SP3", "HIF1A"))
p <- DotPlot(integrated, genes, assay = "RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_Genes_RNA_celltypes.pdf"), height = 5.5, width = 9)
print(p)
dev.off()

TFs <- c("MA0079.4", "MA0746.2", "MA1515.1", "MA0162.4", "MA0732.1", "MA1106.1")
## Arranged as SP1,SP2, KLF2, EGR1, EGR3, HIF1A
p <- DotPlot(integrated, TFs, assay = "chromvar", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_TFs_chromvar_celltypes.pdf"), height = 4, width = 9)
print(p)
dev.off()


DefaultAssay(snATAC) <- "chromvar"
snATAC -> integrated
paired_cells <- row.names(integrated@meta.data[grep(c("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931"), integrated@meta.data$Sample), ])
integrated_paired <- subset(integrated, cells = paired_cells)

grep("Endo-Pericytes|Oligo_OPC|Astrocytes|Immune|Neuron", integrated_paired@meta.data$celltypes, invert = TRUE)
cellnames <- rownames(integrated_paired@meta.data[grep("Endo-Pericytes|Oligo_OPC|Astrocytes|Immune|Neuron", integrated_paired@meta.data$celltypes, invert = TRUE), ])
integreated_paired_mal <- subset(integrated_paired, cells = cellnames)

DefaultAssay(integreated_paired_mal) <- "chromvar"
chromvar_marker <- FindMarkers(integreated_paired_mal, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type")

p <- DoHeatmap(integreated_paired_mal, features = head(rownames(chromvar_marker[order(-chromvar_marker$avg_log2FC), ]), 20), group.by = "tumor_type")

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/all_Rec_vs_Prim_tumor_type.pdf")
p
dev.off()

p <- DoHeatmap(integreated_paired_mal, features = head(rownames(chromvar_marker[order(-chromvar_marker$avg_log2FC), ]), 20), group.by = "Sample_pt_tumortype")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/all_Rec_vs_Prim_sample.pdf")
p
dev.off()

p <- DoHeatmap(integreated_paired_mal, features = TFs, group.by = "Sample_pt_tumortype")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/all_Rec_vs_Prim_ZFs.pdf")
p
dev.off()

p <- DoHeatmap(integreated_paired_mal, features = TFs, group.by = "tumor_type")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/all_Rec_vs_Prim_TFs_tumor_type.pdf")
p
dev.off()


### Magic RNA require lot of memory
library(Seurat)
library(reticulate)
library(Rmagic)
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

DefaultAssay(integrated) <- "RNA"
integrated <- magic(integrated, npca = 30) ## imputing the RNA data as for RNA PCs are 20

### Since we will be performing the cell type annotation as well
GCP <- c("PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L", "ATOH1", "GLI2", "GLI1", "PTCH1")
GNearly <- c("NHLH2", "NHLH1", "STMN2", "TUBB3", "NRG1", "GRIN2B", "KIF21B", "DISP3")
GNlate <- c("NEUROD1", "NEUROD2", "OXR1", "NRXN2", "RBFOX3", "BRINP1", "GRIN2C")
Neuronal <- c("TRPM3", "OTX2", "ROBO3", "UNC5D", "ITGA3", "ADRA1A", "OTX2-AS1", "PRDM6")
Stem_cells <- c("NES", "SOX2", "MKI67", "TOP2A")
HSP <- c("HSPH1", "HSPA1B", "HSPA1A", "DNAJB1", "HSPB1", "VEGFA", "HIF1A")
RB <- c("RPL35", "RPS5", "RPL4", "RPL26", "RPL29", "RPS17")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[, 1]
genes <- c(GCP, Stem_cells, HSP, RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

# dir.create(paste0(savedir, "/featureplot/"), showWarnings = FALSE, recursive = TRUE)
# pdf(paste0(savedir, "/featureplot/integrated_marker_genes.pdf"), width = 10, height = 8)
# print(p)
# dev.off()

p <- DotPlot(integrated, genes, assay = "RNA", group.by = "sub_celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/chromvar_Genes_RNA_celltypes.pdf"), height = 15, width = 12)
print(p)
dev.off()

clusters <- levels(integrated@meta.data$seurat_clusters)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
    cellnames <- integrated@meta.data[grep(paste("^", clusters[i], "$", sep = ""), integrated@meta.data$seurat_clusters), ] %>% rownames()
    p <- DimPlot(integrated,
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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "UMAP/cluster_splitted.pdf", sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

scATAC@reductions$ATACharmony <- combined_harmony@reductions$ATACharmony
scATAC@reductions$ATACharmonyumap <- combined_harmony@reductions$ATACharmonyumap

clusters <- unique(scATAC@meta.data$celltypes)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
    cellnames <- scATAC@meta.data[grep(paste("^", clusters[i], "$", sep = ""), scATAC@meta.data$celltypes), ] %>% rownames()
    p <- DimPlot(scATAC,
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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "UMAP/ATAC_celltypes_splitted.pdf", sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

clusters <- unique(scATAC@meta.data$celltypes)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
    cellnames <- scATAC@meta.data[grep(paste("^", clusters[i], "$", sep = ""), scATAC@meta.data$celltypes), ] %>% rownames()
    p <- DimPlot(scATAC,
        cells.highlight = cellnames,
        reduction = "ATACharmonyumap",
        label = FALSE, cols.highlight = "deeppink2",
        sizes.highlight = 0.5,
        cols = "gray92"
    ) +
        ggtitle(paste(clusters[i])) +
        NoLegend() +
        theme(plot.title = element_text(hjust = 0.2))
    plot_list[[i]] <- p
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "UMAP/ATACharmonyumap_celltypes_splitted.pdf", sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

clusters <- unique(integrated@meta.data$predicted_celltypes)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
    cellnames <- integrated@meta.data[grep(paste("^", clusters[i], "$", sep = ""), integrated@meta.data$predicted_celltypes), ] %>% rownames()
    p <- DimPlot(integrated,
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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "UMAP/celltype_splitted.pdf", sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

clusters <- unique(integrated@meta.data$predicted.id)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
    cellnames <- integrated@meta.data[grep(paste("^", clusters[i], "$", sep = ""), integrated@meta.data$predicted.id), ] %>% rownames()
    p <- DimPlot(integrated,
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

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "UMAP/sub_celltype_splitted.pdf", sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

### Add Module Score
library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(integrated) <- "MAGIC_RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/resources/gene_list/SHH_Gold_Taylor_2024_Nat_comm/", pattern = ".txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(integrated)[match(Tcellsubset, rownames(integrated), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

integrated <- AddModuleScore(integrated, features = markers, slot = "data")

index <- grep("Cluster1$", colnames(integrated@meta.data))
colnames(integrated@meta.data)[index:(index + 34)] <- filename

rm(plot_list)
plot_list <- list()
for (j in 1:length(filename)) {
    plot_list[[j]] <- FeaturePlot(integrated,
        features = filename[j],
        reduction = "umap", raster = TRUE
    ) +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra)
}


dir.create(paste(savedir, "featureplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "/featureplot/SHH_MB_het_gene_score_combined.pdf", sep = ""))
print(plot_list)
dev.off()

res <- c(1, 1.2, 1.4, 1.6, 1.8, 2, 2.2)
for (j in 2:length(res)) {
    DefaultAssay(integrated) <- "ATAC"
    integrated <- FindClusters(integrated, verbose = FALSE, resolution = res[j])
    pdf(paste0(savedir, "UMAP/scRNA_integration_res_", res[j], ".pdf"))
    print(DimPlot(integrated, reduction = "umap", label = TRUE))
    dev.off()

    p <- DotPlot(scATAC, genes, assay = "RNA", group.by = "seurat_clusters2") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/integrated_marker_genes_RNA_divided_33_", res[j], ".pdf"), height = 32, width = 19)
    print(p)
    dev.off()

    p <- DotPlot(integrated, genes, assay = "MAGIC_RNA") +
        scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
        coord_flip() +
        # scale_size(range = c(1, 10)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "/dotplot/integrated_marker_genes_MAGIC_RNA_", res[j], ".pdf"), height = 32, width = 19)
    print(p)
    dev.off()

    clusters <- unique(scATAC@meta.data$seurat_clusters2)
    clusters <- clusters[order(clusters)]
    rm(plot_list)
    plot_list <- list()

    # clusters <- c("33_GCP", "33_GN")
    for (i in 1:length(clusters)) {
        cellnames <- scATAC@meta.data[grep(paste("^", clusters[i], "$", sep = ""), scATAC@meta.data$seurat_clusters2), ] %>% rownames()
        p <- DimPlot(scATAC,
            cells.highlight = cellnames,
            reduction = "ATACharmonyumap",
            label = FALSE, cols.highlight = "deeppink2",
            sizes.highlight = 0.5,
            cols = "gray92"
        ) +
            ggtitle(paste(clusters[i])) +
            NoLegend() +
            theme(plot.title = element_text(hjust = 0.2))
        plot_list[[i]] <- p
    }

    pdf(paste(savedir, "UMAP/cluster_splitted_res_ATAC_harmony_", res[j], ".pdf", sep = ""), width = 5.5, height = 5)
    print(plot_list)
    dev.off()
}

### Adding the celltypes
clus33 <- rownames(scATAC@meta.data[grep("33", scATAC@meta.data$seurat_clusters), ])
umap_embed <- as.data.frame(scATAC@reductions$umap@cell.embeddings)
umap_embed_clus33 <- umap_embed[match(clus33, rownames(umap_embed)), ]
umap_embed_clus33_gt_0 <- rownames(umap_embed_clus33[umap_embed_clus33$umap_2 > 0, ])
umap_embed_clus33_ls_0 <- rownames(umap_embed_clus33[umap_embed_clus33$umap_2 < 0, ])

scATAC@meta.data$seurat_clusters2 <- scATAC@meta.data$seurat_clusters
scATAC@meta.data$seurat_clusters2 <- as.character(scATAC@meta.data$seurat_clusters2)
scATAC@meta.data[match(umap_embed_clus33_gt_0, rownames(scATAC@meta.data)), "seurat_clusters2"] <- "33_GCP"
scATAC@meta.data[match(umap_embed_clus33_ls_0, rownames(scATAC@meta.data)), "seurat_clusters2"] <- "33_GN"

clus_CT <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/clus_celltypes.txt", header = TRUE, sep = "\t")
clus_CT$Celltypes <- gsub("GN premigratory", "GN Premigratory", clus_CT$Celltypes) %>%
    gsub("GN mgratory", "GN migratory", .) %>%
    gsub("GCP-GN ", "GCP-GN", .)


patterns <- clus_CT$Cluster
replacements <- clus_CT$Celltypes
# names(replacement) <- patterns

integrated@meta.data$celltypes <- integrated@meta.data$seurat_clusters

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    integrated@meta.data$celltypes <- str_replace_all(integrated@meta.data$celltypes, pattern, replacements[i])
}

integrated@meta.data$celltypes <- gsub("GN premigratory", "GN Premigratory", integrated@meta.data$celltypes) %>%
    gsub("GN mgratory", "GN migratory", .) %>%
    gsub("GCP-GN ", "GCP-GN", .)

p <- DimPlot(scATAC,
    group.by = "celltypes",
    reduction = "chromumap",
    label = TRUE,
    cols = c(
        "GCP cycling" = "#1a69a1", "GN cycling" = "#436957", "GN Premigratory" = "#5fe3a7", "GN HSP" = "#216647", "GCP HSP" = "#3b95d4",
        "GCP-GN" = "#5bcbcb",
        "GN migratory" = "#30b377", "GN postmigratory" = "#107044", "GCP" = "#6fb9ed", "GCP NES" = "#4593cb", "GCP Ribosomal" = "#155e92",
        "Pericyte" = "#aa40fc",
        "Astrocytes" = "grey", "Immune" = "black", "Oligo_OPC" = "hotpink",
        "Endo-Pericytes" = "#17becf", "Neuron" = "gray33"
    )
)

pdf(paste0(savedir, "UMAP/chromvar_Celltypes_2.pdf"), width = 7, height = 5.5)
p
dev.off()

p <- DimPlot(scATAC,
    group.by = "Sample",
    reduction = "umap",
    label = TRUE
)

pdf(paste0(savedir, "UMAP/integrated_samples.pdf"), width = 7, height = 5.5)
p
dev.off()

p <- DimPlot(scATAC,
    # group.by = "Sample",
    reduction = "ATACharmonyumap",
    label = TRUE
)

pdf(paste0(savedir, "UMAP/integrated_harmony_clustering.pdf"), width = 7, height = 5.5)
p
dev.off()


### Since we will be performing the cell type annotation as well
GCP <- c("PCNT", "BOC", "NTRK3", "SHROOM2", "SFRP1", "PTCH2", "TCERG1L", "ATOH1", "GLI2", "GLI1", "PTCH1")
GNearly <- c("NHLH2", "NHLH1", "STMN2", "TUBB3", "NRG1", "GRIN2B", "KIF21B", "DISP3")
GNlate <- c("NEUROD1", "NEUROD2", "OXR1", "NRXN2", "RBFOX3", "BRINP1", "GRIN2C")
Neuronal <- c("TRPM3", "OTX2", "ROBO3", "UNC5D", "ITGA3", "ADRA1A", "OTX2-AS1", "PRDM6")
Stem_cells <- c("NES", "SOX2", "MKI67", "TOP2A")
HSP <- c("HSPH1", "HSPA1B", "HSPA1A", "DNAJB1", "HSPB1", "VEGFA", "HIF1A")
RB <- c("RPL35", "RPS5", "RPL4", "RPL26", "RPL29", "RPS17")
nonmalignant_cells <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/non_mal_genes")[, 1]
genes <- c(GCP, Stem_cells, HSP, RB, GNearly, GNlate, Neuronal, nonmalignant_cells)

# dir.create(paste0(savedir, "/featureplot/"), showWarnings = FALSE, recursive = TRUE)
# pdf(paste0(savedir, "/featureplot/integrated_marker_genes.pdf"), width = 10, height = 8)
# print(p)
# dev.off()

p <- DotPlot(integrated, genes, assay = "RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 8)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/integrated_marker_genes_celltypes.pdf"), height = 30, width = 8)
print(p)
dev.off()

p <- DotPlot(integrated, genes, assay = "MAGIC_RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    # scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/integrated_marker_genes_celltypes_MAGIC_RNA.pdf"), height = 30, width = 8)
print(p)
dev.off()

p <- DotPlot(integrated, nonmalignant_cells, assay = "RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 8)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/integrated_marker_genes_non_malignant.pdf"), height = 18, width = 8)
print(p)
dev.off()

DefaultAssay(integrated) <- "RNA"
RNA_markers <- FindAllMarkers(integrated)
write.table(RNA_markers, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/RNA_markers_res2.2.txt", quote = F, row.names = T, col.names = T, sep = "\t")

saveRDS(integrated, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_2.RDS")

DefaultAssay(integrated) <- "ATAC"
ATAC_markers <- FindAllMarkers(integrated)
write.table(ATAC_markers, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/ATAC_markers_res2.2.txt", quote = F, row.names = T, col.names = T, sep = "\t")

# Find for the recurrent Peaks
# get top differentially accessible peaks
ATAC_clus <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/ATAC_markers_res2.2.txt")
dir.create(paste0(savedir, "Table/ATAC_markers"), showWarnings = FALSE, recursive = TRUE)
clus <- unique(ATAC_clus$cluster)

for (i in 7:length(clus)) {
    try({
        ATAC_clus_subset <- ATAC_clus[grep(paste0("^", clus[i], "$"), ATAC_clus$cluster), ]
        da_peaks <- ATAC_clus_subset[ATAC_clus_subset$avg_log2FC > 1 & ATAC_clus_subset$p_val_adj < 0.05, "gene"]
        top.da.peak <- grep("^GL|^KI", da_peaks, invert = TRUE, value = TRUE)
        # test enrichment

        enriched.motifs <- FindMotifs(
            object = integrated,
            assay = "ATAC",
            features = top.da.peak
        )
        assign(paste0("ATAC_clus_", clus[i]), enriched.motifs)
        write.table(enriched.motifs, paste0(savedir, "Table/ATAC_markers/nocutoff/ATAC_clus_", clus[i], "_nocutoff.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

        enriched.motifs_sig <- enriched.motifs[enriched.motifs$p.adjust < 0.05 & enriched.motifs$fold.enrichment > 1, ]
        assign(paste0("ATAC_clus_sig_", clus[i]), enriched.motifs_sig)
        write.table(enriched.motifs_sig, paste0(savedir, "Table/ATAC_markers/lfc_1_sig_motif/ATAC_clus_", clus[i], "_lFC_1_sig.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
    })
}

genes <- c(
    "PDGFRB", "RGS5", "ACTA2", "CDH6", "COL1A1", "COL1A2", "COL6A2", "COL3A1", "EDNRA", "DCN", "BGN", "NR2F2", "PRKG1",
    "FLT1", "PTPRB", "EGFL7", "CLDN5", "VWF", "CDH5", "CD93", "TIE1", "ERG", "ABCB1", "ENG", "EGR1", "CD34", "PLVAP"
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
p <- DotPlot(integrated, genes, assay = "MAGIC_RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    # scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/integrated_marker_genes_pericyte_endo.pdf"), height = 10, width = 19)
print(p)
dev.off()

p <- DotPlot(integrated, genes, assay = "RNA") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/integrated_marker_genes_pericyte_endo_RNA.pdf"), height = 10, width = 19)
print(p)
dev.off()

### Cell type annotation
# Integrating with scRNA-seq data
# To help interpret the scATAC-seq data, we can classify cells based on an scRNA-seq experiment from the same biological system (human PBMC). We utilize methods for
# cross-modality integration and label transfer, described https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub, with a more in-depth
# tutorial here. We aim to identify shared correlation patterns in the gene activity matrix and scRNA-seq dataset to identify matched biological states across
# the two modalities. This procedure returns a classification score for each cell for each scRNA-seq-defined cluster label.

# Here we load a pre-processed scRNA-seq dataset for human PBMCs, also provided by 10x Genomics. You can download the raw data for this experiment from the 10x website,
# and view the code used to construct this object on GitHub. Alternatively, you can download the pre-processed Seurat object here.

### Loading the snRNA seq data
## Since it is very memory intensive so ran using
# sbatch /diazlab/data3/.abhinav/projects/SHH/scATAC/merging/ATAC_projection.sh

### Motif Analysis
# In this tutorial, we will perform DNA sequence motif analysis in Signac. We will explore two complementary options for performing motif analysis:
# one by finding overrepresented motifs in a set of differentially accessible peaks, one method performing differential motif activity analysis between groups of cells.
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

# add motif information
DefaultAssay(integrated) <- "ATAC"
integrated <- AddMotifs(
    object = integrated,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
)

# To facilitate motif analysis in Signac, we have created the Motif class to store all the required information, including a list of position weight matrices
# (PWMs) or position frequency matrices (PFMs) and a motif occurrence matrix. Here, the AddMotifs() function constructs a Motif object and adds it to our mouse
# brain dataset, along with other information such as the base composition of each peak. A motif object can be added to any Seurat assay using the SetAssayData()
# function. See the object interaction vignette for more information.

### Performing Differential in scATAC for RNA
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

### Patient ID
patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns)
patterns <- gsub("-", "_", patterns)
replacements <- long_sample$Patient.ID
# names(replacement) <- patterns

integrated@meta.data$Sample <- gsub("sample_", "", integrated@meta.data$dataset) %>% gsub("_128.*.", "", .)
integrated@meta.data$Patient.ID <- integrated@meta.data$Sample

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    integrated@meta.data$Patient.ID <- str_replace_all(integrated@meta.data$Patient.ID, pattern, replacements[i])
}

replacements <- long_sample$Tumor.type
integrated@meta.data$tumor_type <- integrated@meta.data$Sample

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    integrated@meta.data$tumor_type <- str_replace_all(integrated@meta.data$tumor_type, pattern, replacements[i])
}

# Finding overrepresented motifs
# To identify potentially important cell-type-specific regulatory sequences, we can search for DNA motifs that are overrepresented in a set of peaks that are
# differentially accessible between cell types.
# Here, we find differentially accessible peaks between Pvalb and Sst inhibitory interneurons. For sparse data (such as scATAC-seq), we find it is often
# necessary to lower the min.pct threshold in FindMarkers() from the default (0.1, which was designed for scRNA-seq data).
# We then perform a hypergeometric test to test the probability of observing the motif at the given frequency by chance, comparing with a background set of
# peaks matched for GC content.
DefaultAssay(integrated) <- "ATAC"
integrated@meta.data$Patient_tumortype <- paste(integrated@meta.data$Patient.ID, integrated@meta.data$tumor_type, sep = "_")

pdf(paste0(savedir, "UMAP/integrated_patient_tumortype_splitted.pdf"), width = 15, height = 12)
DimPlot(integrated, group.by = "Patient_tumortype", split.by = "Patient_tumortype", ncol = 4, reduction = "umap") + NoLegend()
dev.off()

paired_samples <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"

for (i in 1:length(paired_samples)) {
    da_peaks_sample <- FindMarkers(
        object = integrated,
        ident.1 = paste0(paired_samples[i], "_Recurrent"),
        ident.2 = paste0(paired_samples[i], "_Primary"),
        only.pos = FALSE,
        test.use = "LR",
        min.pct = 0.05,
        group.by = "Patient_tumortype",
        latent.vars = "nCount_ATAC"
    )
    assign(paste0(paired_samples[i], "_rec_vs_pri"), da_peaks_sample)
    write.table(da_peaks_sample, paste0(savedir, "scdifferential/", paired_samples[i], "_rec_vs_pri_ATAC_both.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
}

### Finding the closest genes to the differential peaks
filenames <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential", pattern = "_rec_vs_pri_ATAC_both.txt", full.names = TRUE)
samplename <- gsub("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/", "", filenames) %>% gsub("_rec_vs_pri_ATAC_both.txt", "", .)

DefaultAssay(integrated) <- "ATAC"
for (i in 1:length(samplename)) {
    try({
        da_peaks <- read.table(filenames[i])
        # da_peaks_rec <- rownames(da_peaks[da_peaks$avg_log2FC > 0.1 & da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
        # da_peaks_pri <- rownames(da_peaks[da_peaks$avg_log2FC < -0.1 & da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2, ])

        da_peaks_rec <- rownames(da_peaks[da_peaks$avg_log2FC > 0.1, ])
        da_peaks_pri <- rownames(da_peaks[da_peaks$avg_log2FC < -0.1, ])


        closest_genes_rec <- ClosestFeature(integrated, regions = da_peaks_rec)
        closest_genes_prim <- ClosestFeature(integrated, regions = da_peaks_pri)

        write.table(closest_genes_rec, paste0("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/close_features/", samplename[i], "_close_gene_rec_lfc_0.1.txt"),
            quote = F, row.names = T, col.names = T, sep = "\t"
        )
        write.table(closest_genes_prim, paste0("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/close_features/", samplename[i], "_close_gene_prim_lfc_0.1.txt"),
            quote = F, row.names = T, col.names = T, sep = "\t"
        )
    })
}

### Performing the Peak Annotation for the differential Genes
#### Identifying the common peaks between the samples
library(UpSetR)
library(dplyr)

filename <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/",
    pattern = "_rec_vs_pri_ATAC_both.txt",
    recursive = TRUE, full.names = TRUE
)

markersname <- gsub("_rec_vs_pri_ATAC_both.txt", "", basename(filename))

# extracting the primary and recurrent markers
for (i in 1:length(markersname)) {
    PTmarkers <- read.table(filename[i], sep = "\t", header = TRUE)
    recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val < 0.005 & PTmarkers$avg_log2FC > 0, ]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
    assign(paste(markersname[i], "_peaks_recurrent", sep = ""), recurrentmarkers)
    primarymarkers <- rownames(PTmarkers[PTmarkers$p_val < 0.005 & PTmarkers$avg_log2FC < 0, ])
    assign(paste(markersname[i], "_peaks_primary", sep = ""), primarymarkers)
}

recurrentname <- ls(pattern = "_peaks_recurrent")

H <- list(
    "PT_7WYPEC3Q_peaks_recurrent" = PT_7WYPEC3Q_peaks_recurrent,
    "PT_9S6WMQ92_peaks_recurrent" = PT_9S6WMQ92_peaks_recurrent,
    "PT_H3WWDMW9_peaks_recurrent" = PT_H3WWDMW9_peaks_recurrent,
    "PT_XA98HG1C_peaks_recurrent" = PT_XA98HG1C_peaks_recurrent
)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H, show_elements = FALSE, stroke_color = "Black", stroke_linetype = "solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/"
dir.create(paste(savedir, "Table/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "Table/Peaks_recurrent_high_venny_lfc_0_p_0.005.pdf", sep = ""))
a
dev.off()

### Making an upset plot
n <- max(
    length(PT_7WYPEC3Q_peaks_recurrent), length(PT_9S6WMQ92_peaks_recurrent),
    length(PT_H3WWDMW9_peaks_recurrent), length(PT_XA98HG1C_peaks_recurrent)
)

length(PT_7WYPEC3Q_peaks_recurrent) <- n
length(PT_9S6WMQ92_peaks_recurrent) <- n
length(PT_H3WWDMW9_peaks_recurrent) <- n
length(PT_XA98HG1C_peaks_recurrent) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q_peaks_recurrent), (PT_9S6WMQ92_peaks_recurrent),
    (PT_H3WWDMW9_peaks_recurrent), (PT_H3WWDMW9_peaks_recurrent)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")

### To find the interseected genes
list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q_peaks_recurrent,
    "PT_9S6WMQ92" = PT_9S6WMQ92_peaks_recurrent,
    "PT_H3WWDMW9" = PT_H3WWDMW9_peaks_recurrent,
    "PT_XA98HG1C" = PT_H3WWDMW9_peaks_recurrent
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
library(stringr)
final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

df_int$int2 <- gsub("\\|", "__", df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3), 1]
morethan_3 <- df_int[(df_int$sample_num >= 3), 1]

write.table(df_int, paste(savedir, "Table/recurrent_peaks_intersect_lfc_0_pval_0.005.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "Table/recurrent_intersected_peak_number_lfc_0_pval_0.005.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir, "Table/peak_recurrent_all_samples_lfc_0_pval_0.005.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir, "Table/peak_recurrent_morethan_3_lfc_0_pval_0.005.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

DefaultAssay(integrated) <- "ATAC"
da_peaks <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/Table/peak_recurrent_morethan_3_lfc_0.txt")[, 1]

closest_genes_rec <- ClosestFeature(integrated, regions = da_peaks)
# closest_genes_prim <- ClosestFeature(integrated, regions = da_peaks_pri)

write.table(closest_genes_rec,
    paste0("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/close_features/Close_gene_rec_lfc_0_pvalue_0.005.txt"),
    quote = F, row.names = F, col.names = T, sep = "\t"
)


# Find the motifs for the recurrent Peaks
# get top differentially accessible peaks
dir.create(paste0(savedir, "motif/"))
scdiff_peaks <- ls(pattern = "_rec_vs_pri")
sample <- gsub("_rec_vs_pri", "", scdiff_peaks)
for (i in 1:length(scdiff_peaks)) {
    da_peaks <- get(scdiff_peaks[i])
    # since we are interested in primary as we already did for the recurrent for recurrent I have used p_val < 0.005
    top.da.peak <- row.names(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2 & da_peaks$avg_log2FC < 0, ])
    top.da.peak <- grep("^GL|^KI", top.da.peak, invert = TRUE, value = TRUE)
    # test enrichment
    enriched.motifs <- FindMotifs(
        object = integrated,
        features = top.da.peak
    )
    assign(paste0(sample[i], "_primary_high_motif"), enriched.motifs)
    write.table(enriched.motifs, paste0(savedir, "motif/", paired_samples[i], "_primary_high_ATAC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
}


sig_motif <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/", pattern = "_primary_high_ATAC.txt", full.names = TRUE)
samplename <- gsub("_primary_high_ATAC.txt", "", basename(sig_motif))
for (i in 1:length(sig_motif)) {
    TFs <- read.table(sig_motif[i], header = TRUE)
    TFs_sig <- TFs[TFs$p.adjust < 0.05, ]
    assign(paste0(samplename[i], "_motif_sig"), TFs_sig)
}


### Performing upset plot on these samples
#### Making an upset plot for the GCP markers
### Making an Upset plot ####
motif_sig <- ls(pattern = "_motif_sig")
PT_H3WWDMW9_motif_sig <- ""

library(data.table)
PT_7WYPEC3Q <- get(motif_sig[1])["motif.name"][, 1]
PT_9S6WMQ92 <- get(motif_sig[2])["motif.name"][, 1]
# PT_H3WWDMW9 <- get(motif_sig[3])["motif.name"][, 1]
PT_XA98HG1C <- get(motif_sig[4])["motif.name"][, 1]
PT_H3WWDMW9 <- "" # since in primary does not have any motif

list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q, "PT_9S6WMQ92" = PT_9S6WMQ92,
    "PT_H3WWDMW9" = PT_H3WWDMW9, "PT_XA98HG1C" = PT_XA98HG1C
)

n <- max(
    length(PT_7WYPEC3Q), length(PT_9S6WMQ92), length(PT_H3WWDMW9), length(PT_XA98HG1C)
)

length(PT_7WYPEC3Q) <- n
length(PT_9S6WMQ92) <- n
length(PT_H3WWDMW9) <- n
length(PT_XA98HG1C) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q), (PT_9S6WMQ92), (PT_H3WWDMW9), (PT_XA98HG1C)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_motif_sig", "", motif_sig)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "motif/Primary_motif_enriched_significant.pdf", sep = ""), width = 12, height = 7)
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
library(stringr)
df_int$sample_num <- str_count(df_int$int, "PT_")
motif_filter_common <- df_int[(df_int$sample_num > 1), 1]
motif_filter_all <- df_int[(df_int$sample_num > 2), 1]

write.table(df_int, paste(savedir, "motif/primary_motif_intersect_name.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "motif/primary_motif_intersect_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(motif_filter_common, paste(savedir, "motif/primary_motif_intersect_markers_morethan2.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
write.table(motif_filter_all, paste(savedir, "motif/primary_motif_intersect_markers_all.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)


### VEGFA visualization

# find DA peaks overlapping gene of interest
# regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(pbmc, "CD4"))

paired_cells <- row.names(integrated@meta.data[grep(c("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931"), integrated@meta.data$Sample), ])
integrated_paired <- subset(integrated, cells = paired_cells)

### Malignant cells
grep("Endo-Pericytes|Oligo_OPC|Astrocytes|Immune|Neuron", integrated_paired@meta.data$celltypes, invert = TRUE)
cellnames <- rownames(integrated_paired@meta.data[grep("Endo-Pericytes|Oligo_OPC|Astrocytes|Immune|Neuron", integrated_paired@meta.data$celltypes, invert = TRUE), ])
integreated_paired_mal <- subset(integrated_paired, cells = cellnames)

cov_plot <- CoveragePlot(
    object = integrated,
    region = "ANGPT2",
    # extend.upstream = 500000,
    # extend.downstream = 500000,
    group.by = "Patient_tumortype"
)

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/coverage_plot/ANGPT2_link_peaks_integrated.pdf")
cov_plot
dev.off()

expr_plot <- ExpressionPlot(
    object = integrated_paired,
    features = "VEGFA",
    assay = "MAGIC_RNA",
    group.by = "Patient_tumortype"
) + geom_boxplot()

p <- CombineTracks(
    plotlist = list(cov_plot),
    expression.plot = expr_plot,
)

dir.create(paste0(savedir, "coverage_plot"), showWarnings = FALSE)
pdf(paste0(savedir, "coverage_plot/VEGFA_tumortype_patient_paired_2.pdf"))
p
dev.off()

pdf(paste0(savedir, "coverage_plot/VEGFA_tumortype_patient_paired_Expr_plot.pdf"))
expr_plot
dev.off()


### Removed unpaired and non-malignant cells
genes <- unique(c("VEGFA", "HSPH1", "DNAJB1", "HSPA1B", "EGR1", "EGR3", "KLF2", "SP1", "SP3", "HIF1A", "HSF1"))
p <- DotPlot(integreated_paired_mal, genes, assay = "RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_Genes_RNA_celltypes.pdf"), height = 5.5, width = 9)
print(p)
dev.off()

genomic_region <- "chr6-43769103-43770414"
p <- DotPlot(integreated_paired_mal, genomic_region, assay = "ATAC", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_Genes_ATAC_VEGFA_celltypes.pdf"), height = 5.5, width = 9)
print(p)
dev.off()


TFs <- c("MA0079.4", "MA0746.2", "MA1515.1", "MA0162.4", "MA0732.1", "MA1106.1", "MA0486.2")
## Arranged as SP1,SP2, KLF2, EGR1, EGR3, HIF1A
p <- DotPlot(integreated_paired_mal, TFs, assay = "chromvar", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_TFs_chromvar_celltypes.pdf"), height = 4, width = 9)
print(p)
dev.off()

#### Testing using the deviation score for the plotting
integreated_paired_mal2 <- integreated_paired_mal
integreated_paired_mal2@assays$chromvar@scale.data <- integreated_paired_mal2@assays$chromvar@data

genes <- unique(c("VEGFA", "HSPH1", "DNAJB1", "HSPA1B", "EGR1", "EGR3", "KLF2", "SP1", "SP3", "HIF1A", "HSF1"))
p <- DotPlot(integreated_paired_mal2, genes, assay = "RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_Genes_RNA_celltypes_data_slot.pdf"), height = 5.5, width = 9)
print(p)
dev.off()

TFs <- c("MA0079.4", "MA0746.2", "MA1515.1", "MA0162.4", "MA0732.1", "MA1106.1", "MA0486.2")
## Arranged as SP1,SP2, KLF2, EGR1, EGR3, HIF1A
p <- DotPlot(integreated_paired_mal2, TFs, assay = "chromvar", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_TFs_chromvar_celltypes_data_slot.pdf"), height = 4, width = 9)
print(p)
dev.off()

### Removed 2118 sample since it is skewing the whole analysis
cells_no2118 <- rownames(integreated_paired_mal2@meta.data[grep("7316_2118", integreated_paired_mal2@meta.data$Sample, invert = TRUE), ])
integreated_paired_mal2_no2118 <- subset(integreated_paired_mal2, cells = cells_no2118)

genes <- unique(c("VEGFA", "HSPH1", "DNAJB1", "HSPA1B", "EGR1", "EGR3", "KLF2", "SP1", "SP3", "HIF1A", "HSF1"))
p <- DotPlot(integreated_paired_mal2_no2118, genes, assay = "RNA", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_Genes_RNA_celltypes_data_slot_no2118.pdf"), height = 5.5, width = 9)
print(p)
dev.off()

TFs <- c("MA0079.4", "MA0746.2", "MA1515.1", "MA0162.4", "MA0732.1", "MA1106.1", "MA0486.2")
## Arranged as SP1,SP2, KLF2, EGR1, EGR3, HIF1A
p <- DotPlot(integreated_paired_mal2_no2118, TFs, assay = "chromvar", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/Interested_paired_mal_TFs_chromvar_celltypes_data_slot_no2118.pdf"), height = 4, width = 9)
print(p)
dev.off()


# Computing motif activities
# We can also compute a per-cell motif activity score by running chromVAR. This allows us to visualize motif activities per cell, and also provides an alternative method
#  of identifying differentially-active motifs between cell types.

# ChromVAR identifies motifs associated with variability in chromatin accessibility between cells. See the chromVAR paper for a complete description of the method.
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(integrated))) %in% main.chroms)
integrated[["ATAC"]] <- subset(integrated[["ATAC"]], features = rownames(integrated[["ATAC"]])[keep.peaks])

integrated <- RunChromVAR(
    object = integrated,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(integrated) <- "chromvar"

# look at the activity of Mef2c
rm(plot_list)
plot_list <- list()
motifs <- c("MA1106.1", "MA0079.4", "MA0746.2", "MA0685.1", "MA0747.1", "MA1564.1", "MA0162.4", "MA0472.2", "MA0732.1", "MA0733.1")
motifs <- c("MA0848.1", "MA0849.1", "MA0850.1", "MA1515.1", "MA0162.4", "MA0732.1", "MA0079.4", "MA0486.2")
motifs <- unique(c("MA0848.1", "MA0849.1", "MA0157.2", "MA0850.1", "MA1515.1", "MA0162.4", "MA0732.1", "MA0079.4", "MA0486.2", "MA1106.1", "MA0079.4", "MA0746.2", "MA0685.1", "MA0747.1", "MA1564.1", "MA0162.4", "MA0472.2", "MA0732.1", "MA0733.1"))
for (i in 1:length(motifs)) {
    p <- FeaturePlot(
        object = integrated,
        features = motifs[i],
        min.cutoff = "q10",
        max.cutoff = "q90",
        raster = TRUE
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}

pdf(paste0(savedir, "motif/FOXO4_FOXO6_FOXP3_SP1_KLF2_EGR1_EGR3_SP1.pdf"))
plot_list
dev.off()

rm(plot_list)
plot_list <- list()
motifs <- c("MA1106.1", "MA0079.4", "MA0746.2", "MA0162.4")
for (i in 1:length(motifs)) {
    p <- FeaturePlot(
        object = integrated,
        features = motifs[i],
        # min.cutoff = "q10",
        # max.cutoff = "q90",
        raster = TRUE,
        slot = "data"
    ) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
    plot_list[[i]] <- p
}


pdf(paste0(savedir, "motif/HIF1A_SP1_SP3_EGR1_nocutoff.pdf"))
plot_list
dev.off()

### Making an Motif plot
p <- MotifPlot(
    object = integrated,
    motifs = motifs,
    assay = "ATAC"
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste0(savedir, "motif/motifplot.pdf"), width = 12, height = 10)
p
dev.off()

chromvar_markers <- FindAllMarkers(integrated)
write.table(chromvar_markers, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/chromvar_markers.txt", sep = "\t", row.names = T, col.names = T, quote = F)

### Adding the transcription factors gene into it.
chromvar_markers <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/chromvar_markers.txt", header = TRUE)
jaspar_id_name <- read.table("/diazlab/data3/.abhinav/resources/pfms/JASPAR2020_CORE_non-redundant_pfms_jaspar/final_JASPAR_MAid_TFnames.txt")
colnames(jaspar_id_name) <- c("motif_id", "motif_name")

ls patterns <- jaspar_id_name$motif_id
replacements <- jaspar_id_name$motif_name

chromvar_markers$TF_name <- rownames(chromvar_markers)

library(stringr)
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    chromvar_markers$TF_name <- str_replace_all(chromvar_markers$TF_name, pattern, replacements[i])
}

write.table(
    chromvar_markers,
    "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/chromvar_markers_with_TFnames.txt",
    quote = F, row.names = T, col.names = T, sep = "\t"
)

### Finding all the TFs based on the celltypes
Idents(integreated_paired_mal) <- integreated_paired_mal@meta.data$celltypes
chromvar_markers_CT <- FindAllMarkers(integreated_paired_mal, assay = "chromvar")
write.table(chromvar_markers_CT, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/chromvar_markers_paired_mal_celltypes.txt", sep = "\t", row.names = T, col.names = T, quote = F)

### Adding the transcription factors gene into it.
# chromvar_markers <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/chromvar_markers.txt", header = TRUE)
jaspar_id_name <- read.table("/diazlab/data3/.abhinav/resources/pfms/JASPAR2020_CORE_non-redundant_pfms_jaspar/final_JASPAR_MAid_TFnames.txt")
colnames(jaspar_id_name) <- c("motif_id", "motif_name")

patterns <- jaspar_id_name$motif_id
replacements <- jaspar_id_name$motif_name

chromvar_markers_CT$TF_name <- rownames(chromvar_markers_CT)

library(stringr)
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    chromvar_markers_CT$TF_name <- str_replace_all(chromvar_markers_CT$TF_name, pattern, replacements[i])
}

write.table(
    chromvar_markers_CT,
    "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/chromvar_markers_CT_with_TFnames.txt",
    quote = F, row.names = T, col.names = T, sep = "\t"
)

### Identifying the quantitative changes
scATAC@meta.data$Sample_pt_tumortype <- paste(scATAC@meta.data$Patient_tumortype, scATAC@meta.data$Sample, sep = "_")
data2_melted <- as.data.frame(table(scATAC@meta.data$celltypes, scATAC@meta.data$Sample_pt_tumortype))
colnames(data2_melted) <- c("celltype", "samples", "percentage")
library(ggplot2)
# create a dataset
# Stacked + percent

p <- ggplot(data2_melted, aes(fill = celltype, y = percentage, x = samples)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = c(
        "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc",
        "#8c564b", "#e377c2", "#b5bd61", "#17becf",
        "#aec7e8", "#ffbb78", "blanchedalmond", "darkseagreen3", "plum2",
        "yellow", "hotpink", "black", "blue", "white"
    )) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste0(savedir, "Table/sample_celltype_percentage.pdf"))
p
dev.off()

DefaultAssay(integrated) <- "ATAC"
integrated <- FindClusters(integrated, verbose = FALSE, resolution = 2.2)

DefaultAssay(integrated) <- "chromvar"
p <- DotPlot(integrated, motifs, assay = "chromvar") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "dotplot/FOXO4_FOXO6_FOXP3_SP1_KLF2_EGR1_EGR3_SP1_HSF1_clus.pdf"), height = 5, width = 26)
print(p)
dev.off()

### Since cluster 42 has been skewing the all the Transcription factor removing that
cells_no_42 <- rownames(integrated@meta.data[grep("42", integrated@meta.data$seurat_clusters, invert = TRUE), ])
integrated_no_42 <- subset(integrated, cells = cells_no_42)

# DefaultAssay(integrated_no_42) <- "chromvar"
# p <- DotPlot(integrated_no_42, motifs, assay = "chromvar", group.by = "celltypes") +
#     scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
#     coord_flip() +
#     scale_size(range = c(1, 10)) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# dir.create(paste0(savedir, "motif/dotplot/"), showWarnings = FALSE, recursive = TRUE)
# pdf(paste0(savedir, "/motif/dotplot/HIF1A_SP1_SP3_EGR1_dotplot_noclus_42.pdf"), height = 5, width = 15)
# print(p)
# dev.off()
# integrated <- FindVariableFeatures(integrated, assay = "chromvar", selection.method = "dispersion")
### FindVariable features is not working taking all the motifs whicha re differential between primary and recurrent
diff_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/scdifferential",
    pattern = "rec_vs_pri_motif_both.txt",
    full.names = TRUE
)

samplename <- gsub("_rec_vs_pri_motif_both.txt", "", basename(diff_files))
for (i in 1:length(diff_files)) {
    motif_diff <- read.table(diff_files[i])
    motif_diff_sig <- motif_diff[motif_diff$p_val_adj < 0.05 & (motif_diff$avg_diff > 1 | motif_diff$avg_diff < -0.5), ]
    assign(paste0(samplename[i], "_motif_sig_3"), rownames(motif_diff_sig))
}

integrated <- RunPCA(
    integrated,
    assay = "chromvar",
    features = unique(c(PT_7WYPEC3Q_motif_sig_3, PT_9S6WMQ92_motif_sig_3, PT_H3WWDMW9_motif_sig_3, PT_XA98HG1C_motif_sig_3)),
    reduction.name = "chrompca",
    reduction.key = "chromPC_"
)

p <- ElbowPlot(integrated, ndims = 50, reduction = "chrompca")

dir.create(paste0(savedir, "motif/elbow_plot/", sep = ""), showWarnings = FALSE)
pdf(paste0(savedir, "motif/elbow_plot/chrompca.pdf"))
p
dev.off()

integrated <- RunUMAP(
    integrated,
    assay = "chromvar",
    reduction.key = "chromUMAP_",
    reduction.name = "chromumap",
    reduction = "chrompca",
    dims = 1:30
)

p <- DimPlot(integrated, reduction = "chromumap", group.by = "celltypes")

dir.create(paste0(savedir, "motif/UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/motif_umap_wid_celltypes.pdf"))
p
dev.off()

p <- DimPlot(integrated, reduction = "chromumap", group.by = "Patient_tumortype")
dir.create(paste0(savedir, "motif/UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/motif_umap_wid_Patient_tumortype.pdf"))
p
dev.off()

p <- DimPlot(integrated, reduction = "chromumap", group.by = "tumor_type")
dir.create(paste0(savedir, "motif/UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/motif_umap_wid_tumortype.pdf"))
p
dev.off()

p <- DimPlot(integrated, reduction = "chromumap", group.by = "celltypes", split.by = "celltypes", ncol = 3)

dir.create(paste0(savedir, "motif/UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/motif_umap_wid_celltypes_splitted.pdf"))
p
dev.off()

p <- DimPlot(integrated, reduction = "chromumap", group.by = "Patient_tumortype", split.by = "Patient_tumortype", ncol = 3)
dir.create(paste0(savedir, "motif/UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/motif_umap_wid_Patient_tumortype_splitted.pdf"))
p
dev.off()

p <- DimPlot(integrated, reduction = "chromumap", group.by = "tumor_type", split.by = "tumortype", ncol = 2)
dir.create(paste0(savedir, "motif/UMAP/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/motif_umap_wid_tumortype_splitted.pdf"))
p
dev.off()

## Adding the batch information
sam_batch <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/sample_batches.txt", header = TRUE, )
sam_batch$samplename <- gsub("-128.*.", "", sam_batch$Sample) %>% gsub("-", "_", .)

patterns <- sam_batch$samplename
replacements <- sam_batch$Batch

integrated@meta.data$Batch <- integrated@meta.data$Sample

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    integrated@meta.data$Batch <- str_replace_all(integrated@meta.data$Batch, pattern, replacements[i])
}

### Trying to perform integration on chromvar matrix
# Install Harmony if not already installed
library(harmony)
scATAC <- ScaleData(scATAC, features = rownames(scATAC))

scATAC <- RunPCA(
    scATAC,
    assay = "chromvar",
    features = rownames(scATAC),
    reduction.name = "chrompca",
    reduction.key = "chromPC_"
)

scATAC <- RunHarmony(
    scATAC,
    group.by.vars = "Batch",
    reduction.use = "chrompca",
    assay.use = "chromvar",
    reduction.save = "chromharmony"
)

scATAC <- RunUMAP(
    scATAC,
    assay = "chromvar",
    reduction.key = "chromUMAP_",
    reduction = "chromharmony",
    reduction.name = "chromharmonyumap",
    dims = 2:30
)

stopifnot(all(rownames(integrated@meta.data) == row.names(scATAC@meta.data)))

p1 <- DimPlot(scATAC, reduction = "chromharmonyumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(scATAC, reduction = "chromharmonyumap", label = TRUE, group.by = "celltypes")
p3 <- DimPlot(scATAC, reduction = "chromharmonyumap", label = TRUE, group.by = "Batch")

pdf(paste0(savedir, "motif/chromvar/UMAP/batches_chromharmony_sample_celltypes_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

p1 <- DimPlot(scATAC, reduction = "chromumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(scATAC, reduction = "chromumap", label = TRUE, group.by = "celltypes")
p3 <- DimPlot(scATAC, reduction = "chromumap", label = TRUE, group.by = "Batch")

pdf(paste0(savedir, "motif/chromvar/UMAP/batches_chrom_sample_celltypes_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

### Integrating the scATAC
### Harmony Integration
library(Signac)
library(Seurat)
library(harmony)

combined <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/combined.RDS")
combined@meta.data$Sample <- gsub("_128.*.", "", rownames(combined@meta.data))

DefaultAssay(combined) <- "ATAC"

combined_harmony <- RunHarmony(
    combined,
    group.by.vars = "Sample",
    reduction.use = "lsi",
    assay.use = "ATAC",
    reduction.save = "ATACharmony",
    project.dim = FALSE
)

combined_harmony <- RunUMAP(
    combined_harmony,
    assay = "ATAC",
    reduction.key = "ATACharmonyUMAP_",
    reduction = "ATACharmony",
    reduction.name = "ATACharmonyumap",
    dims = 2:30
)

stopifnot(all(rownames(combined_harmony@meta.data) == rownames(scATAC@meta.data)))
combined_harmony@meta.data$celltypes <- scATAC@meta.data$celltypes
combined_harmony@meta.data$Batch <- scATAC@meta.data$Batch

p1 <- DimPlot(combined_harmony, reduction = "ATACharmonyumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(combined_harmony, reduction = "ATACharmonyumap", label = TRUE, group.by = "celltypes")
p3 <- DimPlot(combined_harmony, reduction = "ATACharmonyumap", label = TRUE, group.by = "Batch")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste0(savedir, "UMAP/Sample_ATACharmony_sample_celltypes_batch.pdf"), width = 19, height = 5.5)
p1 + p2 + p3
dev.off()

p4 <- DimPlot(combined_harmony, reduction = "ATACharmonyumap", label = FALSE, group.by = "Sample", split.by = "Sample", ncol = 4)
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste0(savedir, "UMAP/Sample_ATACharmony_sample_splitted.pdf"), width = 15, height = 15)
p4
dev.off()

# saveRDS(combined_harmony, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/ATAC_harmony.RDS")
saveRDS(scATAC, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_chromvar_3.RDS")

### Performing based on the batches
stopifnot(all(rownames(combined@meta.data) == rownames(scATAC@meta.data)))
combined@meta.data$Batch <- scATAC@meta.data$Batch

harmony_batch <- RunHarmony(
    combined,
    group.by.vars = "Batch",
    reduction.use = "lsi",
    assay.use = "ATAC",
    reduction.save = "ATACBatchharmony",
    project.dim = FALSE
)

harmony_batch <- RunUMAP(
    harmony_batch,
    assay = "ATAC",
    reduction.key = "ATACBatchharmonyUMAP_",
    reduction = "ATACBatchharmony",
    reduction.name = "ATACBatchharmonyumap",
    dims = 2:30
)

stopifnot(all(rownames(harmony_batch@meta.data) == rownames(scATAC@meta.data)))
harmony_batch@meta.data$celltypes <- scATAC@meta.data$celltypes
harmony_batch@meta.data$Sample <- scATAC@meta.data$Sample

p1 <- DimPlot(harmony_batch, reduction = "ATACBatchharmonyumap", label = TRUE, group.by = "Sample")
p2 <- DimPlot(harmony_batch, reduction = "ATACBatchharmonyumap", label = TRUE, group.by = "celltypes")
p3 <- DimPlot(harmony_batch, reduction = "ATACBatchharmonyumap", label = TRUE, group.by = "Batch")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste0(savedir, "UMAP/ATAC_batches_harmony_sample_celltypes_batch.pdf"), width = 15, height = 5.5)
p1 + p2 + p3
dev.off()

p4 <- DimPlot(harmony_batch, reduction = "ATACBatchharmonyumap", label = FALSE, group.by = "Sample", split.by = "Sample", ncol = 3)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste0(savedir, "UMAP/ATAC_batches_hamony_batch_splitted.pdf"), width = 15, height = 15)
p4
dev.off()

### Harmony Batch integration show high batch effect between samples

### Footpringting
# gather the footprinting information for sets of motifs
integrated <- Footprint(
    object = integrated,
    assay = "ATAC",
    motif.name = c("SP1", "SP3", "KLF2", "HIF1A", "KLF2", "EGR1", "HSF1"),
    genome = BSgenome.Hsapiens.UCSC.hg38,
    # regions = StringToGRanges(peaks[1]),
    # key = "VEGFA_Footprint_"
)

# plot the footprint data for each group of cells
p <- PlotFootprint(integrated, features = c("SP1", "SP3", "KLF2", "HIF1A", "KLF2", "EGR1", "HSF1"), assay = "ATAC", label.top = 5, show.expected = FALSE, group.by = "tumor_type")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
dir.create(paste0(savedir, "motif/footprint/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/footprint/DNA_footprint_interested_TFs_.pdf"), width = 20, height = 20)
p
dev.off()

integrated <- Footprint(
    object = integrated,
    assay = "ATAC",
    motif.name = c("NFYA", "NFYB", "NFYC", "KLF14", "YY2", "NFIA", "NFIC", "HIC2", "SREBF1"),
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p <- PlotFootprint(integrated,
    features = c("NFYA", "NFYB", "NFYC", "KLF14", "YY2", "NFIA", "NFIC", "HIC2", "SREBF1"),
    assay = "ATAC", label.top = 5, show.expected = TRUE, group.by = "celltypes"
)

dir.create(paste0(savedir, "motif/footprint/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/footprint/primary_DNA_footprint_interested_NFY.pdf"), width = 20, height = 20)
p
dev.off()

integreated_paired_mal <- Footprint(
    object = integreated_paired_mal,
    assay = "ATAC",
    motif.name = c("NFYA", "NFYB", "NFYC", "KLF14", "YY2", "NFIA", "NFIC", "HIC2", "SREBF1"),
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p <- PlotFootprint(integreated_paired_mal, features = c("NFYA", "NFYB", "NFYC", "KLF14", "YY2", "NFIA", "NFIC", "HIC2", "SREBF1"), assay = "ATAC", label.top = 5, show.expected = FALSE, group.by = "celltypes")

dir.create(paste0(savedir, "motif/footprint/"), showWarnings = FALSE)
pdf(paste0(savedir, "motif/footprint/primary_paired_malignant_DNA_footprint_interested_TFs.pdf"), width = 20, height = 20)
p
dev.off()

# We can also directly test for differential activity scores between cell types. This tends to give similar results as performing an enrichment test on differentially accessible peaks between the cell types (shown above).
# When performing differential testing on the chromVAR z-score, we can set mean.fxn=rowMeans and fc.name="avg_diff" in the FindMarkers() function so that the fold-change calculation computes the average difference in z-score between the groups.

DefaultAssay(integrated) <- "chromvar"
integrated@meta.data$Patient_tumortype <- paste(integrated@meta.data$Patient.ID, integrated@meta.data$tumor_type, sep = "_")

paired_samples <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"

dir.create(paste0(savedir, "motif/scdifferential/"), showWarnings = FALSE)

for (i in 1:length(paired_samples)) {
    differential.activity <- FindMarkers(
        object = integrated,
        ident.1 = paste0(paired_samples[i], "_Recurrent"),
        ident.2 = paste0(paired_samples[i], "_Primary"),
        only.pos = FALSE,
        mean.fxn = rowMeans,
        fc.name = "avg_diff",
        group.by = "Patient_tumortype"
    )
    assign(paste0(paired_samples[i], "_motif_rec_vs_pri"), differential.activity)
    write.table(
        differential.activity,
        paste0(
            savedir, "motif/scdifferential/",
            paired_samples[i], "_rec_vs_pri_motif_both.txt"
        ),
        quote = F, row.names = T,
        col.names = T, sep = "\t"
    )
}

## We have identify the chromvar motif id and TF name
jaspar_id_name <- read.table("/diazlab/data3/.abhinav/resources/pfms/JASPAR2020_CORE_non-redundant_pfms_jaspar/final_JASPAR_MAid_TFnames.txt")
colnames(jaspar_id_name) <- c("motif_id", "motif_name")

diff_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar",
    pattern = "rec_vs_pri_motif_both.txt",
    full.names = TRUE,
    recursive = TRUE
)

samplename <- gsub("_rec_vs_pri_motif_both.txt", "", basename(diff_files))
for (i in 1:length(diff_files)) {
    motif_diff <- read.table(diff_files[i])
    motif_diff_sig <- motif_diff[motif_diff$p_val_adj < 0.05 & (motif_diff$avg_diff > 0.5 | motif_diff$avg_diff < -0.5), ]
    motif_diff_sig$TF_name <- jaspar_id_name[match(rownames(motif_diff_sig), jaspar_id_name$motif_id), "motif_name"]
    assign(paste0(samplename[i], "_motif_sig"), motif_diff_sig)

    motif_diff_sig_rec <- motif_diff_sig[motif_diff_sig$avg_diff > 0.5, ]
    assign(paste0(samplename[i], "_motif_sig_rec"), motif_diff_sig_rec)

    motif_diff_sig_pri <- motif_diff_sig[motif_diff_sig$avg_diff < -0.5, ]
    assign(paste0(samplename[i], "_motif_sig_pri"), motif_diff_sig_pri)

    write.table(
        motif_diff_sig,
        paste0(savedir, "motif/scdifferential/", samplename[i], "_rec_vs_pri_motif_sig_avg_diff_0.5.txt"),
        quote = F, row.names = T, col.names = T, sep = "\t"
    )
}

### Performing to make the Heatmap
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/"
jaspar_id_name <- read.table("/diazlab/data3/.abhinav/resources/pfms/JASPAR2020_CORE_non-redundant_pfms_jaspar/final_JASPAR_MAid_TFnames.txt")
colnames(jaspar_id_name) <- c("motif_id", "motif_name")

diff_files <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar",
    pattern = "rec_vs_pri_motif_both.txt",
    full.names = TRUE,
    recursive = TRUE
)

samplename <- gsub("_rec_vs_pri_motif_both.txt", "", basename(diff_files))
for (i in 1:length(diff_files)) {
    motif_diff <- read.table(diff_files[i])
    motif_diff_sig <- motif_diff[motif_diff$p_val_adj < 0.05 & (motif_diff$avg_diff > 0.2 | motif_diff$avg_diff < -0.2), ] ## Keeping it very lenient for heatmap
    motif_diff_sig$TF_name <- jaspar_id_name[match(rownames(motif_diff_sig), jaspar_id_name$motif_id), "motif_name"]
    assign(paste0(samplename[i], "_motif_sig"), motif_diff_sig)

    motif_diff_sig_rec <- motif_diff_sig[motif_diff_sig$avg_diff > 0.2, ]
    assign(paste0(samplename[i], "_motif_sig_rec"), motif_diff_sig_rec)

    motif_diff_sig_pri <- motif_diff_sig[motif_diff_sig$avg_diff < -0.2, ]
    assign(paste0(samplename[i], "_motif_sig_pri"), motif_diff_sig_pri)

    write.table(
        motif_diff_sig,
        paste0(savedir, samplename[i], "_rec_vs_pri_motif_sig_avg_diff_0.2_padj_0.5.txt"),
        quote = F, row.names = T, col.names = T, sep = "\t"
    )
}

motif_sig <- ls(pattern = "_motif_sig_rec")
library(data.table)
PT_7WYPEC3Q <- get(motif_sig[1])[, "TF_name"]
PT_9S6WMQ92 <- get(motif_sig[2])[, "TF_name"]
PT_H3WWDMW9 <- get(motif_sig[3])[, "TF_name"]
PT_XA98HG1C <- get(motif_sig[4])[, "TF_name"]

list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q, "PT_9S6WMQ92" = PT_9S6WMQ92,
    "PT_H3WWDMW9" = PT_H3WWDMW9, "PT_XA98HG1C" = PT_XA98HG1C
)

n <- max(
    length(PT_7WYPEC3Q), length(PT_9S6WMQ92), length(PT_H3WWDMW9), length(PT_XA98HG1C)
)

length(PT_7WYPEC3Q) <- n
length(PT_9S6WMQ92) <- n
length(PT_H3WWDMW9) <- n
length(PT_XA98HG1C) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q), (PT_9S6WMQ92), (PT_H3WWDMW9), (PT_XA98HG1C)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_motif_sig", "", motif_sig)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/"
pdf(paste(savedir, "Recurrent_motif_enriched_significant_avg_diff_0.2_padj_ls_0.05.pdf", sep = ""), width = 12, height = 7)
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
motif_filter_common <- df_int[(df_int$sample_num > 2), 1]
motif_filter_all <- df_int[(df_int$sample_num > 3), 1]

write.table(df_int, paste(savedir, "rec_motif_intersect_name_padj_0.05_lfc_0.2.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "rec_motif_intersect_gene_number_padj_0.05_lfc_0.2.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(motif_filter_common, paste(savedir, "rec_motif_intersect_markers_morethan2_padj_0.05_lfc_0.2.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
write.table(motif_filter_all, paste(savedir, "rec_motif_intersect_markers_all_padj_0.05_lfc_0.2.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

motif_high_recurrent_all <- jaspar_id_name[match(motif_filter_all,jaspar_id_name$motif_name),"motif_id"]
p <- DoHeatmap(integreated_paired_mal, features = motif_high_recurrent_all, group.by = "tumor_type") 
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_Rec_high_in_all_samples_motif_tumor_type_padj_0.05.pdf")
p
dev.off()

motif_high_recurrent_all <- jaspar_id_name[match(motif_filter_all,jaspar_id_name$motif_name),"motif_id"]
p <- DoHeatmap(integreated_paired_mal, features = motif_high_recurrent_all, group.by = "Patient_tumortype")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_Rec_all_high_motif_Patient_tumor_type_padj_0.05.pdf")
p
dev.off()

motif_high_recurrent <- jaspar_id_name[match(motif_filter_common,jaspar_id_name$motif_name),"motif_id"]
p <- DoHeatmap(integreated_paired_mal, features = motif_high_recurrent, group.by = "tumor_type")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_Rec_high_motif_tumor_type_padj_0.05.pdf")
p
dev.off()

motif_high_recurrent <- jaspar_id_name[match(motif_filter_common,jaspar_id_name$motif_name),"motif_id"]
p <- DoHeatmap(integreated_paired_mal, features = motif_high_recurrent, group.by = "Patient_tumortype")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_Rec_high_motif_Patient_tumor_type_padj_0.05.pdf")
p
dev.off()

motif_high_recurrent <- jaspar_id_name[match(motif_filter_common,jaspar_id_name$motif_name),"motif_id"]

p <- DotPlot(integreated_paired_mal, 
            motif_high_recurrent,
            assay = "chromvar",
            group.by = "Patient_tumortype") +
            scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/dotplot_common_Rec_high_motif_Patient_tumor_type_padj_0.05.pdf", height = 15, width = 6)
p
dev.off()

#### SP and KLF
# limits = c(-1.5,1.5)
motif_filter_common_SP_KLF <- grep("SP|KLF|EGR",motif_filter_common,value=TRUE)
motif_filter_common_SP_KLF_order <- motif_filter_common_SP_KLF[order(motif_filter_common_SP_KLF)]
motif_recurrent_SP_KLF <- jaspar_id_name[match(motif_filter_common_SP_KLF_order, jaspar_id_name$motif_name),]

p <- DotPlot(integreated_paired_mal, 
            c(motif_recurrent_SP_KLF$motif_id,"MA1106.1"), 
            assay = "chromvar",
            # layer = "data",
            group.by = "Patient_tumortype") +
            scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue2","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/dotplot/common_Rec_high_SP_KLF_EGR_Patient_tumor_type_padj_0.05.pdf", height = 8)
p
dev.off()

motif_filter_common_FOX <- motif_filter_common[grep("FOX",motif_filter_common)]
motif_FOX_id <- jaspar_id_name[match(motif_filter_common_FOX, jaspar_id_name$motif_name),]

p <- DotPlot(integreated_paired_mal, 
            motif_FOX_id$motif_id, 
            assay = "chromvar",
            # layer = "data",
            group.by = "Patient_tumortype") +
            scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
            # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/dotplot/common_Rec_high_FOX_Patient_tumor_type_padj_0.05.pdf")
p
dev.off()

p <- DotPlot(integreated_paired_mal, 
            c(motif_recurrent_SP_KLF$motif_id,"MA1106.1"), 
            assay = "chromvar",
            # layer = "data",
            group.by = "celltypes") +
            scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) + 
            # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/dotplot/celltypes_KLF_SP_common_Rec_high.pdf")
p
dev.off()

p <- DotPlot(integreated_paired_mal, 
            motif_FOX_id$motif_id, 
            assay = "chromvar",
            # layer = "data",
            group.by = "predicted.id") +
            scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
            # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/dotplot/predicted_common_Rec_high_FOX.pdf")
p
dev.off()

genes <- c("VEGFA","RICTOR","AKT3","PIK3CA","FOXO1","CCND2","FOXO3",
           "EIF4EBP1","EIF2AK3","EIF2A","EIF3E","EIF2AK4","XBP1","ATF6","HSPB1","CDKN1B")
p <- DotPlot(integreated_paired_mal,
            genes, 
            assay = "RNA",
            # layer = "data",
            group.by = "celltypes") +
            scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
            # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/dotplot/Gene_RNA_VEGFA_ER.pdf")
p
dev.off()

#### Making a Motifplot
DefaultAssay(integreated_paired_mal) <- "ATAC"
p <- MotifPlot(
  object = integreated_paired_mal,
  motifs = c(motif_recurrent_SP_KLF$motif_id,"MA1106.1")
)

dir.create("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/Motif_plot/",showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/Motif_plot/SP_KLF_EGR_common_Rec_high.pdf", width = 9)
p
dev.off()


motif_filter_common_FOX <- motif_filter_common[grep("FOX",motif_filter_common)]
motif_all_pos <- jaspar_id_name[match(motif_filter_all, jaspar_id_name$motif_name),]

p <- DotPlot(integreated_paired_mal, 
            motif_all_pos$motif_id, 
            assay = "chromvar",
            # layer = "data",
            group.by = "Patient_tumortype") +
            scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
            # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            coord_flip() +
            scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/dotplot/common_Rec_high_all_high_Patient_tumor_type_padj_0.05_all.pdf")
p
dev.off()


rm(plot_list)
plot_list <- list()
KLF_SP_motifs <- motif_recurrent_SP_KLF$motif_id
for(i in 1:length(KLF_SP_motifs)){
    plot_list[[i]] <- VlnPlot(integreated_paired_mal, 
            KLF_SP_motifs[i], 
            assay = "chromvar",
            pt.size = 0,
            slot = "scale.data",
            group.by = "Patient_tumortype") + 
            geom_boxplot()
            # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
            # coord_flip() +
            # scale_size(range = c(1, 10)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

}

dir.create("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/vlnplot",showWarnings = FALSE)
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/vlnplot/common_Rec_high_SP_KLF_Patient_tumor_type_padj_0.05_scaledata.pdf")
plot_list
dev.off()



### for Primary TF
motif_sig <- ls(pattern = "_motif_sig_pri")
library(data.table)
PT_7WYPEC3Q <- get(motif_sig[1])[, "TF_name"]
PT_9S6WMQ92 <- get(motif_sig[2])[, "TF_name"]
PT_H3WWDMW9 <- get(motif_sig[3])[, "TF_name"]
PT_XA98HG1C <- get(motif_sig[4])[, "TF_name"]

list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q, "PT_9S6WMQ92" = PT_9S6WMQ92,
    "PT_H3WWDMW9" = PT_H3WWDMW9, "PT_XA98HG1C" = PT_XA98HG1C
)

n <- max(
    length(PT_7WYPEC3Q), length(PT_9S6WMQ92), length(PT_H3WWDMW9), length(PT_XA98HG1C)
)

length(PT_7WYPEC3Q) <- n
length(PT_9S6WMQ92) <- n
length(PT_H3WWDMW9) <- n
length(PT_XA98HG1C) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q), (PT_9S6WMQ92), (PT_H3WWDMW9), (PT_XA98HG1C)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_motif_sig", "", motif_sig)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/"
pdf(paste(savedir, "Primary_motif_enriched_pvalue_0.5_lfc_0.2.pdf", sep = ""), width = 12, height = 7)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
library(dplyr)
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

final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

library(stringr)
df_int$int2 <- gsub("\\|", "__", df_int$int)
# GCP_common <- df_int[grep("mar_7316_2118__mar_7316_278__mar_7316_2978__mar_7316_3023__mar_7316_311__mar_7316_333__mar_7316_4529__mar_7316_5881__mar_7316_737__mar_7316_931__mar_DOD4182__mar_SF10961__mar_SF12930__mar_SF7994__mar_SF8368__mar_SF8539",df_int$int2),1]

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
motif_filter_common <- df_int[(df_int$sample_num > 2), 1]
motif_filter_all <- df_int[(df_int$sample_num > 3), 1]

write.table(df_int, paste(savedir, "pri_motif_intersect_name_lfc_0.2.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "pri_motif_intersect_gene_number_lfc_0.2.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(motif_filter_common, paste(savedir, "pri_motif_intersect_markers_morethan2_lfc_0.2.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
write.table(motif_filter_all, paste(savedir, "pri_motif_intersect_markers_all_lfc_0.2.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

motif_high_pri <- jaspar_id_name[match(motif_filter_common,jaspar_id_name$motif_name),"motif_id"]
p <- DoHeatmap(integreated_paired_mal, features = motif_high_pri, group.by = "tumor_type") + scale_fill_gradientn(colors = c("blue", "white", "red"))
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_prim_high_motif_tumor_type.pdf")
p
dev.off()

motif_high_pri_all <- jaspar_id_name[match(motif_filter_all,jaspar_id_name$motif_name),"motif_id"]
p <- DoHeatmap(integreated_paired_mal, features = motif_high_pri_all, group.by = "tumor_type") + scale_fill_gradientn(colors = c("blue", "white", "red"))
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_prim_in_all_samples_motif_tumor_type.pdf")
p
dev.off()

### combining all the primary and recurrent motifs
recurrent_all_motif <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/rec_motif_intersect_markers_all_padj_0.5_lfc_0.txt")[,1]
recurrent_all_motif_id <- jaspar_id_name[match(recurrent_all_motif,jaspar_id_name$motif_name),"motif_id"]

recurrent_3_motif <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/rec_motif_intersect_markers_morethan2_padj_0.5_lfc_0.txt")[,1]
recurrent_3_motif_id <- jaspar_id_name[match(recurrent_3_motif,jaspar_id_name$motif_name),"motif_id"]

recurrent_3_motif_id_spec <- recurrent_3_motif_id[grep(paste(recurrent_all_motif_id, sep ="",collapse = "|"), recurrent_3_motif_id,invert=TRUE)]

primary_all_motif <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/pri_motif_intersect_markers_all_lfc_0.txt")[,1]
primary_all_motif_id <- jaspar_id_name[match(primary_all_motif,jaspar_id_name$motif_name),"motif_id"]

primary_3_motif <- read.table("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/pri_motif_intersect_markers_morethan2_lfc_0.txt")[,1]
primary_3_motif_id <- jaspar_id_name[match(primary_3_motif,jaspar_id_name$motif_name),"motif_id"]

primary_3_motif_id_spec <- primary_3_motif_id[grep(paste(primary_all_motif_id, sep ="",collapse = "|"), primary_3_motif_id,invert=TRUE)]

combined_motif <- c(recurrent_all_motif_id, recurrent_3_motif_id_spec,primary_all_motif_id, primary_3_motif_id_spec)

DefaultAssay(integreated_paired_mal) <- "chromvar"
p <- DoHeatmap(integreated_paired_mal, features = combined_motif, group.by = "tumor_type") + scale_fill_gradientn(colors = c("blue", "white", "red"))
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_prim_and_rec_in_all_samples_motif_tumor_type.pdf")
p
dev.off()

p <- DoHeatmap(integreated_paired_mal, features = combined_motif, group.by = "tumor_type")
pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/chromvar/common_prim_and_rec_in_all_samples_motif_tumor_type_nocolor.pdf")
p
dev.off()


### Finding the common motif in recurrent vs primary
### Common Motifs
motif_sig <- ls(pattern = "_motif_sig")

library(data.table)
PT_7WYPEC3Q <- rownames(get(motif_sig[1]))
PT_9S6WMQ92 <- rownames(get(motif_sig[2]))
PT_H3WWDMW9 <- rownames(get(motif_sig[3]))
PT_XA98HG1C <- rownames(get(motif_sig[4]))

list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q, "PT_9S6WMQ92" = PT_9S6WMQ92,
    "PT_H3WWDMW9" = PT_H3WWDMW9, "PT_XA98HG1C" = PT_XA98HG1C
)

n <- max(
    length(PT_7WYPEC3Q), length(PT_9S6WMQ92), length(PT_H3WWDMW9), length(PT_XA98HG1C)
)

length(PT_7WYPEC3Q) <- n
length(PT_9S6WMQ92) <- n
length(PT_H3WWDMW9) <- n
length(PT_XA98HG1C) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q), (PT_9S6WMQ92), (PT_H3WWDMW9), (PT_XA98HG1C)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function

colnames(fld) <- gsub("_motif_sig", "", motif_sig)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "motif/Recurrent_motif_enriched_significant.pdf", sep = ""), width = 12, height = 7)
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
motif_filter_common <- df_int[(df_int$sample_num > 2), 1]
motif_filter_all <- df_int[(df_int$sample_num > 3), 1]

write.table(df_int, paste(savedir, "motif/motif_intersect_name.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "motif/motif_intersect_gene_number.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(motif_filter_common, paste(savedir, "motif/motif_intersect_markers_morethan2.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)
write.table(motif_filter_all, paste(savedir, "motif/motif_intersect_markers_all.txt", sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)


### Quick Checking of VEGFA and HIF1A
# Priamry and Recurrent
genes <- c(
    "EIF3E", "EIF4EBP1", "EIF2A", "EIF2AK3", "ATF4", "ATF6", "NFE2L2", "EIF2S1", "AKT1",
    "FOXO1", "RICTOR", "CDKN1B", "AKT3", "VEGFA", "HSPB1", "CCND2", "TP53", "RPTOR", "PTEN",
    "FLT1", "FLT4", "KDR", "FOXO4", "FOXO6", "SP3", "SP1", "KLF2", "EGR1", "EGR3", "MMP9"
)

integrated@meta.data$patient_tumortype <- paste(integrated@meta.data$Patient.ID, integrated@meta.data$tumor_type, sep = "_")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(integrated, genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/ER_stress_VEGF_pathway_recurrent_genes_malignant.pdf", sep = ""))
plot_list
dev.off()

paired_cells <- row.names(integrated@meta.data[grep(c("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931"), integrated@meta.data$Sample), ])
integrated_paired <- subset(integrated, cells = paired_cells)
DefaultAssay(integrated_paired) <- "RNA"
integrated_paired <- NormalizeData(integrated_paired)

genes <- c(
    "EIF3E", "EIF4EBP1", "EIF2A", "EIF2AK3", "ATF4", "ATF6", "NFE2L2", "EIF2S1", "AKT1",
    "FOXO1", "RICTOR", "CDKN1B", "AKT3", "VEGFA", "HIF1A", "HSPB1", "CCND2", "TP53", "RPTOR", "PTEN",
    "FLT1", "FLT4", "KDR"
)

# integrated@meta.data$patient_tumortype <- paste(integrated@meta.data$Patient.ID, integrated@meta.data$tumor_type, sep = "_")
DefaultAssay(integrated_paired) <- "MAGIC_RNA"
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
    plot_list[[i]] <- VlnPlot(integrated_paired, genes[i], pt.size = 0, group.by = "patient_tumortype", assay = "MAGIC_RNA") + geom_boxplot()
}

dir.create(paste(savedir, "vlnplot", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "vlnplot/ER_stress_VEGF_pathway_recurrent_genes_paired.pdf", sep = ""))
plot_list
dev.off()

### Performing Differential in scATAC for RNA
library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header = TRUE, sep = "\t")

### Patient ID
patterns <- long_sample$Sample.ID
patterns <- gsub("SF7994R", "SF7994", patterns)
patterns <- gsub("-", "_", patterns)
replacements <- long_sample$Patient.ID
# names(replacement) <- patterns

integrated@meta.data$Sample <- gsub("sample_", "", integrated@meta.data$dataset) %>% gsub("_128.*.", "", .)
integrated@meta.data$Patient.ID <- integrated@meta.data$Sample

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    integrated@meta.data$Patient.ID <- str_replace_all(integrated@meta.data$Patient.ID, pattern, replacements[i])
}

replacements <- long_sample$Tumor.type
integrated@meta.data$tumor_type <- integrated@meta.data$Sample

for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    integrated@meta.data$tumor_type <- str_replace_all(integrated@meta.data$tumor_type, pattern, replacements[i])
}

R_vs_P_all <- FindMarkers(integrated_paired, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type")
dir.create(paste(savedir, "scdifferential", sep = ""), showWarnings = FALSE)
write.table(R_vs_P_all, paste(savedir, "scdifferential/R_vs_P_all_lognorm.txt", sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)

#### USING MAST
R_vs_P_MAST <- FindMarkers(integrated_paired, ident.1 = "Recurrent", ident.2 = "Primary", group.by = "tumor_type", test.use = "MAST")
write.table(R_vs_P_all, paste(savedir, "scdifferential/R_vs_P_all_MAST.txt", sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)

### Single cell Peak Calling
library(Signac)
library(Seurat)

pbmc <- readRDS("../vignette_data/pbmc.rds")
DimPlot(pbmc)

peaks <- CallPeaks(
    object = pbmc,
    group.by = "predicted.id"
)


### Linking Peaks to Genes
# For each gene, we can find the set of peaks that may regulate the gene by by computing the correlation between gene expression and accessibility at nearby
# peaks, and correcting for bias due to GC content, overall accessibility, and peak size. See the Signac paper for a full description of the method we use to
# link peaks to genes.

# Running this step on the whole genome can be time consuming, so here we demonstrate peak-gene links for a subset of genes as an example. The same function can
#  be used to find links for all genes by omitting the genes.use parameter:

DefaultAssay(integrated) <- "ATAC"

# first compute the GC content for each peak
integrated <- RegionStats(integrated, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
integrated <- LinkPeaks(
    object = integrated,
    peak.assay = "ATAC",
    expression.assay = "SCT",
)

p1 <- CoveragePlot(
    object = integrated,
    region = "VEGFA",
    features = "VEGFA",
    assay = "ATAC",
    expression.assay = "MAGIC_RNA",
    group.by = "Patient_tumortype",
    peaks = TRUE,
    links = TRUE,
    extend.upstream = 50000,
    extend.downstream = 50000
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_MAGIC_RNA_specific_50K.pdf", sep = ""))
p1
dev.off()


VEGFA_integrated <- LinkPeaks(
    object = integrated,
    peak.assay = "ATAC",
    expression.assay = "SCT",
    genes.use = "VEGFA"
)

p1 <- CoveragePlot(
    object = integrated,
    region = "VEGFA",
    features = "VEGFA",
    assay = "ATAC",
    expression.assay = "MAGIC_RNA",
    group.by = "Patient_tumortype",
    peaks = TRUE,
    links = TRUE,
    extend.upstream = 50000,
    extend.downstream = 50000
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_MAGIC_RNA_specific_50K.pdf", sep = ""))
p1
dev.off()


# RNA imputation

# Above we transferred categorical information (the cell labels) and mapped the query data onto an existing reference UMAP. We can also transfer continuous data
# from the reference to the query in the same way. Here we demonstrate transferring the gene expression values from the PBMC multiome dataset (that measured DNA
# accessibility and gene expression in the same cells) to the PBMC scATAC-seq dataset (that measured DNA accessibility only). Note that we could also transfer
# these values using the MapQuery() function call above by setting the refdata parameter to a list of values.
library(Signac)
library(Seurat)
library(reticulate)
library(Rmagic)

ATAC_integrated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_chromvar_4.RDS")
RNA_integraated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed.RDS")
transfer.anchors <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/anchors_snRNA_integrated_cca.RDS")

# map query onto the reference dataset
# compute UMAP and store the UMAP model
RNA_integraated <- RunUMAP(RNA_integraated, reduction = "pca", dims = 1:30, return.model = TRUE)

ATAC_integrated <- MapQuery(
    anchorset = transfer.anchors,
    reference = RNA_integraated,
    query = ATAC_integrated,
    refdata = RNA_integraated$sub_celltypes,
    reference.reduction = "pca",
    new.reduction.name = "ref.pca",
    reduction.model = "umap"
)

# predict gene expression values
imputed_rna <- TransferData(
    anchorset = transfer.anchors,
    refdata = LayerData(RNA_integraated, assay = "RNA", layer = "data"),
    weight.reduction = ATAC_integrated[["integrated_lsi"]],
    dims = 2:30
)

ATAC_integrated[["imputed_rna"]] <- imputed_rna

p1 <- DimPlot(RNA_integraated, reduction = "umap", group.by = "sub_celltypes", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("snRNA")
p2 <- DimPlot(ATAC_integrated, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scATAC")
p3 <- DimPlot(ATAC_integrated, reduction = "ref.umap", group.by = "celltypes", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scATAC Celltypes")

pdf(paste(savedir, "UMAP/scRNA_scATAC_imputed.pdf", sep = ""), width = 16, height = 5.5)
p1 | p2 | p3
dev.off()

DefaultAssay(ATAC_integrated) <- "imputed_rna"
use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

DefaultAssay(ATAC_integrated) <- "RNA"
ATAC_integrated <- magic(ATAC_integrated, npca = 30) ## imputing the RNA data as for RNA PCs are 20

saveRDS(ATAC_integrated, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_chromvar.RDS")

### Performing Coembeddding
ATAC_integrated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_QC_snRNA_ref_mapped_chromvar.RDS")
coembed <- merge(x = RNA_integrated, y = ATAC)

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Coembed/coembedding.pdf", width = 12, height = 5.5)
DimPlot(coembed, group.by = c("orig.ident", "celltypes")) + NoLegend()
dev.off()

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))

integrated <- ATAC_integrated
integrated <- LinkPeaks(
    object = integrated,
    peak.assay = "ATAC",
    expression.assay = "imputed_rna",
)

p1 <- CoveragePlot(
    object = integrated,
    region = "VEGFA",
    features = "VEGFA",
    assay = "ATAC",
    expression.assay = "imputed_rna",
    group.by = "Patient_tumortype",
    extend.upstream = 200000,
    extend.downstream = 200000
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_imputed_pseudo_RNA.pdf", sep = ""))
p1
dev.off()

p1 <- CoveragePlot(
    object = integrated,
    region = "chr6-44072937-44073865",
    # features = "VEGFA",
    assay = "ATAC",
    expression.assay = "MAGIC_imputed_rna",
    group.by = "celltypes",
    extend.upstream = 100,
    extend.downstream = 100
)

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_magic_imputed_pseudo_VEGFA.pdf", sep = ""))
p1
dev.off()

as.data.frame(integrated@assays$ATAC@links[grep("VEGFA", integrated@assays$ATAC@links$gene), ])

### RUnning Cicero
# convert to CellDataSet format and make the cicero object
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)

DefaultAssay(integrated) <- "ATAC"
integrated.cds <- as.cell_data_set(x = integrated)
integrated.cicero <- make_cicero_cds(integrated.cds, reduced_coordinates = reducedDims(integrated.cds)$UMAP)

# Identifying which motif binding to the peaks
library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome)

### The motif is already present in the object
motif_df <- as.data.frame(integrated@assays$ATAC@motifs@data)
motif_names <- as.data.frame(as.matrix(integrated@assays$ATAC@motifs@motif.names))
motif_names

stopifnot(all(rownames(motif_names) == colnames(motif_df)))
colnames(motif_df) <- paste(motif_names$V1, rownames(motif_names), sep = "_")
VEGFA_links <- as.data.frame(integrated@assays$ATAC@links)
VEGFA_links <- as.data.frame(integrated@assays$ATAC@links[grep("^SP1$", integrated@assays$ATAC@links$gene), ])

VEGFA_required_df <- motif_df[match(VEGFA_links$peak, rownames(motif_df)), ]

as.matrix(VEGFA_required_df) -> VEGFA_required_mat

# Convert TRUE to 1 and FALSE to 0
VEGFA_required_df_num <- matrix(as.integer(VEGFA_required_mat), nrow = nrow(VEGFA_required_mat), ncol = ncol(VEGFA_required_mat))
rownames(VEGFA_required_df_num) <- rownames(VEGFA_required_mat)
colnames(VEGFA_required_df_num) <- colnames(VEGFA_required_mat)

VEGFA_columnSum <- as.data.frame(colSums(VEGFA_required_df_num))
VEGFA_columnSum$TF <- gsub("_.*.","",rownames(VEGFA_columnSum))
colnames(VEGFA_columnSum) <- c("peak_number","TF")

# VEGFA_required_df_num_ZF <- VEGFA_required_df_num[,grep("SP[0-9]|KLF|FOX|EGR|HIF1A",colnames(VEGFA_required_df_num))]
VEGFA_required_df_num_ZF <- VEGFA_required_df_num[,match(rownames(VEGFA_columnSum[VEGFA_columnSum$peak_number > 15,]),colnames(VEGFA_required_df_num))]
colnames(VEGFA_required_df_num_ZF) <- gsub("_.*.","",colnames(VEGFA_required_df_num_ZF))
VEGFA_required_df_num_ZF_ordered <- VEGFA_required_df_num_ZF[,order(colnames(VEGFA_required_df_num_ZF))]
# VEGFA_required_df_num_ZF_ordered_2 <- VEGFA_required_df_num_ZF_ordered[,-1] ### Removing Arnh:HIF1A

library(ComplexHeatmap)
p <- Heatmap(t(as.matrix(VEGFA_required_df_num_ZF_ordered)),
        name = "Peak Number",
        col = colorRampPalette(c("white", "mediumturquoise"))(50),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title_position = "topleft", legend_direction = "horizontal"))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/motifmatchrbinding_SP1.pdf", height = 7, width = 7)
draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

rm(plot_list)
plot_list <- list()
peaks <- VEGFA_links$peak
for (i in 1:length(peaks)) {
    plot_list[[i]] <- CoveragePlot(
        object = integrated,
        region = peaks[i],
        # features = "VEGFA",
        assay = "ATAC",
        expression.assay = "imputed_rna",
        group.by = "Patient_tumortype",
        extend.upstream = 1000,
        extend.downstream = 1000
    )
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_imputed_all.pdf", sep = ""))
plot_list
dev.off()

rm(plot_list)
plot_list <- list()
peaks <- VEGFA_links$peak
for (i in 1:length(peaks)) {
    plot_list[[i]] <- CoveragePlot(
        object = integrated,
        region = peaks[i],
        # features = "VEGFA",
        assay = "ATAC",
        expression.assay = "imputed_rna",
        group.by = "celltypes",
        extend.upstream = 1000,
        extend.downstream = 1000
    )
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_celltypes.pdf", sep = ""))
plot_list
dev.off()


rm(plot_list)
plot_list <- list()
peaks <- VEGFA_links$peak
for (i in 1:length(peaks)) {
    plot_list[[i]] <- CoveragePlot(
        object = integrated,
        region = peaks[i],
        # features = "VEGFA",
        assay = "ATAC",
        expression.assay = "imputed_rna",
        group.by = "predicted.id",
        extend.upstream = 1000,
        extend.downstream = 1000
    )
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_imputed_predicted_id.pdf", sep = ""))
plot_list
dev.off()

# ChipSeq Unibind TF DNA interactions
motif_TFs <- read.table("/diazlab/data3/.abhinav/projects/SHH/public_ChipSeq/all_bed/bedtools_intersected/Transcription_factor_peaks.txt", header = FALSE, sep = "\t")
pivot_data <- dcast(motif_TFs, V1 ~ V2, length)
# pivot_data_ordered <- pivot_data[c(1,7,4,2,3,6,5,8),]
rownames(pivot_data) <- pivot_data$V1
pivot_data2 <- pivot_data[,-1]

library(ComplexHeatmap)
p <- Heatmap(as.matrix(pivot_data2), 
        name = "Gene Number", 
        col = colorRampPalette(c("white", "red"))(50),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title_position = "topleft", legend_direction = "horizontal"))

pdf("/diazlab/data3/.abhinav/projects/SHH/public_ChipSeq/all_bed/bedtools_intersected/ChipSeq_TF_motifbinding.pdf", height = 12, width = 7)
draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

## celltype testing
peaks <- c(
    "chr1-31181339-31182772", "chr1-43526867-43528176", "chr1-15153284-15154367",
    "chr1-21638466-21640630", "chr1-110618932-110621379", "chr1-1040372-1041433",
    "chr1-35659906-35661383", "chr1-31413487-31414721", "chr1-226637519-226638884", "chr1-203692700-203693668"
)

rm(plot_list)
plot_list <- list()
# peaks <- VEGFA_links$peak
for (i in 1:length(peaks)) {
    plot_list[[i]] <- CoveragePlot(
        object = integrated,
        region = peaks[i],
        # features = "VEGFA",
        assay = "ATAC",
        expression.assay = "imputed_rna",
        group.by = "predicted.id",
        extend.upstream = 1000,
        extend.downstream = 1000
    )
}

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir, "coverage_plot/VEGFA_link_peaks_imputed_predicted_id_GCP_HSP_test_2.pdf", sep = ""))
plot_list
dev.off()

#### Performing the differential expression from the gene body and promotor peaks
DefaultAssay(integrated) <- "imputed_rna"
patientid <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")

for (i in 1:length(patientid)) {
    group1 <- paste(patientid[i], "_Recurrent", sep = "")
    group2 <- paste(patientid[i], "_Primary", sep = "")

    markers <- FindMarkers(integrated,
        ident.1 = group1,
        ident.2 = group2,
        group.by = "Patient_tumortype",
        test.use = "MAST",
        min.pct = 0,
        logfc.threshold = 0
    )

    assign(paste(patientid[i], "markers_MAST_imputed_rna", sep = "_"), markers)

    write.table(markers,
        paste(savedir, "Table/", patientid[i], "_imputed_rna_recurrent_vs_primary_markers_MAST_nocutoff.txt", sep = ""),
        quote = F,
        row.names = T,
        col.names = T,
        sep = "\t"
    )
}

#### Identifying the common genes in recurrent and Primary
library(UpSetR)
library(dplyr)

filename <- list.files("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/Table",
    pattern = "predicted_rna_recurrent_vs_primary_markers_MAST_nocutoff.txt",
    recursive = TRUE, full.names = TRUE
)

markersname <- gsub("_predicted_rna_recurrent_vs_primary_markers_MAST_nocutoff.txt", "", basename(filename))

### Performing the Peak Annotation for the differential Genes
#### Identifying the common peaks between the samples


# extracting the primary and recurrent markers
for (i in 1:length(markersname)) {
    PTmarkers <- read.table(filename[i], sep = "\t", header = TRUE)
    recurrentmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0.4, ]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
    assign(paste(markersname[i], "_impute_RNA_recurrent", sep = ""), recurrentmarkers)
    primarymarkers <- rownames(PTmarkers[PTmarkers$p_val < 0.05 & PTmarkers$avg_log2FC < -0.4, ])
    assign(paste(markersname[i], "_impute_RNA_primary", sep = ""), primarymarkers)
}

recurrentname <- ls(pattern = "_impute_RNA_recurrent")

H <- list(
    "PT_7WYPEC3Q_impute_RNA_recurrent" = PT_7WYPEC3Q_impute_RNA_recurrent,
    "PT_9S6WMQ92_impute_RNA_recurrent" = PT_9S6WMQ92_impute_RNA_recurrent,
    "PT_H3WWDMW9_impute_RNA_recurrent" = PT_H3WWDMW9_impute_RNA_recurrent,
    "PT_XA98HG1C_impute_RNA_recurrent" = PT_XA98HG1C_impute_RNA_recurrent
)

# create customised venn diagram
library(ggvenn)
a <- ggvenn(H, show_elements = FALSE, stroke_color = "Black", stroke_linetype = "solid")
savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/"
dir.create(paste(savedir, "Table/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "Table/Predicted_RNA_recurrent_high_venny_lfc_0.4_p_0.05.pdf", sep = ""))
a
dev.off()

### Making an upset plot
n <- max(
    length(PT_7WYPEC3Q_impute_RNA_recurrent), length(PT_9S6WMQ92_impute_RNA_recurrent),
    length(PT_H3WWDMW9_impute_RNA_recurrent), length(PT_XA98HG1C_impute_RNA_recurrent)
)

length(PT_7WYPEC3Q_impute_RNA_recurrent) <- n
length(PT_9S6WMQ92_impute_RNA_recurrent) <- n
length(PT_H3WWDMW9_impute_RNA_recurrent) <- n
length(PT_XA98HG1C_impute_RNA_recurrent) <- n

df <- as.data.frame(cbind(
    (PT_7WYPEC3Q_impute_RNA_recurrent), (PT_9S6WMQ92_impute_RNA_recurrent),
    (PT_H3WWDMW9_impute_RNA_recurrent), (PT_XA98HG1C_impute_RNA_recurrent)
))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA
# df$V1 <- NULL
library(UpSetR)
fld <- fromList(as.list(df)) ## this is an intersect function
colnames(fld) <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/scdifferential/"
pdf(paste(savedir, "Table/Upset_Predicted_RNA_recurrent_high_venny_lfc_0.5_p_0.05.pdf", sep = ""), width = 12, height = 7)
upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
dev.off()

### To find the interseected genes
list_filter <- list(
    "PT_7WYPEC3Q" = PT_7WYPEC3Q_impute_RNA_recurrent,
    "PT_9S6WMQ92" = PT_9S6WMQ92_impute_RNA_recurrent,
    "PT_H3WWDMW9" = PT_H3WWDMW9_impute_RNA_recurrent,
    "PT_XA98HG1C" = PT_XA98HG1C_impute_RNA_recurrent
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
library(stringr)
final <- df_int %>%
    group_by(int) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n))

df_int$int2 <- gsub("\\|", "__", df_int$int)

### We will be considering all the genes that are coming in 14,15 and 16 samples
df_int$sample_num <- str_count(df_int$int, "PT_")
all_sample <- df_int[(df_int$sample_num > 3), 1]
morethan_3 <- df_int[(df_int$sample_num >= 3), 1]

write.table(df_int, paste(savedir, "Table/recurrent_predicted_rna_intersect_lfc_0.5_pval_adj_0.05.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir, "Table/recurrent_intersected_predicted_rna_number_lfc_0.5_pval_adj_0.05.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

write.table(all_sample, paste(savedir, "Table/predicted_rna_recurrent_all_samples_lfc_0.5_pval_adj_0.05.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(morethan_3, paste(savedir, "Table/predicted_rna_recurrent_morethan_3_lfc_0.5_pval_adj_0.05.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")


#### Performing the differential expression from the gene body and promotor peaks
DefaultAssay(integrated) <- "RNA"
patientid <- c("PT_7WYPEC3Q", "PT_9S6WMQ92", "PT_H3WWDMW9", "PT_XA98HG1C")

for (i in 1:length(patientid)) {
    group1 <- paste(patientid[i], "_Recurrent", sep = "")
    group2 <- paste(patientid[i], "_Primary", sep = "")

    markers <- FindMarkers(integrated,
        ident.1 = group1,
        ident.2 = group2,
        group.by = "Patient_tumortype",
        test.use = "MAST",
        min.pct = 0,
        logfc.threshold = 0
    )

    assign(paste(patientid[i], "markers_MAST_predicted_rna", sep = "_"), markers)

    write.table(markers,
        paste(savedir, "Table/", patientid[i], "_predicted_rna_recurrent_vs_primary_markers_MAST_nocutoff.txt", sep = ""),
        quote = F,
        row.names = T,
        col.names = T,
        sep = "\t"
    )
}

### Checking the gene expression for the transcription factors
c("SP9","SP8","SP4","SP3","SP2","SP1",
  "KLF6","KLF5","KLF3","KLF2","KLF16","KLF15","KLF11",
  "EGR4","EGR3","EGR2","EGR1","FOXN3","FOXO3","FOXG1","FOXD1",
  "FOXP2","FOXA2","FOXA3","FOXK1","FOXP1","FOXC2","FOXK2","FOXP3",
  "FOXO6","FOXO4","FOXL1","FOXI1","FOXF2","HIF1A") -> genes

p <- DotPlot(integrated, genes, assay = "imputed_rna", group.by = "celltypes") +
    scale_color_gradientn(colors = c("royalblue2","royalblue1","lightblue3","lightblue2","lightblue1","lightpink2","indianred1","indianred3","indianred4","firebrick4")) +
    # scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/SP_KLF_EGR_FOX_Genes_imputed_RNA_celltypes.pdf"), height = 14, width = 9)
print(p)
dev.off()

# HIF1A
HIF1A_link_peaks <- as.data.frame(integrated@assays$ATAC@links[grep("HIF1A",integrated@assays$ATAC@links$gene),])

### The motif is already present in the object
motif_names <- as.data.frame(as.matrix(integrated@assays$ATAC@motifs@motif.names))
stopifnot(all(rownames(motif_names) == colnames(motif_df)))
colnames(motif_df) <- paste(motif_names$V1, rownames(motif_names), sep = "_")

HIF1A_required_df <- motif_df[match(HIF1A_link_peaks$peak, rownames(motif_df)), ]

as.matrix(HIF1A_required_df) -> HIF1A_required_mat

# Convert TRUE to 1 and FALSE to 0
HIF1A_required_df_num <- matrix(as.integer(HIF1A_required_mat), nrow = nrow(HIF1A_required_mat), ncol = ncol(HIF1A_required_mat))
rownames(HIF1A_required_df_num) <- rownames(HIF1A_required_mat)
colnames(HIF1A_required_df_num) <- colnames(HIF1A_required_mat)

HIF1A_columnSum <- as.data.frame(colSums(HIF1A_required_df_num))
HIF1A_columnSum$TF <- gsub("_.*.","",rownames(HIF1A_columnSum))
colnames(HIF1A_columnSum) <- c("peak_number","TF")

HIF1A_required_df_num_ZF <- HIF1A_required_df_num[,grep("SP[0-9]|KLF|FOX|EGR|HIF1A",colnames(HIF1A_required_df_num))]
colnames(HIF1A_required_df_num_ZF) <- gsub("_.*.","",colnames(HIF1A_required_df_num_ZF))
HIF1A_required_df_num_ZF_ordered <- HIF1A_required_df_num_ZF[,order(colnames(HIF1A_required_df_num_ZF))]
HIF1A_required_df_num_ZF_ordered_2 <- HIF1A_required_df_num_ZF_ordered[,-1] ### Removing Arnh:HIF1A

library(ComplexHeatmap)
p <- Heatmap(t(as.matrix(HIF1A_required_df_num_ZF_ordered_2)),
        name = "Peak Number",
        col = colorRampPalette(c("white", "mediumturquoise"))(50),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title_position = "topleft", legend_direction = "horizontal"))

pdf("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/motif/HIF1A_binding.pdf", height = 12, width = 7)
draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

### Making an heatmap for the Zinc Finger Pseudo Count gene expression
### Malignant cells
cellnames <- rownames(scATAC@meta.data[grep("Endo-Pericytes|Oligo_OPC|Astrocytes|Immune|Neuron", scATAC@meta.data$celltypes, invert = TRUE), ])
scATAC_mal <- subset(scATAC, cells = cellnames)

genes <- c("EGR1","EGR2","EGR3","EGR4",
            "KLF11","KLF15","KLF16","KLF2","KLF3","KLF5","KLF6",
            "SP1","SP2","SP3","SP4","SP8","SP9","HIF1A")

genes <- c("VEGFA","RICTOR","AKT3","PIK3CA","FOXO1","CCND2",
            "EIF4EBP1","EIF2AK3","EIF2A","EIF3E","EIF2AK4","XBP1","ATF6","HSPB1", "SPP1")

p <- DotPlot(scATAC_mal, genes, assay = "imputed_rna", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/snRNA_recurrent_high_pseudoCount_RNA_celltypes.pdf"), height = 7, width = 7)
print(p)
dev.off()

genes <- c("VEGFA","RICTOR","AKT3","PIK3CA","FOXO1","CCND2",
            "EIF4EBP1","EIF2AK3","EIF2A","EIF3E","EIF2AK4","XBP1","ATF6","HSPB1", "SPP1")

p <- DotPlot(scATAC_mal, genes, assay = "imputed_rna", group.by = "celltypes") +
    scale_color_gradientn(colors = ArchRPalettes$solarExtra) +
    coord_flip() +
    scale_size(range = c(1, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(paste0(savedir, "/dotplot/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "/dotplot/snRNA_recurrent_high_pseudoCount_RNA_celltypes.pdf"), height = 7, width = 7)
print(p)
dev.off()

### Performing GRN
library(Pando)
DefaultAssay(scATAC_mal) <- "imputed_rna"
scATAC_mal <- FindVariableFeatures(scATAC_mal)
DefaultAssay(scATAC_mal) <- "ATAC"
scATAC_mal <- initiate_grn(scATAC_mal, peak_assay = "ATAC", rna_assay = "imputed_rna")
library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs)

seurat_object <- find_motifs(
    scATAC_mal,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# This uses motifmatchr to pair up TFs with their putative binding sites. Pando provides a custom motif database (motifs) compiled from JASPAR and CIS-BP, but
# in principle any PFMatrixList object can be provided here. A data frame with motif-to-TF assignments can be provided in the motif_tfs argument.

# Inferring the GRN

# Now everything should be ready to infer the GRN by fitting regression models for the expression of each gene. In Pando, this can be done by using the 
# function infer_grn():

seurat_object <- infer_grn(
    seurat_object,
    peak_to_gene_method = 'Signac',
    method = 'glm'
)

#### Identifying the difference b/w GCP and GN in primary and recurrent
scATAC@meta.data$broad_CT <- scATAC@meta.data$predicted.id

 patterns <- c("Oligo","GCP", "GCP_cycling_1","Astro","Pericyte","OPC","Neuron_other","GN_Premigratory",
 "GCP_cycling_2","GCP_HSP","GN_postmigratory","GN_migratory","Immune","RL_like","Endo","GN_cycling")

 replacements <- c("non-malign","GCP","GCP","non-malign","non-malign","non-malign","non-malign","GN",
 "GCP","GCP","GN","GN","non-malign","non-malign","non-malign","GN")

for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  scATAC@meta.data$broad_CT <- str_replace_all(scATAC@meta.data$broad_CT, pattern, replacements[i])
}

#### Longitudinal analysis
celltypes_ind <- table(scATAC@meta.data$Sample,scATAC@meta.data$broad_CT)
write.table(celltypes_ind, paste("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/celltype_broad_individual.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/Table/celltype_broad_individual.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
  df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample","celltype","percentage")
library(ggplot2)

library(ggpubr)
library(rstatix)
library(stringr)

long_sample <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/Table/sample_longitude.txt", header= TRUE,sep="\t")

patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Tumor.type
# names(replacement) <- patterns
df_melted$timepts <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$timepts <- str_replace_all(df_melted$timepts, pattern, replacements[i])
}

### Patient ID
patterns = long_sample$Sample.ID
patterns <- gsub("SF7994R","SF7994",patterns)
patterns<- gsub("-","_",patterns)
replacements = long_sample$Patient.ID
# names(replacement) <- patterns
df_melted$Patient.ID <- df_melted$sample
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  df_melted$Patient.ID <- str_replace_all(df_melted$Patient.ID, pattern, replacements[i])
}

# df_melted$celltype2 <- gsub("GCP_cycling","GCP",df_melted$celltype) %>% gsub("GNearly|GNlate","GN",.)
df_melted_paired <- df_melted[grep("7316_4529|7316_5881|7316_278|7316_333|7316_2118|7316_737|7316_3023|7316_931",df_melted$sample),]
df_melted_paired$Patient_Timepts <- paste(df_melted_paired$Patient.ID, df_melted_paired$timepts, sep = "_")
# df_wide <- dcast(df_melted_paired, celltype ~ Patient_Timepts, value.var = "percentage", fun.aggregate=mean)
# df_melted -> df_melted_paired

stat.test <- df_melted_paired %>%
  group_by(celltype) %>%
  t_test(percentage ~ timepts, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test$p <- (stat.test$p)/2

bxp <- ggboxplot(
  df_melted_paired, x = "celltype", y = "percentage", 
  color = "timepts", palette = c("#00AFBB", "#E7B800")
) 

bxp2 <- bxp +   geom_dotplot(
  aes(fill = timepts, color = timepts), trim = FALSE,
  binaxis='y', stackdir='center', dotsize = 0.15,
  position = position_dodge(0.8)
)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +   theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white"))

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

savedir <- "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/"
pdf(paste(savedir,"Table/t_test_paired_primary_recurrent_samples_bxplot_two_sided_celltype_all.pdf",sep = ""), width = 10, height = 5)
bxp4
dev.off()


