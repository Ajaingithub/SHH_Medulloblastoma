
# Load necessary libraries
library(parallel)
library(doParallel)
library(dplyr)
library(Seurat)
library(infercnv)

# Detect the number of cores
# num_cores <- detectCores()
# cat("Number of cores detected:", num_cores, "
")

# Register the parallel backend
# cl <- makeCluster(10)
# registerDoParallel(cl)

count <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_7316_1666/CNV/obj_7316_1666_raw_counts.txt", header = TRUE, row.names = 1)
gene_order <- read.table("/diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/genes/final_gene_position_only_chr.txt", header = FALSE, row.names = 1)
annot_file <- read.table("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_7316_1666/CNV/obj_7316_1666_tum_norm_annotation_marker.txt", header = FALSE, row.names = 1)
rownames(annot_file) <- paste("X",gsub("-", ".", row.names(annot_file)),sep="")
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=count,
                                               gene_order_file=gene_order,
                                               annotations_file=annot_file,
                                               ref_group_names=c("normal"))

# Perform infercnv operations to reveal CNV signal
dir.create("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_7316_1666/CNV/output_marker",showWarnings=FALSE,recursive = TRUE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_7316_1666/CNV/output_marker",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             BayesMaxPNormal=0.4,
                             HMM=T,
                             no_prelim_plot=TRUE,
                             png_res=300)

saveRDS(infercnv_obj, "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/obj_7316_1666/CNV/output_marker/obj_7316_1666_infercnv.RDS")

# Stop the cluster
# stopCluster(cl)

