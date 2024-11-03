library(Seurat)
library(reticulate)
library(Rmagic)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
py_discover_config("magic") # to check

integrated <- readRDS("/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_with_RNA.RDS")

DefaultAssay(integrated) <- "RNA"
integrated <- magic(integrated, npca = 30) ## imputing the RNA data as for RNA PCs are 20

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

saveRDS(integrated, "/diazlab/data3/.abhinav/projects/SHH/scATAC/merging/saveRDS_obj/integrated_with_RNA_with_imputation_and_motifs.RDS")
