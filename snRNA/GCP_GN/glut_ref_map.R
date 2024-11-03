#### Mapping the patient object on the reference data
library(Seurat)
library(Matrix)
library(dplyr)
reference_subset <- readRDS("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/human_cerebellar_glut_subset.RDS")
maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
filepath = list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", recursive = TRUE, pattern = ".RDS", full.names = TRUE)
objname = gsub(".RDS","",basename(filepath))

for (i in 1:length(filepath)) {
  try({
    object <- readRDS(filepath[i])
    DefaultAssay(object) <- "SCT"

    anchors <- FindTransferAnchors(reference = reference_subset,
                                query = object,
                                normalization.method = "SCT",
                                dims = 1:50)
    
    predictions.assay <- TransferData(anchorset = anchors,
                           refdata = reference_subset$figure_clusters,
                           prediction.assay = TRUE,
                           weight.reduction = object[["pca"]],
                           dims = 1:30, k.weight = 50)
    
    object[["predictions"]] <- predictions.assay

    celltypes <- unique(reference_subset@meta.data$figure_clusters)
    
    DefaultAssay(object) <- "predictions"
    pdf(paste(maindir,objname[i],"/featureplot/",objname[i],"_reference_glut_mapped.pdf",sep = ""), width = 10, height = 9)
    FeaturePlot(object, features = c(celltypes), ncol = 3)
    dev.off()
    
    prediction <- TransferData(anchorset = anchors,
                              refdata = reference_subset$figure_clusters,
                              prediction.assay = FALSE,
                              weight.reduction = object[["pca"]],
                              dims = 1:30)
                              
    object <- AddMetaData(object, metadata = prediction)
    
    p <- DimPlot(object,
                 group.by = "predicted.id", reduction = "umap",
                 cols = c("10-Glia"='#1f77b4',"02-RL"='#ff7f0e',"03-GCP"='#279e68',"01-PC"='#d62728',"15-Meninges"='#aa40fc',
                          "05-eCN/UBC"='#8c564b',"08-BG"='#e377c2',"09-Ast"='#b5bd61',"19-Ast/Ependymal"='#b5bd61',"11-OPC"='#17becf',
                          "04-GN"='#aec7e8',"07-PIP"='#ffbb78',"14-Microglia"='cadetblue3',"06-iCN"='cornsilk4',"20-Choroid"='plum2',
                          "18-MLI"='yellow',"13-Endothelial"='hotpink',"16-Pericytes" = "black","12-Committed OPC" = "cyan"))
    
    dir.create(paste(maindir,objname[i],"/UMAP/",sep = ""),showWarnings = FALSE)
    pdf(paste(maindir,objname[i],"/UMAP/",objname[i],"_celltypes_umap_glut.pdf",sep = ""))
    print(p)
    dev.off()

    ### Imputing the data
    DefaultAssay(object) <- "RNA"
    library(reticulate)
    library(Rmagic)
    use_python("/diazlab/data3/.abhinav/tools/miniconda3/envs/py39/bin/python")
    py_discover_config("magic") # to check
    object <- magic(object, npca=20) ## imputing the RNA data as for RNA PCs are 20
    DefaultAssay(object) <- "MAGIC_RNA"
    object <- ScaleData(object,features=rownames(object)) # Sacling data on all the features as we need it during the feature plot scaled data for gene that are not variable

    saveRDS(object, paste(maindir,objname[i],"/saveRDS_obj/",objname[i],"_refglut_mapped_imputed.rds",sep = ""))
    })
}
