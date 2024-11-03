#### Performing reference mapping
### Unbiased
library(Seurat)
library(dplyr)
reference <- readRDS("/diazlab/data3/.abhinav/resources/reference/human_cerebellar/human_cerebellar_SCT.RDS")
maindir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/"
filepath = list.files("/diazlab/data3/.abhinav/projects/SHH/snRNA/subset/", recursive = TRUE, pattern = ".RDS", full.names = TRUE)
objname = gsub(".RDS","",basename(filepath))

for (i in 1:length(filepath)) {
  try({
    object <- readRDS(filepath[i])
    DefaultAssay(object) <- "SCT"

    anchors <- FindTransferAnchors(reference = reference,
                                query = object,
                                normalization.method = "SCT",
                                dims = 1:50)
    
    predictions.assay <- TransferData(anchorset = anchors,
                           refdata = reference$figure_clusters,
                           prediction.assay = TRUE,
                           weight.reduction = object[["pca"]],
                           dims = 1:30, k.weight = 50)
    
    object[["predictions"]] <- predictions.assay

    celltypes <- unique(reference@meta.data$figure_clusters)
    
    DefaultAssay(object) <- "predictions"
    pdf(paste(maindir,objname[i],"/featureplot/",objname[i],"_refall_mapped.pdf",sep = ""), width = 10, height = 9)
    print(FeaturePlot(object, features = c(celltypes), ncol = 3))
    dev.off()
    
    prediction <- TransferData(anchorset = anchors,
                              refdata = reference$figure_clusters,
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
    pdf(paste(maindir,objname[i],"/UMAP/",objname[i],"_refall_celltypes_umap.pdf",sep = ""))
    print(p)
    dev.off()


    assign(paste(objname[i],"_CT",sep = ""), object)
    saveRDS(object, paste(maindir,objname[i],"/saveRDS_obj/",objname[i],"_refall_mapped.rds",sep = ""))
    })
}
