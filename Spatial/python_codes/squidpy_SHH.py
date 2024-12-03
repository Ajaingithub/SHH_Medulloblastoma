#### Squidpy 
# https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html
### Image Feature Extraction ####

# Visium datasets contain high-resolution images of the tissue that was used for the gene extraction. Using the function squidpy.im.calculate_image_features()
# you can calculate image features for each Visium spot and create a obs x features matrix in adata that can then be analyzed together with the obs x 
# gene gene expression matrix.
# 
# By extracting image features we are aiming to get both similar and complementary information to the gene expression values. Similar information is 
# for example present in the case of a tissue with two different cell types whose morphology is different. Such cell type information is then contained 
# in both the gene expression values and the tissue image features.
# 
# Squidpy contains several feature extractors and a flexible pipeline of calculating features of different scales and sizes. There are several detailed 
# examples of how to use squidpy.im.calculate_image_features(). Extract image features provides a good starting point for learning more.
# 
# Here, we will extract summary features at different crop sizes and scales to allow the calculation of multi-scale features and segmentation features. 
# For more information on the summary features, also refer to Extract summary features.

# C4 server UCSF
# conda activate squidpy
# python3

import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq
import imagecodecs

savedir = "/diazlab/data3/.abhinav/projects/SHH/Visium/scanpy/"
SF10210 = sc.read(savedir + "SHH_SF10210.h5ad")

#### Image Feature does not really worked since the image quality is not really good. 
img=sq.im.ImageContainer("/diazlab/data2/bhyu0217/projects/DoD/visium/image_files/V11T09-094_A1_SF10210/V11T09-094_A1_SF10210.jpg")
# sq.pl.spatial_scatter(SF10210, color="clusters", img_alpha=0, size=1.3)
# plt.savefig(savedir + "SHH_SF10210/spatial_cluster_scatter.png")

# calculate features for different scales (higher value means more context)
for scale in [1.0, 2.0]:
    feature_name = f"features_summary_scale{scale}"
    sq.im.calculate_image_features(
        SF10210,
        img.compute(),
        features="summary",
        key_added=feature_name,
        n_jobs=4,
        scale=scale,
    )

# combine features in one dataframe
SF10210.obsm["features"] = pd.concat(
  [SF10210.obsm[f] for f in SF10210.obsm.keys() if "features_summary" in f],
  axis="columns",
)

# make sure that we have no duplicated feature names in the combined table
SF10210.obsm["features"].columns = ad.utils.make_index_unique(
    SF10210.obsm["features"].columns
)

# We can use the extracted image features to compute a new cluster annotation. 
# This could be useful to gain insights in similarities across spots based on image morphology.
# helper function returning a clustering
def cluster_features(features: pd.DataFrame, like=None, res = 1) -> pd.Series:
  """
  Calculate leiden clustering of features.
  Specify filter of features using `like`.
  """

  import anndata as ad
  
  # filter features
  if like is not None:
    features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution = res)
    
    return adata.obs["leiden"]

# calculate feature clusters
SF10210.obs["features_cluster"] = cluster_features(SF10210.obsm["features"], like="summary")

# compare feature and gene clusters
# fig, ax = plt.subplots(figsize=(10, 6))
# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
sq.pl.spatial_scatter(SF10210, color=["features_cluster"], size=1.3)
plt.tight_layout()

plt.savefig(savedir + "SHH_SF10210/feature_extraction.png")

# Spatial statistics and graph analysis
#### Neighborhood enrichment ####

# Computing a neighborhood enrichment can help us identify spots clusters that share a common neighborhood structure across the tissue. We can compute 
# such score with the following function: squidpy.gr.nhood_enrichment(). In short, it’s an enrichment score on spatial proximity of clusters: if spots
# belonging to two different clusters are often close to each other, then they will have a high score and can be defined as being enriched. On the other
# hand, if they are far apart, and therefore are seldom a neighborhood, the score will be low and they can be defined as depleted. This score is based
# on a permutation-based test, and you can set the number of permutations with the n_perms argument (default is 1000).

# Since we have the cluster information with the region information. We can use the neighbour enrichment for each of these regions.
# in R
# readRDS("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/R/SF10210/saveRDS_obj/SF10210_impute_2.RDS") -> SF10210
# celltype_cellid <- data.frame(rownames(SF10210@meta.data),SF10210@meta.data$predicted.id) ## 
# colnames(celltype_cellid) <- c("cell_barcode","Predicted_id")
# write.table(celltype_cellid, 
#             "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/R/SF10210/Table/celltype_cellid.txt", 
#             sep = "\t", quote=F, row.names =F)

celltype_clus = pd.read_csv("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/R/SF10210/Table/celltype_cellid.txt", delimiter = "\t")
celltype_clus.index = celltype_clus["cell_barcode"]

### Now we will merge the dataframe
SF10210.obs=pd.merge(SF10210.obs, celltype_clus, left_index=True, right_index=True, how='left')

### We have some with NA due to some quality removal of spots in seurat and scanpy
from collections import Counter
Counter(SF10210.obs.loc[:,"Predicted_id"])

### based on this we have highest of the Media_2 so we will fill the na with the Media_2
SF10210.obs.loc[:,"Predicted_id"]=SF10210.obs.loc[:,"Predicted_id"].fillna("GCP")

CT_colors  = {"GCP_cycling":"#084b62", "GCP_ribo" : "#258fb3",
        "GN_ribo" : "#136e43","GN_premigratory" : "#5fe3a7",
        "GCP_stress" : "#00bfff", "GN_migratory" : "#30b377", 
        "GN_postmigratory" : "#107044", "GCP" : "#6fb9ed",
        "Pericyte" : "#aa40fc", "Astro" : "#8c564b", 
        "Immune" : "#000000", "Oligo" : "#b5bd61", 
        "OPC" : "#FF69B4", "Endo" : "#5C4033", 
        "GCP_GN" : "#00FFFF", "Neuron" : "#0966a9"}

cluster_key = "Predicted_id"
SF10210.uns[f"{cluster_key}_colors"] = [CT_colors[cluster] for cluster in SF10210.obs[cluster_key].cat.categories]

### To check if nan is removed
Counter(SF10210.obs.loc[:,"Predicted_id"])

plt.rcParams["figure.figsize"] = (10, 8)
sc.pl.spatial(SF10210, img_key="hires", color=["Predicted_id"], size=1.5)
plt.tight_layout()
savedir = "/diazlab/data3/.abhinav/projects/SHH/Visium/scanpy/"
plt.savefig(savedir + "SHH_SF10210/spatial_celltypes.pdf", format = "pdf")

sq.gr.spatial_neighbors(SF10210)
sq.gr.nhood_enrichment(SF10210, cluster_key="Predicted_id")
sq.pl.nhood_enrichment(SF10210, cluster_key="Predicted_id")

plt.savefig(savedir + "SHH_SF10210/neigborhood_enrichment_predicted.pdf",bbox_inches="tight",format="pdf")

sq.gr.co_occurrence(SF10210, cluster_key="Predicted_id")
sq.pl.co_occurrence(
    SF10210,
    cluster_key="Predicted_id"
    # clusters=["GCP_stress","GCP_cycling"]
)

plt.savefig(savedir + "SHH_SF10210/neigborhood_enrichment_predicted_all_celltypes.pdf",bbox_inches="tight",format="pdf")

#### Automating the script for all the samples
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import squidpy as sq
from collections import Counter

# Directories
savedir = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/"
image_base_dir = "/diazlab/data2/bhyu0217/projects/DoD/visium/image_files/"
celltype_file_base_dir = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/R/{}/Table/"

# Sample list
all_dir = os.listdir(image_base_dir)
# samples = ["SHH_SF08539", "SHH_SF07994",
#            "SHH_SF03521", "SHH_SF12434", "SHH_SF09782",
#            "SHH_SF08539.2", "SHH_SF09782.2"]

# samples = ['SHH_SF08539.2', 'SHH_SF09782.2']
samples = ['SHH_SF10210']
sample_ids = [sample.split('_')[-1] for sample in samples]

matched_dirs = [dir for sample_id in sample_ids for dir in all_dir if sample_id in dir]
# matched_dirs = ["V11T09-036_A1_SF08539","V11T09-036_D1_SF09782"]

### adding more samples in the list 
matched_dirs.append('V11T09-036_A1_SF08539')
matched_dirs.append('V11T09-036_D1_SF09782')

CT_colors  = {"GCP_cycling":"#084b62", "GCP_ribo" : "#258fb3",
        "GN_ribo" : "#136e43","GN_premigratory" : "#5fe3a7",
        "GCP_stress" : "#00bfff", "GN_migratory" : "#30b377", 
        "GN_postmigratory" : "#107044", "GCP" : "#6fb9ed",
        "Pericyte" : "#aa40fc", "Astro" : "#8c564b", 
        "Immune" : "#000000", "Oligo" : "#b5bd61", 
        "OPC" : "#FF69B4", "Endo" : "#5C4033", 
        "GCP_GN" : "#00FFFF", "Neuron" : "#0966a9"}

import pandas as pd
# Helper function for clustering features
# We can use the extracted image features to compute a new cluster annotation. 
# This could be useful to gain insights in similarities across spots based on image morphology.
# helper function returning a clustering
def cluster_features(features: pd.DataFrame, like=None, res = 1) -> pd.Series:
  """
  Calculate leiden clustering of features.
  Specify filter of features using `like`.
  """
  import anndata as ad  
  # filter features
  if like is not None:
    features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution = res)
    return adata.obs["leiden"]

os.chdir("/diazlab/data3/.abhinav/projects/SHH/Visium/")
import importlib  # Standard library module
import sq_process_sample

# Reload the module after making changes
importlib.reload(sq_process_sample)

# Main loop to process all samples
for i in range(len(samples)):
    sample = samples[i]
    matched_dir = matched_dirs[i]
    sq_process_sample.process_sample(sample,matched_dir,savedir,image_base_dir,celltype_file_base_dir,cluster_features,CT_colors)

### Ligand Receptor Interactions 
# We are continuing the analysis showing couple of feature-level methods that are very relevant for the analysis of spatial molecular data. For instance, after 
# quantification of cluster co-occurrence, we might be interested in finding molecular instances that could potentially drive cellular communication. This naturally
# translates in a ligand-receptor interaction analysis. In Squidpy, we provide a fast re-implementation the popular method CellPhoneDB [Efremova et al., 2020] 
# (code <https://github.com/Teichlab/cellphonedb>_ ) and extended its database of annotated ligand-receptor interaction pairs with the popular database Omnipath 
# [Türei et al., 2016]. You can run the analysis for all clusters pairs, and all genes (in seconds, without leaving this notebook), with squidpy.gr.ligrec(). 
# Furthermore, we’ll directly visualize the results, filtering out lowly-expressed genes (with the means_range argument) and increasing the threshold for the adjusted
# p-value (with the alpha argument).

sq.gr.ligrec(SF10210,
n_perms=1000,
cluster_key="Predicted_id",
use_raw=False,
threshold=0.001,
transmitter_params={"categories": "ligand"},
receiver_params={"categories": "receptor"},)

sq.pl.ligrec(SF10210,
cluster_key="Predicted_id",
means_range=(0.5, np.inf),
source_groups=["GCP_ribo"],
pvalue_threshold=0.001,
target_groups=["Endo"],
# target_groups=["Immature_GAM"],
# alpha=0.05,
swap_axes=False)

plt.savefig(savedir + "SHH_SF10210/LR_cellphoneDB_celltypes_GCP_ribo_Endo.pdf",bbox_inches="tight", format = "pdf")
# plt.savefig(savedir + "slide_075_D1/LR_cellphoneDB_celltypes_per_1000_clus4_Media_Pure_GAM.pdf",bbox_inches="tight", format = "pdf")

sq.pl.ligrec(slide_075_D1,
cluster_key="celltypes",
means_range=(0.5, np.inf), 
# source_groups="clus4_Media_Pure_GAM",
source_groups=["clus4_Media_Pure_GAM","clus8_Media_Mixed_Mac_dominated","clus0_Media_Mixed_SMC_dominated","clus3_Media_Pure_SMC_intact_1",
"clus5_Media_Mixed_GAM_SMC","clus9_Media_Pure_SMC_intact_2","clus7_Media_Mixed_Fibrosis_Healed"],
# target_groups=["Immature_GAM"],
alpha=0.05,
swap_axes=True)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Yuki/spatial/Run01/downstream/python_scanpy/"
plt.savefig(savedir + "slide_075_D1/slide_075_D1_LR_cellphoneDB_celltypes_per_1000_whole_media.pdf",bbox_inches="tight", format = "pdf")


### to get the loose interaction
sq.gr.ligrec(slide_075_D1,
n_perms=100, 
cluster_key="celltypes",
use_raw=False,
threshold=0.05,
transmitter_params={"categories": "ligand"},
receiver_params={"categories": "receptor"},)

sq.pl.ligrec(slide_075_D1,
cluster_key="celltypes",
means_range=(0.5, np.inf), 
# source_groups="clus4_Media_Pure_GAM",
# target_groups=["clus8_Media_Mixed_Mac_dominated","clus0_Media_Mixed_SMC_dominated","clus3_Media_Pure_SMC_intact_1",
# "clus5_Media_Mixed_GAM_SMC","clus9_Media_Pure_SMC_intact_2","clus7_Media_Mixed_Fibrosis_Healed"],
# target_groups=["Immature_GAM"],
alpha=0.05,
swap_axes=False)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Yuki/spatial/Run01/downstream/python_scanpy/"
plt.savefig(savedir + "slide_075_D1/slide_075_D1_LR_cellphoneDB_loose_interactions.pdf",bbox_inches="tight", format = "pdf")

slide_075_D1.write(savedir + "slide_075_D1/adata_save.h5ad")

#### slide_075_A1 neigborhood 
celltype_clus_subset=celltype_clus[celltype_clus.index.str.contains(list(samplename)[2])]
celltype_clus_subset.loc[:,"barcode"] = celltype_clus_subset.index.str.split("_").str[-1].tolist()
celltype_clus_subset.index = celltype_clus_subset.loc[:,"barcode"]

### Now we will merge the dataframe
slide_075_A1.obs=pd.merge(slide_075_A1.obs, celltype_clus_subset, left_index=True, right_index=True, how='left')

### We have some with NA due to some quality removal of spots in seurat and scanpy
from collections import Counter
Counter(slide_075_A1.obs.loc[:,"celltypes_y"])

### based on this we have highest of the Media_2 so we will fill the na with the Media_2
slide_075_A1.obs.loc[:,"celltypes_y"]=slide_075_A1.obs.loc[:,"celltypes_y"].fillna("clus3_Media_Pure_SMC_intact_1")

plt.rcParams["figure.figsize"] = (10, 8)
sc.pl.spatial(slide_075_A1, img_key="hires", color=["celltypes_y"], bw=False, alpha_img=0, size=1.5)
plt.tight_layout()
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Yuki/spatial/Run01/downstream/python_scanpy/"
plt.savefig(savedir + list(samplename)[1] + "/spatial_celltypes_y.png")
plt.show()

sq.gr.spatial_neighbors(slide_075_A1)
sq.gr.nhood_enrichment(slide_075_A1, cluster_key="celltypes_y")

sq.pl.nhood_enrichment(slide_075_A1, cluster_key="celltypes_y")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Yuki/spatial/Run01/downstream/python_scanpy/"
plt.savefig(savedir + "slide_075_A1/neigborhood_enrichment_seurat_combine_celltypes_y.pdf",bbox_inches="tight", format = "pdf")

sq.gr.ligrec(slide_075_A1, 
n_perms=100, 
cluster_key="celltypes_y",
use_raw=False,
threshold=0.001,
transmitter_params={"categories": "ligand"},
receiver_params={"categories": "receptor"},)

sq.pl.ligrec(slide_075_A1,
cluster_key="celltypes_y",
means_range=(0.5, np.inf), 
source_groups="clus4_Media_Pure_GAM",
target_groups=["clus8_Media_Mixed_Mac_dominated","clus0_Media_Mixed_SMC_dominated","clus3_Media_Pure_SMC_intact_1",
"clus5_Media_Mixed_GAM_SMC","clus9_Media_Pure_SMC_intact_2","clus7_Media_Mixed_Fibrosis_Healed"],
# target_groups=["Immature_GAM"],
alpha=0.05,
swap_axes=True)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Yuki/spatial/Run01/downstream/python_scanpy/"
plt.savefig(savedir + "slide_075_A1/LR_cellphoneDB_celltypes_y_per_1000_clus4_Media_Pure_GAM.pdf",bbox_inches="tight", format = "pdf")
