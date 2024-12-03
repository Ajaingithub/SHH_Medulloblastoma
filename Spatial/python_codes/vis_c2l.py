## This is the tutorial https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html#Loading-packages
# Cell2location is a principled Bayesian model that estimates which combination of cell types in which cell abundance could have given 
# the mRNA counts in the spatial data, while modelling technical effects (platform/technology effect, contaminating RNA, unexplained variance).
# 
## Loading packages
# All the codes are running the conda environment py39
# conda activate spatial_analysis

### Since we need the snRNA reference for cell2location. We are converting it into c2l
"""
### R command
library(Seurat)
library(anndata)
snRNA <- readRDS("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/shh_integrated_RNA_NN_cluster_0.8_imputed.RDS")
snRNA@meta.data$cell_types <- snRNA@meta.data$sub_celltypes
snRNA@meta.data$cell_types <- gsub("GCP_cycling_1", "GCP_cycling", snRNA@meta.data$cell_types) %>%
    gsub("GCP_cycling_2", "GCP_ribo", .) %>%
    gsub("GN_cycling", "GN_ribo", .) %>%
    gsub("GCP_HSP", "GCP_stress", .) %>%
    gsub("Neuron_other", "Neuron", .) %>%
    gsub("RL_like", "GCP_GN", .) %>%
    gsub("GN_Premigratory", "GN_premigratory", .)

# Save as .h5ad file
DefaultAssay(snRNA) <- "RNA"
snRNA_diet <- DietSeurat(snRNA, assay = c("RNA","integrated"), dimreducs = c("pca","umap"), layers = c("counts","data"))
SaveH5Seurat(snRNA_diet, filename = "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/snRNA_diet.h5Seurat", overwrite = TRUE)
Convert("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/snRNA_diet.h5Seurat", dest = "h5ad", overwrite = TRUE)
"""

import scanpy as sc
import squidpy as sq
import cell2location as c2l
import numpy as np
import sys
import anndata
import pandas as pd
import os
import gc
import matplotlib.pyplot as plt
import matplotlib as mpl# Set the default figure size
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

#### Checking if the 

# First, let’s define where we save the results of our analysis
# create paths and names to results folders for reference regression and cell2location models

results_folder = '/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/SHH_SF10210/cell2location/'
ref_run_name = f'{results_folder}reference_signatures'
run_name = f'{results_folder}cell2location_map'

### Loading the spatial adata
SF10210 = sc.read("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/SHH_SF10210/saveh5ad/SHH_SF10210_save.h5ad")
SF10210.obs['sample'] = list(SF10210.uns['spatial'].keys())[0]
adata_vis = SF10210.copy()

# adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

# Checking the anndata visium-
pd.set_option('display.max_columns', None)
print(adata_vis.obs[:2])
print(adata_vis.var[:2])
# print(adata_vis.uns)
print(adata_vis.obsm)
print(adata_vis.var_names[:5])
print(adata_vis.obs_names[:2])
print(adata_vis.n_obs)
print(adata_vis.n_vars)
print(adata_vis.X.shape)
print(type(adata_vis.X))
print(pd.DataFrame(adata_vis.X.toarray()).iloc[1:5, :5])

adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True) # gene ids has been index and the column is dropped

# in place whether to change the original Dataframe
os.makedirs(f"{results_folder}", exist_ok=True)
os.chdir(results_folder)
sc.pl.spatial(adata_vis,color='VEGFA',gene_symbols='SYMBOL', save= 'VEGFA.pdf') ##saving the file

# Mitochondia-encoded genes (gene names start with prefix mt- or MT-) are irrelevant for spatial mapping because their expression represents technical artifacts 
# in the single cell and nucleus data rather than biological abundance of mitochondria. Yet these genes compose 15-40% of mRNA in each location. Hence, to avoid 
# mapping artifacts we strongly recommend removing mitochondrial genes.

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# Published scRNA-seq datasets of lymph nodes have typically lacked an adequate representation of germinal centre-associated immune cell populations due to age of patient donors. 
# We, therefore, include scRNA-seq datasets spanning lymph nodes, spleen and tonsils in our single-cell reference to ensure that we captured the full diversity of immune cell states 
# likely to exist in the spatial transcriptomic dataset.
# Here we download this dataset, import into anndata and change variable names to ENSEMBL gene identifiers.

# Read data
# adata_ref = sc.read(
#     f'./data/sc.h5ad',
#     backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
# )

adata_ref = sc.read('/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/snRNA_diet.h5ad')

# Checking the anndata visium
pd.set_option('display.max_columns', None)
print(adata_ref.obs[:2])
print(adata_ref.var[:2])
# print(adata_ref.uns)
print(adata_ref.obsm)
print(adata_ref.var_names[:5])
print(adata_ref.obs_names[:2])
print(adata_ref.n_obs)
print(adata_ref.n_vars)
print(adata_ref.X.shape)
print(type(adata_ref.X))
print(pd.DataFrame(adata_ref.X.toarray()).iloc[1:5, :5])

# sc.pl.umap(adata_ref,color='cell_types', save= 'snRNA_ref_celltypes.pdf')

# Here we rename genes to ENSEMBL ID for correct matching between single cell and spatial data.
adata_ref.var['SYMBOL'] = adata_ref.var.index
adata_ref.X = adata_ref.raw.X

### Since reference does not have ENSG ids
from pybiomart import Server
server = Server(host='http://www.ensembl.org')
dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

# Query the dataset for mapping gene names to Ensembl IDs
gene_mapping = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
import pandas as pd
gene_mapping_df = pd.DataFrame(gene_mapping)
mapping_dict = gene_mapping_df.set_index('Gene name')['Gene stable ID'].to_dict()

gene_names = adata_ref.var['SYMBOL'].tolist()
mapped_genes = [mapping_dict.get(gene, 'Not found') for gene in gene_names]
adata_ref.var['ENSG'] = mapped_genes
valid_genes_mask = adata_ref.var['ENSG'] != 'Not found'
adata_filtered = adata_ref[:, valid_genes_mask].copy()
adata_filtered.var.set_index('ENSG', drop=True, inplace=True)


# The default parameters cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12 are a good starting point, however, you can increase the cut-off to exclude more genes. To preserve marker genes of rare cell types we recommend low cell_count_cutoff=5, however, cell_percentage_cutoff2 and nonz_mean_cutoff can be increased to select between 8k-16k genes.
# In this 2D histogram, orange rectangle highlights genes excluded based on the combination of number of cells expressing that gene (Y-axis) and average RNA count for cells where the gene was detected (X-axis).
# In this case, the downloaded dataset was already filtered using this method, hence no density under the orange rectangle (to be changed in the future version of the tutorial).
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_filtered, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
# In Python, the .copy() method is used to create a shallow copy of an object. When you have a complex data structure like a list, dictionary, or in your case, possibly a NumPy array (given the indexing syntax), using .copy() ensures that you create a new object with the same values as the original but stored at a different memory location. This is particularly useful to avoid unintended side effects when modifying the copied object, as changes to the copied object do not affect the original one.
adata_filtered2 = adata_filtered[:, selected].copy()

# Next, we subset both datasets to the same gene set which is the baseline for the mapping between the single cell and spatial data.
shared_features = [
    feature for feature in adata_vis.var_names if feature in adata_filtered2.var_names
]
adata_filtered2 = adata_filtered2[:, shared_features].copy()
adata_vis = adata_vis[:, shared_features].copy()

# The signatures are estimated from scRNA-seq data, accounting for batch effect, using a Negative binomial regression model.
# Preparing anndata.
# First, prepare anndata object for the regression model:

# prepare anndata for the regression model
# cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
#                         # 10X reaction / sample / batch
#                         batch_key='Sample',
#                         # cell type, covariate used for constructing signatures
#                         labels_key='Subset',
#                         # multiplicative technical effects (platform, 3' vs 5', donor effect)
#                         categorical_covariate_keys=['Method']
#                        )

# import torch
# from torch.utils.data import Dataset
# import numpy as np

# class AnnDataDataset(Dataset):
#     def __init__(self, adata, batch_size):
#         self.adata = adata
#         self.batch_size = batch_size
#     def __len__(self):
#         return len(self.adata)
#     def __getitem__(self, index):
#         # Extract features and labels from AnnData object
#         # Modify this part based on your needs
#         features = self.adata.X[index]  # Ensure this is the feature matrix (adata.X)
#         labels = self.adata.obs['cell_types'].iloc[index]  # Example of extracting labels
#         # Convert to torch tensors
#         features = torch.tensor(features, dtype=torch.float32)
#         labels = torch.tensor(labels, dtype=torch.long)  # Assuming classification problem
#         return features, labels

# # Create the dataset
# dataset = AnnDataDataset(adata_filtered, batch_size=2500)

# # Initialize the DataLoader with more workers
# dataloader = DataLoader(dataset, batch_size=2500, num_workers=50, shuffle=True)

# # Check if DataLoader is working as expected
# for i, data in enumerate(dataloader):
#     print(f"Batch {i} loaded with shape {data.shape}")
#     if i > 2:  # Just check a few batches
#         break


c2l.models.RegressionModel.setup_anndata(
    adata=adata_filtered2,
    batch_key="SampleID",
    labels_key="cell_types",
    categorical_covariate_keys=["batch"],
)

# create the regression model
model = c2l.models.RegressionModel(adata_filtered2)

# view anndata_setup as a sanity check
model.view_anndata_setup()

# Check if GPU is available
use_gpu = torch.cuda.is_available()

trainer_kwargs = {
    'max_epochs': 250,
    'batch_size': 2500,
    'train_size': 1,
    'lr': 0.002,
    # 'gpus': 1 if use_gpu else 0,
   # 'num_workers': 50,  # Increase the number of workers for data loading
}
model.train(**trainer_kwargs)


# Training model.
# Now we train the model to estimate the reference cell type signatures.
# Note that to achieve convergence on your data (=to get stabilization of the loss) you may need to increase max_epochs=250 (See below).
# Also note that here we are using batch_size=2500 which is much larger than scvi-tools default and perform training on all cells in the data (train_size=1) - both parameters are defaults.
model.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, num_workers = 50)

# Determine if the model needs more training.
# Here, we plot Elbow loss history during training, removing first 20 epochs from the plot. This plot should have a decreasing trend and level off by the end of training. If it is still decreasing, increase max_epochs.

model.plot_history(20)
plt.savefig(f'{results_folder}figures/training_elbow_plot_2.pdf')
plt.close()

model.save("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/cell2location/")

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_filtered2 = model.export_posterior(
    adata_filtered2, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Save anndata object with results
# Rename _index column in adata.var
if '_index' in adata_filtered2.var.columns:
    adata_filtered2.var.rename(columns={'_index': 'index_col'}, inplace=True)

# If adata.raw exists, rename _index in raw.var as well
if adata_filtered2.raw is not None and '_index' in adata_filtered2.raw.var.columns:
    adata_filtered2.raw.var.rename(columns={'_index': 'index_col'}, inplace=True)

adata_file = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/cell2location/snRNA_filtered_adata_trained.h5ad"
adata_filtered2.write(adata_file)

## saving visium data also
adata_file2 = f"{results_folder}/SF10210_vis.h5ad"
adata_vis.write(adata_file2)

# You can compute the 5%, 50% and 95% quantiles of the posterior distribution directly rather than using 1000 samples from the distribution (or any other quantiles).
# This speeds up application on large datasets and requires less memory - however, posterior mean and standard deviation cannot be computed this way.
# adata_filtered2 = model.export_posterior(
#     adata_filtered2, use_quantiles=True,
#     # choose quantiles
#     add_to_varm=["q05","q50", "q95", "q0001"],
#     sample_kwargs={'batch_size': 2500}
# )

# Examine QC plots.
# Reconstruction accuracy to assess if there are any issues with inference. This 2D histogram plot should have most observations along a 
# noisy diagonal.
# The estimated expression signatures are distinct from mean expression in each cluster because of batch effects. For scRNA-seq datasets 
# which do not suffer from batch effect (this dataset does), cluster average expression can be used instead of estimating signatures with a
# model. When this plot is very different from a diagonal plot (e.g. very low values on Y-axis, density everywhere) it indicates problems
# with signature estimation.

# For saving the QC plots separately we have wrote the function plot_spatial.py in the same directory
# plt.close()
import matplotlib.pyplot as plt
model.plot_QC()
plt.savefig(f'{results_folder}figures/training_plot_QC.pdf', bbox_inches='tight')
plt.close()

# The model and output h5ad can be loaded later like this:
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_ref = sc.read_h5ad(adata_file)
# mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# Extracting reference cell types signatures as a pd.DataFrame.
# All parameters of the a Negative Binomial regression model are exported into reference anndata object, however for spatial mapping we just
# need the estimated expression of every gene in every cell type. Here we extract that from standard output:

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_filtered2.varm.keys():
    inf_aver = adata_filtered2.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_filtered2.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_filtered2.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_filtered2.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_filtered2.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# Cell2location: spatial mapping
# Find shared genes and prepare anndata. Subset both anndata and reference signatures:

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
adata_vis.X = adata_vis.layers["raw_counts"]

c2l.models.Cell2location.setup_anndata(adata=adata_vis, 
                                       batch_key="sample")

# Choosing hyperparameter ``N_cells_per_location``!
# It is useful to adapt the expected cell abundance N_cells_per_location to every tissue. This value can be estimated from paired histology
# images and as described in the note above. Change the value presented in this tutorial (N_cells_per_location=30) to the value observed 
# in your your tissue.

# Choosing hyperparameter ``detection_alpha``!
# To improve accuracy & sensitivity on datasets with large technical variability in RNA detection sensitivity within the slide/batch - you 
# need to relax regularisation of per-location normalisation (use detection_alpha=20). High technical variability in RNA detection sensitivity 
# is present in your sample when you observe the spatial distribution of total RNA count per location that doesn’t match expected cell numbers 
# based on histological examination.
# We initially opted for high regularisation (detection_alpha=200) as a default because the mouse brain & human lymph node datasets used in 
# our paper have low technical effects and using high regularisation strenght improves consistencly between total estimated cell abundance 
# per location and the nuclei count quantified from histology (Fig S8F in cell2location paper). However, in many collaborations, we see that 
# Visium experiments on human tissues suffer from technical effects. This motivates the new default value of detection_alpha=20 and the recommendation 
# of testing both settings on your data (detection_alpha=20 and detection_alpha=200).

# create and train the model
mod = c2l.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)

mod.view_anndata_setup()

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
plt.close()
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.savefig(f'{results_folder}/figures/spatial_training_elbow_plot.pdf')
plt.close()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)

# Save model
mod.save(f"{results_folder}/spatial_model/")

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{results_folder}SF10210_vis.h5ad"
adata_vis.write(adata_file)

## The model and output h5ad can be loaded late like this
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

plt.close()
mod.plot_QC()
plt.savefig(f'{results_folder}figures/spatial_plot_QC.pdf')
plt.close()

# When intergrating multiple spatial batches and when working with datasets that have substantial variation of detected RNA 
# within slides (that cannot be explained by high cellular density in the histology), it is important to assess whether cell2location 
# normalised those effects. You expect to see similar total cell abundance across batches but distinct RNA detection sensitivity (both 
# estimated by cell2location). You expect total cell abundance to mirror high cellular density in the histology.

mod.plot_spatial_QC_across_batches()
plt.savefig(f"{results_folder}/figures/spatial_plot_QC_across_batches.pdf")
plt.close()


# We use 5% quantile of the posterior distribution, representing the value of cell abundance that the model has high 
# confidence in (aka ‘at least this amount is present’).
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
# from cell2location.utils import select_slide
# slide = select_slide(adata_vis, 'SHH_SF10210')

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
                     
                     sc.pl.spatial(adata_vis, cmap='viridis',
                     # show first 8 cell types
                     color=['Astro', 'Endo', 'GCP', 'GCP_GN', 'GCP_cycling', 
                            'GCP_ribo', 'GCP_stress', 'GN_migratory', 
                            'GN_postmigratory', 'GN_premigratory', 'GN_ribo', 
                            'Immune', 'Neuron', 'OPC', 'Oligo', 'Pericyte'],
                     ncols=4, size=1.3,
                     img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                    vmin=0, vmax='p99.2'
                 )

plt.savefig(f'{results_folder}figures/spatial_each_cell_featureplots_viridis.pdf')
plt.close()

# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['Endo','Pericyte','GCP_stress']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels
# This is done to ensure that the column names (which might be used for labeling or indexing) are represented as strings.

# slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

with mpl.rc_context({'axes.facecolor':  'black','figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        # labels to show on a plot
        color=clust_col,
        # Color which are present in the slide.obs or adata_vis.obs
        labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=7,
        colorbar_position='right'
    )

plt.savefig(f'{results_folder}figures/non_malign_GCP_stress.pdf')
plt.close()

# Downstream analysis
# Identifying discrete tissue regions by Leiden clustering

# We identify tissue regions that differ in their cell composition by clustering locations using cell abundance estimated by cell2location.
# We find tissue regions by clustering Visium spots using estimated cell abundance each cell type. We constuct a K-nearest neigbour (KNN) 
# graph representing similarity of locations in estimated cell abundance and then apply Leiden clustering. The number of KNN neighbours 
# should be adapted to size of dataset and the size of anatomically defined regions (e.i. hippocampus regions are rather small compared 
# to size of the brain so could be masked by large n_neighbors). This can be done for a range KNN neighbours and Leiden clustering resolutions
# until a clustering matching the anatomical structure of the tissue is obtained.
# The clustering is done jointly across all Visium sections / batches, hence the region identities are directly comparable. When there are 
# strong technical effects between multiple batches (not the case here) sc.external.pp.bbknn can be in principle used to account for those
# effects during the KNN construction.

# compute KNN using the cell2location output stored in adata.obsm
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf', n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")

# We can use the location composition similarity graph to build a joint integrated UMAP representation of all section/Visium batches.
# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)

# show regions in UMAP coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['region_cluster'], size=15,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=10)

plt.savefig(f"{results_folder}/figures/spatial_UMAP_cluster_cellabundance.pdf")
plt.close()


with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20)
    
plt.savefig(f"{results_folder}/figures/spatial_UMAP_cluster_sample.pdf")
plt.close()

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata_vis, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5)
    
plt.savefig(f"{results_folder}/figures/spatial_cluster_high_res_imag.pdf", bbox_inches='tight')
plt.close()

# Identifying cellular compartments / tissue zones using matrix factorisation (NMF)
# Here, we use the cell2location mapping results to identify the spatial co-occurrence of cell types in order to better understand the 
# tissue organisation and predict cellular interactions. We performed non-negative matrix factorization (NMF) of the cell type abundance 
# estimates from cell2location.
# This NMF-based decomposition naturally accounts for the fact that multiple cell types and microenvironments can co-exist at the same 
# Visium locations. 

# In practice, it is better to train NMF for a range of factors 𝑅=5,..,30 and select 𝑅
# as a balance between capturing fine-grained and splitting known well-established tissue zones.
# If you want to find a few most disctinct cellular compartments, use a small number of factors. If you want to find very strong 
# co-location signal and assume that most cell types don’t co-locate, use a lot of factors (> 30 - used here).

# Below we show how to perform this analysis. To aid this analysis, we wrapped the analysis shown the notebook on advanced downstream 
# analysis into a pipeline that automates training of the NMF model with varying number of factors:

from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(7, 15), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    # the hyperparameters of NMF can be also adjusted:
    model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)

# Here we plot the NMF weights (Same as saved to `cell_type_fractions_heatmap`)
# res_dict['n_fact12']['mod'].plot_cell_type_loadings()
# plt.savefig(f"{results_folder}/figures/n_fact12_coloc.pdf", bbox_inches='tight')
# plt.close()

# Here we plot the NMF weights (Same as saved to `cell_type_fractions_heatmap`)
# res_dict['n_fact12']['mod'].plot_cell_type_loadings()

# Estimate cell-type specific expression of every gene in the spatial data (needed for NCEM)
# The cell-type specific expression of every gene at every spatial location in the spatial data enables learning cell 
# communication with NCEM model using Visium data (https://github.com/theislab/ncem).
# To derive this, we adapt the approach of estimating conditional expected expression proposed by RCTD (Cable et al) method.
# With cell2location, we can look at the posterior distribution rather than just point estimates of cell type specific expression 
# (see mod.samples.keys() and next section on using full distribution).
# Note that this analysis requires substantial amount of RAM memory and thefore doesn't work on free Google Colab (12 GB limit).

adata_vis.write("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/SHH_SF10210/cell2location/SF10210_vis.h5ad")

mod.samples = adata_vis.uns['mod']
# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
adata_vis.write("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/SHH_SF10210/cell2location/SF10210_vis.h5ad")

# list cell types and genes for plotting
ctypes =['Astro', 'Endo', 'GCP', 'GCP_GN', 'GCP_cycling',
         'GCP_ribo', 'GCP_stress', 'GN_migratory',
         'GN_postmigratory', 'GN_premigratory', 'GN_ribo',
         'Immune', 'Neuron', 'OPC', 'Oligo', 'Pericyte']
# ctypes = ["Endo","GCP_stress","GCP_ribo", "Pericyte"]
genes = ['VEGFA', 'PIK3CA', "CD34","KDR",
         "RICTOR", "HSPB1", "HIF1A", "EIF2AK2","HSPH1", "HSPB1"]

genes = ['VEGFA', 'CD34', 'KDR']

missing_genes = [g for g in genes if g not in adata_vis.var["SYMBOL"].values]
print("Missing genes:", missing_genes)

    # select one slide
    # slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

with mpl.rc_context({'axes.facecolor':  'black'}):
    from tutorial_utils import plot_genes_per_cell_type
    plot_genes_per_cell_type(adata_vis, genes, ctypes)

plt.savefig(f"{results_folder}/figures/interested_genes_each_celltype_2.pdf")
plt.close()

# Example usage
sample_list = [
    "SHH_SF03521", "SHH_SF07994", "SHH_SF08539", "SHH_SF08539.2",
    "SHH_SF09782", "SHH_SF09782.2", "SHH_SF10067", "SHH_SF10210", "SHH_SF12434"
]
base_directory = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python"
snrna_model_directory = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/cell2location/"
snrna_adata_file = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/snRNA_filtered_adata_trained.h5ad"


### Automation
import os
import scanpy as sc
import pandas as pd
from pybiomart import Server

os.chdir("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/")
from c2l_automate import run_cell2location_workflow

# Example usage
sample_names = [
    "SHH_SF03521", "SHH_SF07994", "SHH_SF08539", "SHH_SF08539.2",
    "SHH_SF09782", "SHH_SF09782.2", "SHH_SF10067", "SHH_SF10210", "SHH_SF12434"
]
base_path = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python"
results_root = "/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/"

# adata_ref = sc.read('/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/saveRDS_obj/snRNA_diet.h5ad')
# adata_ref.var['SYMBOL'] = adata_ref.var.index
# adata_ref.X = adata_ref.raw.X

# server = Server(host='http://www.ensembl.org')
# dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

# gene_mapping = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
# gene_mapping_df = pd.DataFrame(gene_mapping)
# mapping_dict = gene_mapping_df.set_index('Gene name')['Gene stable ID'].to_dict()
        
# gene_names = adata_ref.var['SYMBOL'].tolist()
# mapped_genes = [mapping_dict.get(gene, 'Not found') for gene in gene_names]
# adata_ref.var['ENSG'] = mapped_genes
# valid_genes_mask = adata_ref.var['ENSG'] != 'Not found'
# adata_filtered = adata_ref[:, valid_genes_mask].copy()
# adata_filtered.var.set_index('ENSG', drop=True, inplace=True)

# # The default parameters cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12 are a good starting point, however, you can 
# # increase the cut-off to exclude more genes. To preserve marker genes of rare cell types we recommend low cell_count_cutoff=5, however, cell_percentage_cutoff2
# # and nonz_mean_cutoff can be increased to select between 8k-16k genes.
# # In this 2D histogram, orange rectangle highlights genes excluded based on the combination of number of cells expressing that gene (Y-axis) and average RNA count for cells where the gene was detected (X-axis).
# # In this case, the downloaded dataset was already filtered using this method, hence no density under the orange rectangle (to be changed in the future version of the tutorial).
# from cell2location.utils.filtering import filter_genes
# selected = filter_genes(adata_filtered, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# adata_filtered = adata_filtered[:, selected].copy()
# Save anndata object with results
# # Rename _index column in adata.var
if '_index' in adata_filtered.var.columns:
    adata_filtered.var.rename(columns={'_index': 'index_col'}, inplace=True)
             
# If adata.raw exists, rename _index in raw.var as well
if adata_filtered.raw is not None and '_index' in adata_filtered.raw.var.columns:
     adata_filtered.raw.var.rename(columns={'_index': 'index_col'}, inplace=True)

# adata_filtered.write("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/snRNA_save/adata_filtered_removed_genes.h5ad")
adata_filtered_path = os.path.join(base_path, "snRNA_save/adata_filtered_removed_genes.h5ad")
adata_filtered = sc.read(adata_filtered_path)

run_cell2location_workflow(sample_names, base_path, results_root, adata_filtered)


