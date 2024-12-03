import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

os.chdir("/Users/M256617/Documents/Industry/UCSF/DIaz_lab/projects/SHH/Visium")
### Checking the directory location

# The function datasets.visium_sge() downloads the dataset from 10x Genomics and returns an AnnData object that contains counts, images and spatial 
# coordinates. We will calculate standards QC metrics with pp.calculate_qc_metrics and percentage of mitochondrial read counts per sample.
adata = sc.read("/Users/M256617/Documents/Industry/UCSF/DIaz_lab/projects/SHH/Visium/scanpy/SHH_SF03521.h5ad")
adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

### QC and preprocessing
fig, axs = plt.subplots(2, 3, figsize=(10, 8))
sns.histplot(adata.obs["total_counts"], kde = True, ax = axs[0,0])
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 2000], kde = True, ax = axs[0,1])
sns.histplot(adata.obs["pct_counts_mt"], kde = False, ax = axs[0,2])
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax = axs[1,0])
sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 2500], kde=False, ax = axs[1,1])
sns.histplot(adata.var["n_cells_by_counts"], kde = True, ax = axs[1,2])

savedir = "/Users/M256617/Documents/Industry/UCSF/DIaz_lab/projects/SHH/Visium/"
os.makedirs((savedir+"SF03521"), exist_ok=True)
plt.savefig(savedir + "SF03521/qc_plots.png")  # Replace "path/to/save/plots.png" with your desired file path and format
sc.pl.spatial(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],save='spatial_plot.png')

import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Directory containing the .h5ad files
data_dir = "/Users/M256617/Documents/Industry/UCSF/DIaz_lab/projects/SHH/Visium/scanpy"
output_dir = "/Users/M256617/Documents/Industry/UCSF/DIaz_lab/projects/SHH/Visium/"

# Iterate over all .h5ad files in the directory
for filename in os.listdir(data_dir):
    if filename.endswith(".h5ad"):
        file_path = os.path.join(data_dir, filename)
        adata = sc.read(file_path)
        # Process the AnnData object
        adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        # QC and preprocessing plots
        fig, axs = plt.subplots(3, 3, figsize=(12, 10))
        sns.histplot(adata.obs["total_counts"], kde=True, ax=axs[0, 0])
        sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 2000], bins = 50, kde=True, ax=axs[0, 1])
        sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 1000], bins = 50, kde=True, ax=axs[0, 2])
        # sns.histplot(adata.obs["pct_counts_mt"], kde=False, ax=axs[0, 2])
        sns.histplot(adata.obs["n_genes_by_counts"], kde=True, bins=60, ax=axs[1, 0])
        sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 2000], bins = 50, kde=True, ax=axs[1, 1])
        sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 1000], bins = 50, kde=True, ax=axs[1, 2])
        sns.histplot(adata.var["n_cells_by_counts"], kde=True, ax=axs[2, 0])
        sns.histplot(adata.var["n_cells_by_counts"][adata.var["n_cells_by_counts"] < 800], bins = 50, kde=True, ax=axs[2, 1])
        sns.histplot(adata.var["n_cells_by_counts"][adata.var["n_cells_by_counts"] < 500], bins = 50, kde=True, ax=axs[2, 2])
        # Create output directory for each file
        file_output_dir = os.path.join(output_dir, filename.split('.')[0])
        os.makedirs(file_output_dir, exist_ok=True)
        # Save QC plots
        plt.savefig(os.path.join(file_output_dir, (filename.split('.')[0]+"_qc_plots.png")))
        plt.close(fig)  # Close the figure to free memory
        # Save spatial plots
        os.chdir(file_output_dir)
        sc.pl.spatial(adata, color=['n_genes_by_counts', 'total_counts'], save=(filename.split('.')[0]+'spatial_plot.png'), show=False)
        sc.pl.spatial(adata, color=['n_genes_by_counts', 'total_counts'], vmin =0, vmax=2000, save=(filename.split('.')[0]+'spatial_plot_1.png'), show=False)
        sc.pl.spatial(adata, color=['n_genes_by_counts', 'total_counts'], vmin =0, vmax=1000, save=(filename.split('.')[0]+'spatial_plot_2.png'), show=False)
        sc.pl.spatial(adata, color=['n_genes_by_counts', 'total_counts'], vmin =0, vmax=500, save=(filename.split('.')[0]+'spatial_plot_3.png'), show=False)
        plt.close(fig)


adata.layers["raw_counts"].data