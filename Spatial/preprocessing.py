def process_and_visualize_data(data_path, savedir, filename):
  # Import libraries
  import scanpy as sc
  import pandas as pd
  import matplotlib.pyplot as plt
  import seaborn as sns
  import os
  
  # Load dataset
  # adata = sc.read_visium(data_path)
  adata = sc.read(data_path)
  adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
  adata.var_names_make_unique()
  adata.var["mt"] = adata.var_names.str.startswith("MT-")
  sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
  
  # QC and preprocessing
  fig, axs = plt.subplots(2, 3, figsize=(10, 8))
  sns.histplot(adata.obs["total_counts"], kde=True, ax=axs[0, 0])
  sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 200], kde=False, bins=40, ax=axs[0, 1])
  sns.histplot(adata.obs["pct_counts_mt"], kde=True, ax=axs[0, 2])
  sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1, 0])
  sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[1, 1])
  sns.histplot(adata.var["n_cells_by_counts"], kde=True, ax=axs[1, 2])
  
  os.chdir(savedir)
  os.makedirs(filename, exist_ok=True)
  os.chdir("./"+filename)
  plt.savefig("qc_plots.png")
  
  # sc.pp.filter_cells(adata, min_counts=500) To maintain the same thing with the Seurat filtering
  sc.pp.filter_cells(adata, min_genes=500)
  adata = adata[adata.obs["pct_counts_mt"] < 20]
  print(f"Cells after min genes > 500 and MT filter < 20 filter: {adata.n_obs}")
  sc.pp.filter_genes(adata, min_cells=10, inplace=True)
  print(f"Genes after min cells filter = 10: {adata.n_vars}")
  
  # Normalization and highly variable genes
  sc.pp.normalize_total(adata, inplace=True)
  sc.pp.log1p(adata)
  sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
  
  # Manifold embedding and clustering
  sc.pp.pca(adata)
  sc.pp.neighbors(adata)
  sc.tl.umap(adata)
  sc.tl.leiden(adata, key_added="clusters")
  
  # Plotting UMAP
  plt.rcParams["figure.figsize"] = (4, 4)
  sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
  plt.savefig("UMAP.png")
  # plt.show()
  
  # Visualization in spatial coordinates
  plt.rcParams["figure.figsize"] = (8, 8)
  sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts", "clusters"], bw=False, alpha_img=0, size=1.3)
  plt.savefig("spatial.png")
  # plt.show()
  
  # Cluster-specific visualization
  plt.rcParams["figure.figsize"] = (10, 8)
  sc.pl.spatial(adata, img_key="hires", color="clusters", groups=["1", "3"], alpha_img=0, size=1.3)
  plt.savefig("clus_1_3.png")
  # plt.show()
  
  # Cluster marker genes
  sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
  sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby="clusters", show_gene_labels=True)
  plt.savefig("heatmap_markers.png")
  # plt.show()
  
  adata.write("adata_save.h5ad")
  
  return adata
  
  
