def process_and_visualize_data(data_path, savedir, filename):
  # Import libraries
  import scanpy as sc
  import pandas as pd
  import matplotlib.pyplot as plt
  import seaborn as sns
  import os
  
  # Load dataset
  # Determine if the data path is a file (h5ad) or a directory (Cell Ranger output)
  if data_path.endswith(".h5ad"):
    adata = sc.read(data_path)
  else:
    adata = sc.read_visium(data_path)
  
  adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
  adata.var_names_make_unique()
  adata.var["mt"] = adata.var_names.str.startswith("MT-")
  sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
  
  # QC and preprocessing
  fig, axs = plt.subplots(2, 3, figsize=(10, 8))
  sns.histplot(adata.obs["total_counts"], kde=True, ax=axs[0, 0])
  sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 1000], kde=False, bins=40, ax=axs[0, 1])
  sns.histplot(adata.obs["pct_counts_mt"], kde=True, ax=axs[0, 2])
  sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1, 0])
  sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[1, 1])
  sns.histplot(adata.var["n_cells_by_counts"], kde=True, ax=axs[1, 2])
  
  os.chdir(savedir)
  os.makedirs(filename, exist_ok=True)
  os.chdir("./"+filename)
  os.makedirs("./figures", exist_ok=True)
  plt.savefig("./figures/"+filename+"qc_plots.pdf")
  plt.close()

  # Making a Violin plots
  fig_vln, axs_vln = plt.subplots(2, 3, figsize=(15, 10))
  sns.violinplot(y=adata.obs["total_counts"], ax=axs_vln[0, 0])
  sns.violinplot(y=adata.obs["total_counts"][adata.obs["total_counts"] < 1000], ax=axs_vln[0, 1])
  sns.violinplot(y=adata.obs["pct_counts_mt"], ax=axs_vln[0, 2])
  sns.violinplot(y=adata.obs["n_genes_by_counts"], ax=axs_vln[1, 0])
  sns.violinplot(y=adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], ax=axs_vln[1, 1])
  sns.violinplot(y=adata.var["n_cells_by_counts"], ax=axs_vln[1, 2])
    
  plt.savefig("./figures/"+filename+"qc_plots_vln.pdf")
  plt.close()
  
  ## Before fitlering the cell counts
  cells_before = adata.n_obs
  genes_before = adata.n_vars
  
  # sc.pp.filter_cells(adata, min_counts=500) To maintain the same thing with the Seurat filtering
  sc.pp.filter_cells(adata, min_genes=500)
  adata = adata[adata.obs["pct_counts_mt"] < 20]
  print(f"Cells after min genes > 500 and MT filter < 20 filter: {adata.n_obs}")
  sc.pp.filter_genes(adata, min_cells=10, inplace=True)
  print(f"Genes after min cells filter = 10: {adata.n_vars}")
  cells_after = adata.n_obs
  genes_after = adata.n_vars

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
  sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save="_"+filename+".pdf", show = False)
  # plt.savefig("UMAP.png")
  # plt.show()
  
  # Visualization in spatial coordinates
  plt.rcParams["figure.figsize"] = (8, 8)
  sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts", "clusters"], 
                # bw=False, 
                # alpha_img=0, 
                size=1.3, 
                save = "_spatial_"+filename+".pdf", 
                show = False)
  # plt.savefig("spatial.png")
  # plt.show()
  
  # Cluster-specific visualization
  # plt.rcParams["figure.figsize"] = (10, 8)
  # sc.pl.spatial(adata, img_key="hires", color="clusters", groups=["1", "3"], alpha_img=0, size=1.3)
  # plt.savefig("clus_1_3.png")
  # plt.show()
  
  # Cluster marker genes
  sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
  sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby="clusters", show_gene_labels=True, save = "_"+filename+".pdf", show = False)
  # plt.savefig("heatmap_markers.png")
  # plt.show()
  
  os.makedirs("./saveh5ad", exist_ok=True)
  adata.write("saveh5ad/"+filename+"_save.h5ad")
  
  return adata, cells_before, cells_after, genes_before, genes_after
  
  
