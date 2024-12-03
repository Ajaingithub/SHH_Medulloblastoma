def run_cell2location_workflow(sample_names, base_path, results_root, adata_filtered):
    """
    Automates the Cell2location workflow for a list of samples.

    Args:
        sample_names (list): List of sample names to process.
        base_path (str): Base path containing the input files for each sample.
        results_root (str): Root folder to save the results for each sample.
        snrna_model_path (str): Path to the trained snRNA model.

    Returns:
        None
    """
    for sample in sample_names:
        # loading packages
        import os
        import scanpy as sc
        import cell2location as c2l
        import matplotlib.pyplot as plt
        import numpy as np
        import torch

        print(f"Processing sample: {sample}")
        results_folder = os.path.join(results_root, sample, "cell2location")
        os.makedirs(results_folder, exist_ok=True)

        # Load spatial data
        adata_vis_path = os.path.join(base_path, sample, "saveh5ad", f"{sample}_save_2.h5ad")
        adata_vis = sc.read(adata_vis_path)
        adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
        adata_vis.var['SYMBOL'] = adata_vis.var_names
        adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

        # Save a spatial plot for VEGFA
        os.chdir(results_folder)
        sc.pl.spatial(adata_vis, color='VEGFA', gene_symbols='SYMBOL', save=f'VEGFA_{sample}.pdf')

        # Remove mitochondrial genes
        adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
        adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
        adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

        # Load filtered snRNA data
        shared_features = [feature for feature in adata_vis.var_names if feature in adata_filtered.var_names]
        adata_filtered = adata_filtered[:, shared_features].copy()
        adata_vis = adata_vis[:, shared_features].copy()

        ### Model Training
        c2l.models.RegressionModel.setup_anndata(
            adata=adata_filtered,
            batch_key="SampleID",
            labels_key="cell_types",
            categorical_covariate_keys=["batch"],
        )

        # create the regression model
        model = c2l.models.RegressionModel(adata_filtered)

        # view anndata_setup as a sanity check
        print(model.view_anndata_setup())

        # Check if GPU is available
        use_gpu = torch.cuda.is_available()
        print("GPU available", use_gpu)
        
        # Set precision for float32 matrix multiplication to medium or high
        torch.set_float32_matmul_precision('medium')  # or 'high'

        # Set the number of workers in the datasplitter_kwargs
        datasplitter_kwargs = {
            'num_workers': 50  
            }
        
        trainer_kwargs = {
            'max_epochs': 250,
            'batch_size': 2500,
            'train_size': 1,
            'lr': 0.002,
            'datasplitter_kwargs': datasplitter_kwargs
        }
        
        model.train(**trainer_kwargs)

        model.plot_history(20)
        plt.savefig(f'{results_folder}/figures/{sample}_training_elbow_plot.pdf')
        plt.close()
        
        model.save(f'{results_folder}/{sample}_snRNA_ref_model')

        # Export posterior from snRNA model
        adata_filtered = model.export_posterior(
            adata_filtered, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
        )

        # Save anndata object with results
        # # Rename _index column in adata.var
        if '_index' in adata_filtered.var.columns:
             adata_filtered.var.rename(columns={'_index': 'index_col'}, inplace=True)
             
        # If adata.raw exists, rename _index in raw.var as well
        if adata_filtered.raw is not None and '_index' in adata_filtered.raw.var.columns:
             adata_filtered.raw.var.rename(columns={'_index': 'index_col'}, inplace=True)
             
        # Extracting reference cell types signatures as a pd.DataFrame.
        # All parameters of the a Negative Binomial regression model are exported into reference anndata object, however for spatial mapping we just
        # need the estimated expression of every gene in every cell type. Here we extract that from standard output:
        
        # export estimated expression in each cluster
        if 'means_per_cluster_mu_fg' in adata_filtered.varm.keys():
             inf_aver = adata_filtered.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                                         for i in adata_filtered.uns['mod']['factor_names']]].copy()
        else:
             inf_aver = adata_filtered.var[[f'means_per_cluster_mu_fg_{i}'
                                             for i in adata_filtered.uns['mod']['factor_names']]].copy()
             
        inf_aver.columns = adata_filtered.uns['mod']['factor_names']

        # Prepare spatial data for Cell2location
        intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
        adata_vis = adata_vis[:, intersect].copy()
        inf_aver = inf_aver.loc[intersect, :].copy()
        # cell_signatures = adata_filtered.varm['means_per_cluster_mu_fg'][intersect, :].copy()
        adata_vis.X = adata_vis.layers["raw_counts"]

        c2l.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

        # Train the Cell2location model
        model = c2l.models.Cell2location(
            adata_vis, cell_state_df=inf_aver,
            N_cells_per_location=30, detection_alpha=20
        )
        print(model.view_anndata_setup())

        datasplitter_kwargs = {
            'num_workers': 50  
            }
        
        trainer_kwargs = {
            'max_epochs': 30000,
            'batch_size': None,
            'train_size': 1,
            'datasplitter_kwargs': datasplitter_kwargs
        }

        model.train(**trainer_kwargs)

        # Save model and QC plot
        model.save(os.path.join(results_folder,f"{sample}_spatial_model"))
        model.plot_history(1000)
        plt.savefig(os.path.join(results_folder, "figures", f"{sample}_training_elbow_plot.pdf"))

        # Export posterior for spatial data
        adata_vis = model.export_posterior(
            adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': adata_vis.n_obs}
        )
        adata_vis.write(os.path.join(results_folder, f"{sample}_vis.h5ad"))

        # Plot QC metrics
        plt.close()
        model.plot_QC()
        plt.savefig(os.path.join(results_folder, "figures", f"{sample}_spatial_plot_QC.pdf"))

        print(f"Finished processing sample: {sample}")
        
    # -------------------------------------------
    # Confidence in Cell Abundance: Adding 5% Quantile
    # -------------------------------------------
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

    # Spatial Plots: Visualizing Cell Abundance
    with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(
            adata_vis,
            cmap='viridis',
            color=[
                'Astro', 'Endo', 'GCP', 'GCP_GN', 'GCP_cycling', 'GCP_ribo',
                'GCP_stress', 'GN_migratory', 'GN_postmigratory', 'GN_premigratory',
                'GN_ribo', 'Immune', 'Neuron', 'OPC', 'Oligo', 'Pericyte'
            ],
            ncols=4, size=1.3,
            img_key='hires',
            vmin=0, vmax='p99.2'
        )
    plt.savefig(f'{results_folder}/figures/spatial_each_cell_featureplots_viridis.pdf')
    plt.close()

    # Spatial Visualization for Selected Cell Clusters
    clust_labels = ['Endo', 'Pericyte', 'GCP_stress']
    clust_col = [str(i) for i in clust_labels]

    with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': (15, 15)}):
        fig = plot_spatial(
            adata=adata_vis,
            color=clust_col,
            labels=clust_labels,
            show_img=True,
            style='fast',
            max_color_quantile=0.992,
            circle_diameter=7,
            colorbar_position='right'
        )
    plt.savefig(f'{results_folder}/figures/{sample}_non_malign_GCP_stress.pdf')
    plt.close()

    # -------------------------------------------
    # Identifying Tissue Regions Using Clustering
    # -------------------------------------------
    # Compute KNN and Leiden Clustering
    sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf', n_neighbors=15)
    sc.tl.leiden(adata_vis, resolution=1.1)
    adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")

    # UMAP Representation
    sc.tl.umap(adata_vis, min_dist=0.3, spread=1)

    # UMAP by Region
    with mpl.rc_context({'axes.facecolor': 'white', 'figure.figsize': [8, 8]}):
        sc.pl.umap(
            adata_vis,
            color=['region_cluster'],
            size=15,
            color_map='RdPu',
            ncols=2,
            legend_loc='on data',
            legend_fontsize=10
        )
    plt.savefig(f"{results_folder}/figures/{sample}_spatial_UMAP_cluster_cellabundance.pdf")
    plt.close()

    # UMAP by Sample
    with mpl.rc_context({'axes.facecolor': 'white', 'figure.figsize': [8, 8]}):
        sc.pl.umap(
            adata_vis,
            color=['sample'],
            size=30,
            color_map='RdPu',
            ncols=2,
            legend_fontsize=20
        )
    plt.savefig(f"{results_folder}/figures/{sample}_spatial_UMAP_cluster_sample.pdf")
    plt.close()

    # High-Resolution Spatial Clustering
    with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(
            adata_vis,
            color=['region_cluster'],
            size=1.3,
            img_key='hires',
            alpha=0.5
        )
    plt.savefig(f"{results_folder}/figures/{sample}_spatial_cluster_high_res_imag.pdf", bbox_inches='tight')
    plt.close()

    # -------------------------------------------
    # Identifying Tissue Zones Using NMF
    # -------------------------------------------
    res_dict, adata_vis = run_colocation(
        adata_vis,
        model_name='CoLocatedGroupsSklearnNMF',
        train_args={
            'n_fact': np.arange(7, 15),
            'sample_name_col': 'sample',
            'n_restarts': 3
        },
        model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
        export_args={'path': f'{run_name}/CoLocatedComb/'}
    )

    # -------------------------------------------
    # Gene and Cell Type-Specific Expression
    # -------------------------------------------
    adata_vis.write(f"{results_folder}{sample}_vis.h5ad")

    mod.samples = adata_vis.uns['mod']
    expected_dict = mod.module.model.compute_expected_per_cell_type(
        mod.samples["post_sample_q05"], mod.adata_manager
    )

    # Add Cell-Type Specific Expression Layers
    for i, n in enumerate(mod.factor_names_):
        adata_vis.layers[n] = expected_dict['mu'][i]

    adata_vis.write("/diazlab/data3/.abhinav/projects/SHH/Visium/3rd_strict_cutoff/python/SHH_SF10210/cell2location/SF10210_vis.h5ad")

    # Plotting Interested Genes for Each Cell Type
    ctypes = ['Astro', 'Endo', 'GCP', 'GCP_GN', 'GCP_cycling', 'GCP_ribo', 
            'GCP_stress', 'GN_migratory', 'GN_postmigratory', 'GN_premigratory', 
            'GN_ribo', 'Immune', 'Neuron', 'OPC', 'Oligo', 'Pericyte']
    genes = ['VEGFA', 'CD34', 'KDR']
    missing_genes = [g for g in genes if g not in adata_vis.var["SYMBOL"].values]
    print("Missing genes:", missing_genes)

    with mpl.rc_context({'axes.facecolor': 'black'}):
        from tutorial_utils import plot_genes_per_cell_type
        plot_genes_per_cell_type(adata_vis, genes, ctypes)

    plt.savefig(f"{results_folder}/figures/interested_genes_each_celltype_2.pdf")
    plt.close()
