# Function to process a single sample
def process_sample(sample, matched_dir, savedir, image_base_dir,celltype_file_base_dir,cluster_features,CT_colors):
    import os
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import squidpy as sq
    from collections import Counter
    
    print(f"Processing sample: {sample}")
    
    # Load data
    adata = sc.read(f"{savedir}{sample}/saveh5ad/{sample}_save.h5ad")
    img_path = f"{image_base_dir}{matched_dir}/{matched_dir}.jpg"

    
    """
    calculate_image_features function extracts features from spatial transcriptomics data's associated histology images. This process involves analyzing
    pixel intensities, textures, and other spatial properties of the image within certain areas of interest (e.g., cell neighborhoods).
     
    Feature Extraction Workflow
    The function computes features for each spatial region corresponding to the barcodes in adata. These features could include:
        Intensity features: Mean, median, standard deviation of pixel intensities.
        Texture features: Haralick features (e.g., entropy, contrast) calculated using the Gray-Level Co-occurrence Matrix (GLCM).
        Edge features: Detection of boundaries or gradients in the image.
        Other statistics: Counts, ratios, or distributions of pixel values.
    """

    # Image processing
    if os.path.exists(img_path):
        img = sq.im.ImageContainer(img_path) # The ImageContainer is a versatile object that Squidpy uses to store and manage image data. It can handle raw images, preprocessed layers, and associated metadata.
        feature_dfs = []  # List to store the individual feature dataframes
        for scale in [1.0, 2.0]:
            feature_name = f"features_summary_scale{scale}"
            sq.im.calculate_image_features(
                adata,
                img.compute(),
                features="summary", # Requests summary-level features (default) like intensity, texture, and edge information.
                key_added=feature_name,
                n_jobs=4,
                scale=scale,
            )
            
            feature_df = adata.obsm[feature_name]
            feature_df.columns = [f"{col}_scale{scale}" for col in feature_df.columns]
            feature_dfs.append(feature_df)

        # Combine features into one dataframe
        # adata.obsm["features"] = pd.concat(
        #     [adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f],
        #     axis="columns",
        # )
        # adata.obsm["features"].columns = pd.Index(
        #     adata.obsm["features"].columns.str.replace(" ", "_", regex=True)
        # )
    
    # Combine all features into one dataframe with unique column names
    adata.obsm["features"] = pd.concat(feature_dfs, axis="columns")

    # Clean column names by replacing spaces with underscores (if any)
    adata.obsm["features"].columns = pd.Index(
        adata.obsm["features"].columns.str.replace(" ", "_", regex=True)
    )
    
    # Feature clustering
    os.makedirs(f"{savedir}{sample}/figures", exist_ok=True)
    adata.obs["features_cluster_0.2"] = cluster_features(adata.obsm["features"], like="summary", res = 0.2)
    adata.obs["features_cluster_0.5"] = cluster_features(adata.obsm["features"], like="summary", res = 0.5)
    adata.obs["features_cluster_0.8"] = cluster_features(adata.obsm["features"], like="summary", res = 0.8)
    adata.obs["features_cluster_1"] = cluster_features(adata.obsm["features"], like="summary", res = 1)
    adata.obs["features_cluster_1.2"] = cluster_features(adata.obsm["features"], like="summary", res = 1.2)

    cluster_columns = ["features_cluster_0.2","features_cluster_0.5","features_cluster_0.8","features_cluster_1","features_cluster_1.2"]
    fig, axes = plt.subplots(1, len(cluster_columns), figsize=(20, 5)) ## Creating a figure with multiple subplots
    
    # Plot each resolution in a separate subplot
    for ax, cluster_column in zip(axes, cluster_columns):
        sq.pl.spatial_scatter(adata, color=[cluster_column], size=1.3, ax=ax)
        ax.set_title(f"Resolution {cluster_column.split('_')[-1]}")
    
    # sq.pl.spatial_scatter(adata, color=["features_cluster"], size=1.3)
    plt.tight_layout()
    plt.savefig(f"{savedir}{sample}/figures/{sample}_feature_extraction.pdf", format = "pdf")
    
    # Celltype merging
    celltype_file = f"{celltype_file_base_dir.format(sample.replace('SHH_', '').replace('.2', ''))}celltype_cellid.txt"
    if os.path.exists(celltype_file):
        celltype_clus = pd.read_csv(celltype_file, delimiter="\t")
        celltype_clus.index = celltype_clus["cell_barcode"]
        adata.obs = pd.merge(adata.obs, celltype_clus, left_index=True, right_index=True, how="left")
        adata.obs["Predicted_id"] = adata.obs["Predicted_id"].fillna("GCP")
        cluster_key = "Predicted_id"
        # Ensure that the cluster_key column is of the 'category' dtype
        adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')
        adata.uns[f"{cluster_key}_colors"] = [CT_colors[cluster] for cluster in adata.obs[cluster_key].cat.categories]
        sc.pl.spatial(adata, img_key="hires", color=["Predicted_id"], size=1.5)
        plt.tight_layout()
        plt.savefig(f"{savedir}{sample}/figures/{sample}_spatial_celltypes.pdf", format="pdf")
    
        # Neighborhood enrichment
        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key="Predicted_id")
        sq.pl.nhood_enrichment(adata, cluster_key="Predicted_id")
        plt.savefig(f"{savedir}{sample}/figures/{sample}_neighborhood_enrichment_predicted.pdf", bbox_inches="tight", format="pdf")
        
        # Co-occurence
        sq.gr.co_occurrence(adata, cluster_key="Predicted_id")
        sq.pl.co_occurrence(adata, cluster_key="Predicted_id")
        plt.savefig(f"{savedir}{sample}/figures/{sample}_neighborhood_cooccurence_predicted.pdf", bbox_inches="tight", format="pdf")

    adata.write(f"{savedir}{sample}/saveh5ad/{sample}_save_2.h5ad")


