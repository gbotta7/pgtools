import numpy as np
import gzip
import pandas as pd
import re
import umap

from ./utils import *

# Load data
kc, variant_info, sample_ids = load_kc_matrix(snakemake.input["vcf"])
chrom_vec = np.array([v["chrom"] for v in variant_info])
mapping_prefixes = [id.split('_')[1].replace('pat', '1').replace('hap1', '1').replace('mat', '2').replace('hap2', '2').replace('pri', '0') for id in sample_ids]
chroms_filt = snakemake.params["chroms"]

# Cluster using raw metrics (UMAP + HDBSCAN)
cl_list_raw_umap = {}
for chrom_filt in chroms_filt:
    cl_list_raw_umap[chrom_filt] = {}
    # Filter for a specific chromosome
    chrom_idx = np.where(chrom_vec == chrom_filt)[0]
    chrom_vec_filt = chrom_vec[chrom_idx]
    kc_filt = kc[chrom_idx,:,:]
    for metric in snakemake.params["metrics"]:
        # Compute UMAP embeddings
        cl = compute_umap_hdbscan_clusters(kc = kc_filt,
                                           distance_func = None,
                                           metric = metric,
                                           n_neighbors=snakemake.params["n_neighbors"],
                                           n_components=snakemake.params["n_components"],
                                           random_seed=snakemake.params["random_seed"]
                                           )
        cl_list_raw_umap[chrom_filt][metric] = cl

# Cluster using customized distance function (UMAP + HDBSCAN)
cl_list_custom_umap = {}
for chrom_filt in chroms_filt:
    cl_list_custom_umap[chrom_filt] = {}
    # Filter for a specific chromosome
    chrom_idx = np.where(chrom_vec == chrom_filt)[0]
    chrom_vec_filt = chrom_vec[chrom_idx]
    kc_filt = kc[chrom_idx,:,:]
    for metric in snakemake.params["metrics"]:
        # Compute UMAP embeddings
        cl = compute_umap_hdbscan_clusters(kc = kc_filt,
                                           distance_func = kc_mixed_distance,
                                           inner = metric,
                                           n_neighbors=snakemake.params["n_neighbors"],
                                           n_components=snakemake.params["n_components"],
                                           random_seed=snakemake.params["random_seed"]
                                           )
        cl_list_custom_umap[chrom_filt][metric] = cl

# Cluster using customized distance function (tree + cutree)
cl_list_custom_tree = {}
for chrom_filt in chroms_filt:
    cl_list_custom_tree[chrom_filt] = {}
    # Filter for a specific chromosome
    chrom_idx = np.where(chrom_vec == chrom_filt)[0]
    chrom_vec_filt = chrom_vec[chrom_idx]
    kc_filt = kc[chrom_idx,:,:]
    for metric in snakemake.params["metrics"]:
        for k in range(2, 21):
            cl = compute_tree_clusters(kc,
                                    kc_mixed_distance,
                                    k,
                                    inner = metric,
                                    w_comp=1.0,
                                    w_depth=1.0
                                    )

                                
        
# Create the cluster dataframe
cl_df = pd.DataFrame(
    {f"{chrom}_{kind}_{metric}_umap": cl
     for kind, cl_list in (("raw", cl_list_raw_umap), ("cust", cl_list_custom_umap))
     for chrom, by_metric in cl_list.items()
     for metric, cl in by_metric.items()},
    index=mapping_prefixes
).rename_axis("sample_id").reset_index()

# Add cluster labels to metadata and save to file
META_COLUMNS = ["sample_id", "hap_chrom", "individual", "sex", "continent", "pop_code", "country"]
meta = pd.read_csv(snakemake.input["md_path"], sep="\t", header=None, names=META_COLUMNS)
meta = meta.replace(".", np.nan)
meta['sample_id'] = [id.split('_')[1].replace('pat', '1').replace('hap1', '1').replace('mat', '2').replace('hap2', '2').replace('pri', '0') for id in meta['sample_id']]

# Add counts for each chromosome
kc_tot = kc.sum(axis=2)
log1p_cnt = {f"cnt_{c}": np.log1p(kc_tot[chrom_vec == c].sum(axis=0)) for c in sorted(pd.unique(chrom_vec), key=_chrom_key)}
cnt_df = pd.DataFrame(log1p_cnt, index=mapping_prefixes).rename_axis("sample_id").reset_index()

meta_cnt = meta.merge(cnt_df, on="sample_id", how="left")
meta_cnt = meta_cnt.merge(cl_df, on="sample_id", how="left")