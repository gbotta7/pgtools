import numpy as np
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import umap

def _open_vcf(vcf_path):
    """Transparently open a plain-text or bgzip/gzip-compressed VCF."""
    if vcf_path.endswith(".gz") or vcf_path.endswith(".bgz"):
        return gzip.open(vcf_path, "rt")
    return open(vcf_path, "rt")
    
def load_kc_matrix(vcf_path):
    """
    Parse a VCF with a GT:KC FORMAT field into a 3D numpy array.
    Returns kc[i, j, 0] = first allele count, kc[i, j, 1] = second allele
    count, for variant i, sample j.
    """
    rows = []
    variant_info = []
    sample_ids = None

    with _open_vcf(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_ids = line.lstrip("#").split("\t")[9:]
                continue

            fields = line.split("\t")
            chrom, pos, vid, ref, alt = fields[0:5]
            fmt = fields[8].split(":")
            kc_idx = fmt.index("KC")

            sample_fields = fields[9:]
            counts = np.zeros((len(sample_fields), 2), dtype=np.uint16)
            for j, sf in enumerate(sample_fields):
                kc_str = sf.split(":")[kc_idx]
                a1, a2 = kc_str.split(",")
                counts[j, 0] = int(a1)
                counts[j, 1] = int(a2)

            rows.append(counts)
            variant_info.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt})

    kc = np.stack(rows, axis=0)   # shape: (n_variants, n_samples, 2)
    return kc, variant_info, sample_ids


def _chrom_key(c):
    s = c.replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, s)

def _make_distance(X, distance_fn, inner, w_comp, w_depth):
        """Return a validated square distance matrix."""
        m, n, _ = X.shape
        D = np.asarray(distance_fn(X, inner, w_comp, w_depth), dtype=np.float64)
        assert D.shape == (m, m), f"expected {(m, m)}, got {D.shape}"  
        assert np.all(D >= 0) and np.all(np.isfinite(D)), "bad distances"
        return D


def kc_get_features(X, eps=1.0e-10, dtype=np.float64):
    """
    (n_snps, n_samples, 2) counts -> (n_samples, n_snps, 2) features
    per SNP: [w_swap * alt_fraction, w_depth * log1p(total_depth)]
    alpha: pseudocount shrinking low-depth fractions toward 0.5 (0 = off).
    """
    X = np.asarray(X, dtype=dtype)

    tot = X[:,:,0] + X[:,:, 1]
    comp = X[:,:,0] / (tot + eps)
    log_depth = np.log1p(tot)

    return comp, log_depth


def kc_mixed_distance(X, inner, w_comp=1.0, w_depth=0.01):
    """
    m x m distance over samples, summing per-SNP contributions:
      inner='manhattan': d = sum_k ( w_c*|Δcomp_k| + w_d*|Δdepth_k| )
      inner='hypot'    : d = sum_k sqrt( (w_c*Δcomp_k)^2 + (w_d*Δdepth_k)^2 )
    """
    m, n, _ = X.shape
    comp, log_depth = kc_get_features(X)

    D = np.zeros((m, m), dtype=comp.dtype)
    for k in range(n):
        a = comp[:, k]
        b = log_depth[:, k]
        da = w_comp * np.abs(a[:, None] - a[None, :])
        db = w_depth * np.abs(b[:, None] - b[None, :])
        if inner == "manhattan":
            D += (da + db)
        elif inner == "euclidean":
            D += np.hypot(da, db)

    return D

def compute_umap_hdbscan_clusters(kc,
                                  distance_func = None,
                                  metric = "euclidean",
                                  inner = "euclidean",
                                  w_comp=1.0, 
                                  w_depth=1.0,
                                  n_neighbors=20,
                                  n_components=5,
                                  random_seed=42
                                  ):
    if distance_func:
        D = _make_distance(np.transpose(kc, (1, 0, 2)), distance_func, inner, w_comp, w_depth)
        reducer = umap.UMAP(n_neighbors=min(n_neighbors, D.shape[0] - 1),
                            n_components=n_components,
                            metric=metric,
                            min_dist=0.0,   # helpful for HDBSCAN
                            random_state=random_seed)
        emb = reducer.fit_transform(D)
    else:
        X = np.log1p(np.transpose(kc, (1, 0, 2)))
        Xflat = X.reshape(X.shape[0], -1)
        reducer = umap.UMAP(n_neighbors=min(n_neighbors, Xflat.shape[0] - 1),
                            n_components=n_components,
                            metric=metric,
                            min_dist=0.0,   # helpful for HDBSCAN
                            random_state=random_seed)
        emb = reducer.fit_transform(Xflat)

    # compute clusters with HDBSCAN
    cl = HDBSCAN(min_cluster_size=5, min_samples=5, store_centers="medoid").fit(emb)

    return cl
    

def compute_tree_clusters(kc,
                          distance_fn,
                          k,
                          inner = "euclidean",
                          w_comp=1.0,
                          w_depth=1.0
                          ):


    D = _make_distance(kc, distance_fn, inner, w_comp, w_depth)
    Z = linkage(squareform(D, checks=False), method="ward")

    # compute clusters with tree for max k clusters
    cl = fcluster(Z, t=k, criterion="maxclust")

    return cl