import gzip
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import umap
import sys
import os
from itertools import combinations

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

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
        else:
            print("Distance not known")

    return D


def umap_plot(kc, meta, color_col, legend_title=None,
              distance_func=None, metric="euclidean",
              inner="euclidean", w_comp=1.0, w_depth=0.1,
              n_neighbors=10, min_dist=0.1, n_components=2,
              color_mode="auto", cmap_numeric="viridis", robust=True,
              figsize_per_panel=(5, 5), save_path=None):
    """UMAP of samples from a single-chromosome count array.

    kc : array (n_variants, n_samples, n_features)
    meta : DataFrame with one row per sample, aligned to kc.shape[1]
    color_mode : 'auto' | 'categorical' | 'numeric'

    Draws every pairwise projection of the n_components-dim embedding
    (a single panel when n_components == 2).
    """
    legend_title = legend_title or color_col

    # ---- embedding ---------------------------------------------------
    X = np.transpose(kc, (1, 0, 2))          # samples first

    if distance_func is not None:
        D = np.asarray(distance_func(X, inner, w_comp, w_depth), dtype=np.float64)
        m = X.shape[0]
        assert D.shape == (m, m), f"expected {(m, m)}, got {D.shape}"
        assert np.all(D >= 0) and np.all(np.isfinite(D)), "bad distances"
        data = D
    else:
        data = np.log1p(X).reshape(X.shape[0], -1)

    reducer = umap.UMAP(n_neighbors=min(n_neighbors, data.shape[0] - 1),
                        min_dist=min_dist, n_components=n_components,
                        metric=metric, random_state=42)
    emb = reducer.fit_transform(data)

    # ---- colors ------------------------------------------------------
    vals = meta[color_col]
    if color_mode == "auto":
        color_mode = "numeric" if pd.api.types.is_numeric_dtype(vals) else "categorical"

    if color_mode == "numeric":
        v = pd.to_numeric(vals, errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(v)
        if not finite.any():
            raise ValueError(f"{color_col} has no finite values")
        lo, hi = (np.nanpercentile(v[finite], [2, 98]) if robust
                  else (v[finite].min(), v[finite].max()))
        if lo == hi:
            hi = lo + 1e-9
        norm = mpl.colors.Normalize(vmin=lo, vmax=hi)
        cmap = plt.get_cmap(cmap_numeric)
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    else:
        labels = vals.fillna("Unknown")
        categories = sorted(labels.unique(), key=str)
        cat_cmap = plt.get_cmap("tab20", max(len(categories), 1))
        color_map = {c: cat_cmap(k) for k, c in enumerate(categories)}

    def _draw(ax, x, y):
        if color_mode == "numeric":
            if (~finite).any():
                ax.scatter(x[~finite], y[~finite], color="lightgrey",
                           s=15, alpha=0.6, zorder=1)
            ax.scatter(x[finite], y[finite], c=v[finite], cmap=cmap, norm=norm,
                       s=15, alpha=0.85, linewidths=0, zorder=2)
        else:
            for cat in categories:
                m = (labels == cat).to_numpy()
                ax.scatter(x[m], y[m], color=color_map[cat], s=15, alpha=0.8)

    # ---- figure ------------------------------------------------------
    pairs = list(combinations(range(n_components), 2))
    fig, axes = plt.subplots(1, len(pairs),
                             figsize=(figsize_per_panel[0] * len(pairs),
                                      figsize_per_panel[1]),
                             squeeze=False)
    axes = axes[0]

    for ax, (a, b) in zip(axes, pairs):
        _draw(ax, emb[:, a], emb[:, b])
        ax.set_xlabel(f"UMAP{a + 1}")
        ax.set_ylabel(f"UMAP{b + 1}")
        if len(pairs) > 1:
            ax.set_title(f"{a + 1} vs {b + 1}", fontsize=9)
        ax.set_aspect("equal", adjustable="datalim")

    if color_mode == "numeric":
        cbar = fig.colorbar(sm, ax=list(axes), fraction=0.03, pad=0.02)
        cbar.set_label(legend_title)
        if robust:
            cbar.ax.set_title("2–98%", fontsize=7)
    else:
        handles = [plt.Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=color_map[c], markersize=8, label=c)
                   for c in categories]
        fig.legend(handles=handles, title=legend_title,
                   bbox_to_anchor=(1.02, 0.5), loc="center left", fontsize=8)
        fig.tight_layout()

    fig.suptitle(f"UMAP ({n_components}D, n={kc.shape[0]} variants), "
                 f"colored by {legend_title}", y=1.02)

    if save_path is None:
        save_path = f"umap_{color_col}_{n_components}d.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

CHROM = snakemake.params["chrom"]

kc, variant_info, sample_ids = load_kc_matrix(snakemake.input["vcf"])

for v in variant_info:
    v["chrom"] = CHROM

META_COLUMNS = [
    "sample_id",   # matches VCF sample_ids, e.g. "120001_CN1.pat"
    "hap_chrom",   # sex chromosome carried by this haplotype: X / Y / XY / "."
    "individual",  # individual name shared between the two haplotypes (pat/mat or hap1/hap2)
    "sex",         # M / F / "."
    "superpop",    # continental group: EastAsia, WestEurasia, SouthAsia, America, Africa, "."
    "pop_code",    # 1000 Genomes-style population code (e.g. CHS, PUR, LWK), "." if not in that panel
    "country",     # country of origin, "." if not recorded
]

meta_path = snakemake.input["md_path"]

meta = pd.read_csv(meta_path, sep="\t", header=None, names=META_COLUMNS)
meta = meta.replace(".", np.nan)
meta['sample_id'] = [id.split('_')[1].replace('pat', '1').replace('hap1', '1').replace('mat', '2').replace('hap2', '2').replace('pri', '0') for id in meta['sample_id']]

mapping_prefixes = [id.split('_')[1].replace('pat', '1').replace('hap1', '1').replace('mat', '2').replace('hap2', '2').replace('pri', '0') for id in sample_ids]
chrom_cnt = np.log1p(kc.sum(axis=2).sum(axis=0))
cnt_df = (pd.DataFrame({"log1p_chrom_cnt": chrom_cnt}, index=mapping_prefixes)
            .rename_axis("sample_id")
            .reset_index())

meta = meta.merge(cnt_df, on="sample_id", how="inner")

os.makedirs(snakemake.output["plotdir"], exist_ok=True)
for attr in ["hap_chrom", "sex", "superpop", "log1p_chrom_cnt"]:
    for metric in ["manhattan", "euclidean"]:
        save_path = os.path.join(snakemake.output["plotdir"], f"umap_{attr}_{metric}.png")
        color_mode = "auto" if attr == "log1p_chrom_cnt" else "categorical"
        umap_plot(kc,
                meta,
                color_col=attr,
                legend_title=attr,
                color_mode=color_mode,
                distance_func=kc_mixed_distance,
                w_comp=1.0,
                w_depth=1.0,
                metric="precomputed",
                inner=metric,
                n_components=2,
                n_neighbors=15,
                save_path=save_path)