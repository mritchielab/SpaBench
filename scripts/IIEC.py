#!/usr/bin/env python3
"""
IIEC.py

Compute alpha-vs-overall_prop curves for many clustering results and plot grouped by clustering method.

Usage:
  python IIEC.py --adata /path/to/adata.h5ad \
                 --clustering_result /path/to/results/Mouse_brain_squidpy \
                 --out_dir /path/to/output

Outputs:
  - PNG files: one per clustering method (e.g. stlearn) in out_dir/figs/
  - CSV file: out_dir/alpha_curves.csv containing curve data for each workflow

Notes:
  - The script expects clustering CSVs with barcodes as the first column (index) and a single cluster column.
  - Ground truth labels are expected in adata.obs['cluster']. If missing, specify via --ground_truth.
"""

import argparse
from pathlib import Path
import logging
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity, pairwise_distances
import matplotlib.pyplot as plt
from tqdm import tqdm
import scanpy as sc

# Default selected features (from your example). If not present in adata.obsm['img_feats'], we'll use intersection.
DEFAULT_SELECTED_FEATS = [
    "texture_ch-2_contrast_dist-1_angle-1.57",
    "summary_ch-2_quantile-0.5",
    "summary_ch-0_quantile-0.5",
    "summary_ch-1_quantile-0.5",
    "texture_ch-2_correlation_dist-1_angle-1.57",
    "texture_ch-2_dissimilarity_dist-1_angle-1.57",
    "summary_ch-2_quantile-0.9",
]


def parse_args():
    p = argparse.ArgumentParser(description="Compute alpha curves for many clustering results and plot grouped by clustering method")
    p.add_argument("--adata", required=True, help=".h5ad file with image features in obsm['img_feats'] and spatial coords in obsm['spatial']")
    p.add_argument("--clustering_result", required=True, help="Base folder containing clustering CSV files (searched recursively)")
    p.add_argument("--out_dir", required=True, help="Output directory for figures and CSV")
    p.add_argument("--ground_truth", default="cluster", help="Column name in adata.obs for ground truth labels (default: 'cluster')")
    p.add_argument("--bins", type=int, default=101, help="Number of alpha points (default 101 for 0..100%%)")
    p.add_argument("--selected_feats", nargs="*", default=None, help="List of selected feature names to use (if not provided, defaults are used)")
    p.add_argument("--min_overlap", type=float, default=0.5, help="Minimum fraction of barcodes that must overlap between cluster file and adata (default 0.5)")
    return p.parse_args()


def find_cluster_files(base_dir: Path):
    return list(base_dir.rglob("*clusters.csv"))


def read_cluster_csv(fp: Path):
    try:
        df = pd.read_csv(fp, index_col=0, header=0)
        # If multiple columns, pick the last one
        if df.shape[1] >= 1:
            col = df.columns[-1]
            s = df[col].astype(str)
            return s
    except Exception as e:
        logging.warning(f"Failed to read {fp}: {e}")
    return None


def alpha_curve_for_codes(codes, sims, dists, d_max, alphas, eps=1e-8):
    # codes: integer codes array length n
    n = sims.shape[0]
    props = []
    for alpha in alphas:
        D = alpha * (1 - sims) + (1 - alpha) * (dists / d_max)
        s = np.exp(-D)
        N = np.full(n, np.nan)
        for i in range(n):
            same = np.where((codes == codes[i]) & (np.arange(n) != i))[0]
            other = np.where(codes != codes[i])[0]
            if same.size and other.size:
                sw = s[i, same].mean()
                sb = s[i, other].mean()
                N[i] = (sw - sb) / (sw + sb + eps)
        props.append(np.nanmean(N <= 0))
    return np.array(props)


def main():
    args = parse_args()
    adata_path = Path(args.adata)
    base_dir = Path(args.clustering_result)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figs_dir = out_dir / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    logging.info("Loading AnnData %s", adata_path)
    adata = sc.read_h5ad(str(adata_path))

    if "img_feats" not in adata.obsm:
        logging.error("adata.obsm does not contain 'img_feats'. Exiting.")
        sys.exit(1)
    if "spatial" not in adata.obsm:
        logging.error("adata.obsm does not contain 'spatial'. Exiting.")
        sys.exit(1)

    # choose features
    if args.selected_feats:
        selected_feats = args.selected_feats
    else:
        selected_feats = DEFAULT_SELECTED_FEATS

    img_feats = adata.obsm["img_feats"]
    # if it's a DataFrame-like object
    if hasattr(img_feats, "columns"):
        feat_names = list(img_feats.columns)
        use_feats = [f for f in selected_feats if f in feat_names]
        if len(use_feats) == 0:
            logging.warning("None of the default selected_feats found in img_feats; using all available features.")
            Xmat = img_feats.values
        else:
            logging.info("Using %d selected image features (intersection)", len(use_feats))
            Xmat = img_feats[use_feats].values
    else:
        # assume numpy array
        logging.info("img_feats is array-like; will use all columns")
        Xmat = np.asarray(img_feats)

    # standardize, similarities
    X = StandardScaler().fit_transform(Xmat)
    sims = cosine_similarity(X)
    coords = np.asarray(adata.obsm["spatial"])
    dists = pairwise_distances(coords)
    d_max = dists.max()

    # alpha values 0..1
    alphas = np.linspace(0, 1, args.bins)

    # find cluster files
    files = find_cluster_files(base_dir)
    logging.info("Found %d cluster files", len(files))

    records = []
    curves_map = {}  # key: workflow_id -> props array

    # compute ground truth first
    if args.ground_truth not in adata.obs:
        logging.error("Ground truth column '%s' not found in adata.obs. Exiting.", args.ground_truth)
        sys.exit(1)
    gt_series = adata.obs[args.ground_truth].astype(str)
    gt_codes = pd.Categorical(gt_series).codes
    logging.info("Computing ground truth curve")
    gt_props = alpha_curve_for_codes(gt_codes, sims, dists, d_max, alphas)

    for fp in tqdm(files, desc="Processing cluster files"):
        rel = None
        try:
            rel = fp.relative_to(base_dir)
            parts = rel.parts
            # expect structure like Normalization/FeatureSelection/Clustering/filename
            if len(parts) >= 4:
                normalization, feat_sel, clustering = parts[0], parts[1], parts[2]
            elif len(parts) >= 3:
                normalization, feat_sel, clustering = parts[0], parts[1], parts[2]
            else:
                normalization, feat_sel, clustering = "unknown", "unknown", parts[-2] if len(parts)>1 else fp.stem
        except Exception:
            normalization, feat_sel, clustering = "unknown", "unknown", fp.parent.name

        s = read_cluster_csv(fp)
        if s is None:
            continue

        # align to adata.obs index
        idx_inter = s.index.intersection(adata.obs.index)
        if len(idx_inter) / adata.n_obs < args.min_overlap:
            logging.warning("File %s overlaps only %.1f%% of adata barcodes; skipping", fp, 100 * len(idx_inter) / adata.n_obs)
            continue

        # produce codes for overlapping spots only; to integrate with sims/dists which are for full adata, we will create codes array of length n with -1 for missing
        codes = np.full(adata.n_obs, -1, dtype=int)
        common_idx = adata.obs.index.get_indexer(idx_inter)
        # convert labels to categorical using order in s.loc[idx_inter]
        cat = pd.Categorical(s.loc[idx_inter].astype(str))
        codes_common = cat.codes
        codes[common_idx] = codes_common

        # if too many -1 (missing), we already checked overlap, but ensure codes has at least some valid
        valid_mask = codes != -1
        if valid_mask.sum() == 0:
            logging.warning("No valid codes after alignment for %s", fp)
            continue

        # For alpha_curve function we expect codes for all n, but it handles categories; we'll map -1 as its own category but it will have no 'same' neighbors for those entries so okay
        # However better to transform to a contiguous code set where missing are unique
        cat_full = pd.Categorical(codes)
        codes_full = cat_full.codes

        props = alpha_curve_for_codes(codes_full, sims, dists, d_max, alphas)

        workflow_id = f"{normalization}__{feat_sel}__{clustering}__{fp.name}"
        curves_map[workflow_id] = {
            "normalization": normalization,
            "feature_selection": feat_sel,
            "clustering": clustering,
            "file": str(fp),
            "props": props,
        }

    # prepare CSV
    if len(curves_map) == 0:
        logging.error("No valid clustering curves computed. Exiting.")
        sys.exit(1)

    rows = []
    alpha_cols = [f"alpha_{int(a*100)}pct" for a in alphas]
    for wid, info in curves_map.items():
        row = {
            "workflow_id": wid,
            "normalization": info["normalization"],
            "feature_selection": info["feature_selection"],
            "clustering": info["clustering"],
            "file": info["file"],
        }
        for col, val in zip(alpha_cols, info["props"]):
            row[col] = float(val)
        rows.append(row)

    df_out = pd.DataFrame(rows)
    csv_out = out_dir / "alpha_curves.csv"
    df_out.to_csv(csv_out, index=False)
    logging.info("Saved curves CSV to %s", csv_out)

    # group by clustering and plot
    grouped = df_out.groupby("clustering")
    for clustering_name, group in grouped:
        plt.figure(figsize=(8, 5))
        # choose colormap
        cmap = plt.get_cmap("tab20")
        n = len(group)
        for i, (_, row) in enumerate(group.iterrows()):
            vals = row[alpha_cols].values.astype(float)
            plt.plot(np.linspace(0, 100, len(vals)), vals * 100, label=f"{row['normalization']}/{row['feature_selection']}", alpha=0.8, linewidth=1)
        # plot ground truth
        plt.plot(np.linspace(0, 100, len(gt_props)), gt_props * 100, label=f"GROUND_TRUTH ({args.ground_truth})", color="black", linewidth=3)
        plt.xlabel("α (0 = space-only, 1 = feature-only) [%]")
        plt.ylabel("Overall % of spots with N ≤ 0")
        plt.title(f"Clustering method: {clustering_name} (n={len(group)})")
        plt.legend(fontsize=8, loc="best")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        fig_fp = figs_dir / f"{clustering_name}_curves.png"
        plt.savefig(fig_fp, dpi=200)
        plt.close()
        logging.info("Saved figure %s", fig_fp)

    logging.info("Done. Figures in %s, CSV: %s", figs_dir, csv_out)


if __name__ == "__main__":
    main()

