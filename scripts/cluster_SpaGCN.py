#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, random, warnings
import scanpy as sc
import pandas as pd
import numpy as np
import SpaGCN as spg
import matplotlib.pyplot as plt
import cv2
import subprocess
import torch
# norm_data.py lives next to this script; append its directory explicitly because
# `python -P` (used for SpatialLeiden) does not put the script directory on sys.path
import os as _os, sys as _sys
_sys.path.append(_os.path.dirname(_os.path.abspath(__file__)))
from norm_data import read_feature_list, apply_normalized  # noqa: E402

warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_h5ad",   required=True, help="原始 AnnData (.h5ad)")
    parser.add_argument("--positions",  required=True, help="无表头位置 CSV")
    parser.add_argument("--features",   required=True, help="feature list CSV")
    parser.add_argument("--norm_h5",    required=True,
                        help="Normalised expression of this workflow (02_<norm>_normalized.h5, "
                             "written by scripts/export_normalized.R)")
    parser.add_argument("--out_rds",    required=True, help="输出 Seurat RDS 路径")
    parser.add_argument("--out_csv",    required=True, help="输出 cluster CSV 路径")
    parser.add_argument("--out_png",    required=True, help="输出空间图 PNG 路径")
    
    # Flags
    parser.add_argument("--histology",  action="store_true", help="是否使用 histology")
    parser.add_argument("--refine",     action="store_true", help="是否执行 cluster refine")
    
    # Added parameters to match Config/Tutorial
    parser.add_argument("--img",        type=str, default=None, help="Histology image path (required if histology=True)")
    parser.add_argument("--beta",       type=int, default=49,   help="Beta value for histology (default 49)")
    parser.add_argument("--alpha",      type=float, default=1.0, help="Alpha value for histology (default 1)")
    parser.add_argument("--n_clusters", type=int, default=7,    help="Number of clusters (default 7)")
    parser.add_argument("--num_pcs",    type=int, default=50,   help="PCs computed inside SpaGCN.train (package default 50)")
    parser.add_argument("--platform",   type=str, default="Visium", help="Visium (hexagonal grid) or other (square grid)")
    parser.add_argument("--p_res",      type=float, default=0.5, help="p value for searching l (default 0.5)")

    args = parser.parse_args()

    # 1. 读 AnnData (raw; only used for spot barcodes / obs)
    adata = sc.read_h5ad(args.raw_h5ad)
    adata.var_names_make_unique()

    # 2. Replace X by THIS workflow's normalised values for ITS selected genes
    #    (spots restricted to those present in the normalised object).
    #    No further library-size normalisation / log1p is applied below.
    feats = read_feature_list(args.features)
    adata = apply_normalized(adata, args.norm_h5, feats)

    # 3. 读无表头位置文件 & Map Coordinates (after the spot subset, so they line up)
    # 假设列结构: 0=barcode, 1=flag, 2=x_array, 3=y_array, 4=x_pixel, 5=y_pixel
    pos = pd.read_csv(args.positions, header=None, index_col=0)
    pos.index = pos.index.astype(str) # Ensure index is string type
    missing_pos = adata.obs_names.difference(pos.index)
    if len(missing_pos):
        raise ValueError(f"{len(missing_pos)} spots have no entry in {args.positions}")

    # Load Pixel Coordinates (for training)
    adata.obs["x_pixel"] = pos[4].loc[adata.obs_names].values
    adata.obs["y_pixel"] = pos[5].loc[adata.obs_names].values

    # Load Array Coordinates (for refinement - critical per tutorial)
    adata.obs["x_array"] = pos[2].loc[adata.obs_names].values
    adata.obs["y_array"] = pos[3].loc[adata.obs_names].values

    # 4. 计算邻接矩阵 (Training Adjacency)
    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()
    
    # Logic: If histology is requested AND image is provided, use image + beta + alpha
    if args.histology and args.img and os.path.exists(args.img):
        print(f"Calculating adjacency with histology image: {args.img}")
        img = cv2.imread(args.img)
        adj = spg.calculate_adj_matrix(
            x=x_pixel, 
            y=y_pixel,
            x_pixel=x_pixel, 
            y_pixel=y_pixel,
            image=img, 
            beta=args.beta, 
            alpha=args.alpha, 
            histology=True
        )
    else:
        print("Calculating adjacency without histology (coordinates only)")
        adj = spg.calculate_adj_matrix(
            x=x_pixel, 
            y=y_pixel, 
            histology=False
        )

    # 5. No prefilter / normalize_total / log1p: adata.X already holds the
    #    workflow's normalised values (Step 1) for its selected genes (Step 2).

    # 6. 搜索 l (using p from args)
    l = spg.search_l(args.p_res, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    # 7. 搜索 res (using n_clusters from args)
    res = spg.search_res(
        adata, adj, l, args.n_clusters,
        start=0.7, step=0.1, tol=5e-3, lr=0.05,
        max_epochs=20, r_seed=100, t_seed=100, n_seed=100
    )

    # 8. 运行 SpaGCN
    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(100); torch.manual_seed(100); np.random.seed(100)
    clf.train(
        adata, adj,
        num_pcs=min(args.num_pcs, adata.n_vars - 1, adata.n_obs - 1),
        init_spa=True, init="louvain",
        res=res, tol=5e-3, lr=0.05, max_epochs=200
    )
    y_pred, _ = clf.predict()
    adata.obs["pred"] = y_pred.astype(str)
    
    # Default final column is 'pred' unless refined
    final_col = "pred"

    # 9. Refinement (Logic improved based on Tutorial)
    if args.refine:
        print("Performing Cluster Refinement...")
        # Tutorial: Refinement uses 2D adjacency based on array coordinates (not pixel)
        # and histology=False
        x_array = adata.obs["x_array"].tolist()
        y_array = adata.obs["y_array"].tolist()
        
        adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
        
        refined = spg.refine(
            sample_id=adata.obs_names.tolist(),
            pred=adata.obs["pred"].tolist(),
            dis=adj_2d,     # Uses the array-coordinate specific adjacency
            shape="hexagon" if args.platform.lower() == "visium" else "square"  # Visium spots are hexagonal; HD / Stereo-seq bins are square
        )
        adata.obs["refined_pred"] = [str(r) for r in refined]
        final_col = "refined_pred"

    # 10. save cluster CSV
    pd.DataFrame({
        "barcode": adata.obs_names,
        "cluster": adata.obs[final_col]
    }).to_csv(args.out_csv, index=False)

    # 11. plotting
    # Set colors manually or let Scanpy handle it. 
    # Tutorial logic for domains and colors is handled automatically by sc.pl.scatter
    sc.pl.scatter(
        adata,
        x="x_pixel", y="y_pixel",
        color=final_col,
        size=100000/adata.shape[0], # Adjusted size logic from tutorial
        show=False,
        title=f"SpaGCN Clusters ({final_col})"
    )
    plt.gca().invert_yaxis()
    plt.savefig(args.out_png, dpi=300)
    plt.close()

    # 12. export Seurat RDS
    tmp_h5ad = args.out_rds.replace(".rds", ".h5ad")
    adata.write_h5ad(tmp_h5ad)
    print(f"Temporary H5AD written to {tmp_h5ad}")
    try:
        subprocess.check_call([
            "Rscript", "-e",
            (
                "library(reticulate); "
                "library(Seurat); "
                "library(SeuratObject); "
                "library(schard); "
                f"seu_data <- h5ad2seurat('{tmp_h5ad}', assay='Spatial'); "
                f"SaveSeuratRds(seu_data, file = '{args.out_rds}')"
            )
        ])
        print(f"Seurat RDS saved to {args.out_rds}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Rscript failed with exit code {e.returncode}")

if __name__ == "__main__":
    main()
