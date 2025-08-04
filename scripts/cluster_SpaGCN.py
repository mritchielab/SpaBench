#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, random, torch, warnings
import scanpy as sc
import pandas as pd
import numpy as np
import SpaGCN as spg
import matplotlib.pyplot as plt
import cv2
import subprocess
import torch
import random
import warnings

warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_h5ad",   required=True, help="原始 AnnData (.h5ad)")
    parser.add_argument("--positions",  required=True, help="无表头位置 CSV, 格式: barcode,flag,x1,x2,x3,x4,x5")
    parser.add_argument("--features",   required=True, help="feature list CSV, 一列基因名，无表头")
    parser.add_argument("--out_rds",    required=True, help="输出 Seurat RDS 路径")
    parser.add_argument("--out_csv",    required=True, help="输出 cluster CSV 路径")
    parser.add_argument("--out_png",    required=True, help="输出空间图 PNG 路径")
    parser.add_argument("--refine",     action="store_true", help="是否执行 cluster refine")
    args = parser.parse_args()

    # 1. 读 AnnData
    adata = sc.read_h5ad(args.raw_h5ad)
    adata.var_names_make_unique()

    # 2. 读无表头位置文件
    #    列: 0=barcode,1=flag,2=x1,3=x2,4=x3,5=x4,6=x5
    pos = pd.read_csv(args.positions, header=None, index_col=0)
    # Tutorial 映射逻辑：x_pixel<-col5, y_pixel<-col4
    adata.obs["y_pixel"] = pos[5].loc[adata.obs_names].values
    adata.obs["x_pixel"] = pos[4].loc[adata.obs_names].values

    # 3. 读 feature list 并 subset
    feats = pd.read_csv(args.features, header=None, dtype=str)[0]
    keep = [g for g in feats if g in adata.var_names]
    if not keep:
        raise ValueError("No features matched var_names")
    adata = adata[:, keep].copy()

    # 4. 计算邻接矩阵（不含 histology）
    x = adata.obs["x_pixel"].tolist()
    y = adata.obs["y_pixel"].tolist()
    adj = spg.calculate_adj_matrix(x=x, y=y, histology=False)

    # 5. 预过滤 & 归一化
    adata.var_names_make_unique()
    spg.prefilter_genes(adata, min_cells=3)
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    # 6. 搜索 l
    p_target = 0.5
    l = spg.search_l(p_target, adj,
                     start=0.01, end=1000,
                     tol=0.01, max_run=100)

    # 7. 搜索 res
    #    这里默认用 leiden 初始 res=1.0，如果你知道期望簇数可自定义
    n_clusters = int(adata.obs.shape[0]**0.5)  # 示例：√n spots
    res = spg.search_res(adata, adj, l, n_clusters,
                        start=0.5, step=0.1,
                        tol=5e-3, lr=0.05,
                        max_epochs=20,
                        r_seed=100, t_seed=100, n_seed=100)

    # 8. 运行 SpaGCN
    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(100); torch.manual_seed(100); np.random.seed(100)
    clf.train(adata, adj,
              init_spa=True, init="louvain",
              res=res, tol=5e-3,
              lr=0.05, max_epochs=200)
    y_pred, _ = clf.predict()
    adata.obs["pred"] = y_pred.astype(str)

    # 9. if refine
    final_col = "pred"
    if args.refine:
        adj2 = spg.calculate_adj_matrix(x=x, y=y, histology=False)
        refined = spg.refine(
            sample_id=adata.obs_names.tolist(),
            pred     =adata.obs["pred"].tolist(),
            dis      =adj2,
            shape    ="hexagon"
        )
        adata.obs["refined_pred"] = [str(r) for r in refined]
        final_col = "refined_pred"

    # 10. save cluster CSV
    pd.DataFrame({
        "barcode": adata.obs_names,
        "cluster": adata.obs[final_col]
    }).to_csv(args.out_csv, index=False)

    # 11. plotting
    sc.pl.scatter(adata,
                  x="x_pixel", y="y_pixel",
                  color=final_col,
                  size=23, show=False)
    plt.gca().invert_yaxis()
    plt.savefig(args.out_png, dpi=300)
    plt.close()

    # 12. export Seurat RDS via H5AD + Rscript (using schard::h5ad2seurat)
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

