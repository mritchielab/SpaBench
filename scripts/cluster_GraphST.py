#!/usr/bin/env python3
import os
# **二次保险**：把 R_HOME 和 PATH 再次设对
os.environ['R_HOME'] = '/stornext/System/data/software/rhel/9/base/tools/R/4.5.1/lib64/R'
os.environ['PATH'] = os.environ['R_HOME'] + '/bin:' + os.environ.get('PATH', '')

import argparse
import scanpy as sc
import pandas as pd
import torch
from GraphST import GraphST
from GraphST.utils import clustering as gs_clustering

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw_h5ad",    required=True)
    p.add_argument("--positions",   required=True)
    p.add_argument("--features",    required=True)
    p.add_argument("--n_clusters",  type=int, required=True)
    p.add_argument("--radius",      type=int, required=True)
    p.add_argument("--out_rds",     required=True)
    p.add_argument("--out_csv",     required=True)
    p.add_argument("--out_png",     required=True)
    args = p.parse_args()

    # 1. 读原始 anndata
    adata = sc.read_h5ad(args.raw_h5ad)
    # 2. 读坐标文件 (无 header，第一列是 barcode)
    pos = pd.read_csv(args.positions, header=None, index_col=0)
    # 模仿 tutorial 里的列名映射
    adata.obs['array_row'] = pos[1].astype(int)
    adata.obs['array_col'] = pos[2].astype(int)
    # scanpy visium 需要放到 obsm['spatial']
    adata.obsm['spatial'] = adata.obs[['array_row','array_col']].values

    # 3. 只保留特征选择结果里那几列 gene
    feat_df = pd.read_csv(args.features, stringsAsFactors=False)
    genes = feat_df[[feat_df.columns[0]]].iloc[:,0].tolist()
    adata = adata[:, genes].copy()

    # 4. 构建模型 & 训练
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = GraphST.GraphST(adata, device=device)
    model.train()

    # 5. 聚类 (mclust 推荐，带 refinement)
    #    GraphST tutorial 里 “mclust” 用 refinement=True
    from GraphST.utils import clustering
    adata = clustering(
        adata,
        n_clusters=args.n_clusters,
        radius=args.radius,
        method='mclust',
        refinement=True
    )

    # 6. 保存结果
    #    domain 存在 adata.obs['domain']
    out_meta = adata.obs[['domain']].copy()
    out_meta.to_csv(args.out_csv)
    #    保存为 h5ad 以便复查
    adata.write_h5ad(args.out_rds.replace('.rds','.h5ad'))
    #    RDS：可以只存 metadata
    import anndata2ri
    from rpy2.robjects import r, pandas2ri
    pandas2ri.activate()
    anndata2ri.activate()
    # 转到 R 的 SeuratObject
    r(f'''
      library(Seurat)
      df <- read.csv("{args.out_csv}", row.names=1)
      seu <- CreateSeuratObject(counts = matrix(0, nrow=0, ncol=0)) # 占位
      seu@meta.data <- df
      saveRDS(seu, file="{args.out_rds}")
    ''')
    # 7. spatial plot
    import matplotlib.pyplot as plt
    sc.pl.spatial(adata, color='domain', title='GraphST domains',
                  show=False, save=args.out_png.replace('.png',''))
    # sc.pl.spatial 会把文件名写到 fig目录下，–save 只要 basename
    
if __name__ == "__main__":
    main()
