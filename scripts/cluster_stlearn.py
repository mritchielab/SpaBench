#!/usr/bin/env python3
import argparse
import scanpy as sc
import stlearn as st
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-h5ad",       required=True,  dest="raw_h5ad")
    p.add_argument("--features",       required=True,  dest="feat_csv")
    p.add_argument("--out-h5ad",       required=True,  dest="out_h5ad")
    p.add_argument("--out-csv",        required=True,  dest="out_csv")
    p.add_argument("--out-rds",        required=True,  dest="out_rds")
    p.add_argument("--out-png",        required=True,  dest="out_png")
    p.add_argument("--use-morphology", action="store_true",       dest="use_morph")
    p.add_argument("--n-pcs",          type=int,    default=50,   dest="n_pcs")
    p.add_argument("--resolution",     type=float,  default=1.0,  dest="resolution")
    p.add_argument("--n-neighbors",    type=int,    default=10,   dest="n_neighbors")
    args = p.parse_args()

    out_dir = Path(args.out_h5ad).parent
    tiling_dir = out_dir / "tiling"
    tiling_dir.mkdir(parents=True, exist_ok=True)
    print(f"stLearn tiling will write to {tiling_dir}", flush=True)

    # 1. Load raw AnnData
    adata = sc.read(args.raw_h5ad)
    adata.var_names_make_unique()

    # 2. Subset to selected features
    feats = pd.read_csv(args.feat_csv, header=0).iloc[:,0].astype(str).tolist()
    valid_feats = [g for g in feats if g in adata.var_names]

    bad = set(feats) - set(valid_feats)
    if bad:
      print(f"Warning: {len(bad)} genes not found in AnnData and will be skipped:", bad, flush=True)

    adata = adata[:, valid_feats].copy()
    # 3. Convert for stLearn (adds spatial to .obsm)
    adata = st.convert_scanpy(adata)

    # 4. Morphology (if requested)
    if args.use_morph:
        st.pp.tiling(adata, out_path=str(tiling_dir))
        st.pp.extract_feature(adata)

    # 5. Gene preprocessing
    st.pp.filter_genes(adata, min_cells=1)
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)

    # 6. PCA
    st.em.run_pca(adata, n_comps=args.n_pcs)

    # 7. SME normalize (only if morphology is on)
    if args.use_morph:
        st.spatial.SME.SME_normalize(adata, use_data="raw")
        adata.X = adata.obsm["raw_SME_normalized"]
        st.pp.scale(adata)
        st.em.run_pca(adata, n_comps=args.n_pcs)

    # 8. Neighborhood graph and clustering
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=args.n_neighbors)
    sc.tl.louvain(adata, resolution=args.resolution)

    # 10. Save results
    adata.obs['stlearn'] = adata.obs['louvain'].astype(str)
    adata.write_h5ad(args.out_h5ad)
    adata.obs[['stlearn']].to_csv(args.out_csv)

    # 11. Plot spatial clusters
    fig = sc.pl.spatial(
        adata,
        color='stlearn',
        title='stLearn_louvain',
        show=False,
        return_fig=True
    )
    # fig 有可能是单个 Figure 或 list，所以统一处理
    if isinstance(fig, list):
        for f in fig:
            f.savefig(args.out_png, dpi=300, bbox_inches='tight')
            plt.close(f)
    else:
        fig.savefig(args.out_png, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
    # 12. export Seurat RDS via H5AD + Rscript (using schard::h5ad2seurat)
    # 构造临时 h5ad 路径（与 out_rds 只替换后缀）
    tmp_h5ad = args.out_rds.replace(".rds", ".h5ad")
    # 先把当前 AnnData 再写一次为临时 h5ad
    adata.write_h5ad(tmp_h5ad)
    print(f"Temporary H5AD written to {tmp_h5ad}", flush=True)

    # 调用 Rscript 转换并保存为 Seurat RDS
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
        print(f"Seurat RDS saved to {args.out_rds}", flush=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Rscript failed with exit code {e.returncode}")

if __name__ == "__main__":
    main()
