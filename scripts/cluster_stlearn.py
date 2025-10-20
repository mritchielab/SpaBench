#!/usr/bin/env python3
import argparse
import scanpy as sc
import stlearn as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
from matplotlib.figure import Figure
from matplotlib.axes import Axes

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
    p.add_argument("--platform",       required=True, dest="platform")
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
   # 3. Convert for stLearn (adds spatial to .obsm) — ONLY for Visium / VisiumHD
    # If platform is Stereoseq or Slideseq we skip conversion and use adata as-is.
    if str(args.platform).lower() in ("visium", "visiumhd"):
        adata = st.convert_scanpy(adata)
    else:
        # minimal change: skip convert_scanpy for Stereoseq/Slideseq
        print(f"Platform = {args.platform}; skipping st.convert_scanpy() and using input AnnData as-is.", flush=True)
        adata.uns["spatial"] = adata.obsm["spatial"]

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
    if str(args.platform).lower() in ("visium", "visiumhd"):
      fig = sc.pl.spatial(
        adata,
        color='stlearn',
        title='stLearn_louvain',
        show=False,
        return_fig=True
    )    
    else:
      fig = sc.pl.embedding(adata, basis="spatial", color="stlearn",
                      title="stLearn_louvain", show=False, return_fig=True)
    
        # Robust save helper: handle Figure, Axes, list/tuple/ndarray of them, or other duck-typed objects.
    def _save_plot_object(obj, out_path):
        """
        Save a plot object to out_path.
        Accepts matplotlib.figure.Figure, matplotlib.axes.Axes, list/tuple/ndarray of those,
        or other objects supporting savefig/get_figure.
        If multiple objects are given, create numbered files with _0 _1 suffixes.
        """
        p = Path(out_path)
        # If list-like, iterate and save with suffixes
        if isinstance(obj, (list, tuple, np.ndarray)):
            # if only one element, save to original path
            if len(obj) == 1:
                _save_plot_object(obj[0], out_path)
            else:
                stem = p.stem
                ext = p.suffix or ".png"
                parent = p.parent
                for i, sub in enumerate(obj):
                    target = parent / f"{stem}_{i}{ext}"
                    _save_plot_object(sub, str(target))
            return

        # If Figure, save directly
        if isinstance(obj, Figure):
            obj.savefig(str(p), dpi=300, bbox_inches='tight')
            plt.close(obj)
            return

        # If Axes, get its figure then save
        if isinstance(obj, Axes):
            fig = obj.get_figure()
            fig.savefig(str(p), dpi=300, bbox_inches='tight')
            plt.close(fig)
            return

        # duck-typing fallback: try savefig or get_figure
        if hasattr(obj, "savefig") and callable(getattr(obj, "savefig")):
            obj.savefig(str(p), dpi=300, bbox_inches='tight')
            try:
                plt.close(obj)
            except Exception:
                pass
            return
        if hasattr(obj, "get_figure") and callable(getattr(obj, "get_figure")):
            fig = obj.get_figure()
            fig.savefig(str(p), dpi=300, bbox_inches='tight')
            plt.close(fig)
            return

        raise RuntimeError(f"Cannot save plot object of type {type(obj)} to {out_path}")

    # Save the returned figure/axes robustly. If multiple subplots/figs are returned,
    # generate numbered files: out.png -> out_0.png, out_1.png, ...
    _save_plot_object(fig, args.out_png)
    
    # fig can be single Figure or list，所以统一处理
    #if isinstance(fig, list):
    #    for f in fig:
    #        f.savefig(args.out_png, dpi=300, bbox_inches='tight')
    #        plt.close(f)
    #else:
    #    fig.savefig(args.out_png, dpi=300, bbox_inches='tight')
    #    plt.close(fig)
        
    # 12. export Seurat RDS via H5AD + Rscript (using schard::h5ad2seurat)
    tmp_h5ad = args.out_rds.replace(".rds", ".h5ad")
    adata.write_h5ad(tmp_h5ad)
    print(f"Temporary H5AD written to {tmp_h5ad}", flush=True)

    # Rscript transfer and save as Seurat RDS
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
