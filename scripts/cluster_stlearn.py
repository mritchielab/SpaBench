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
import scipy.sparse as sp
from scipy.spatial import cKDTree

# norm_data.py lives next to this script; append its directory explicitly because
# `python -P` (used for SpatialLeiden) does not put the script directory on sys.path
import os as _os, sys as _sys
_sys.path.append(_os.path.dirname(_os.path.abspath(__file__)))
from norm_data import read_feature_list, apply_normalized  # noqa: E402

# SME weight matrices that need histology (morphology) features
MORPH_WEIGHTS = {"weights_matrix_all", "weights_matrix_pd_md",
                 "weights_matrix_gd_md", "morphological_distance"}
NO_MORPH_DEFAULT_WEIGHTS = "physical_distance"


def _spot_unit(adata, platform):
    """Spot-to-spot distance unit in imagecol/imagerow coordinates.

    stLearn (calculate_weight_matrix) regresses image coordinates on the Visium
    array coordinates and uses sqrt(coef_row^2 + coef_col^2) (~ one spot spacing);
    neighbours are spots closer than 3 units. For platforms without a real array
    grid (Slide-seq / Stereo-seq: array_* == image coords) the median
    nearest-neighbour distance is used as the unit instead.
    """
    from sklearn.linear_model import LinearRegression
    obs = adata.obs
    if (str(platform).lower() in ("visium", "visiumhd") and {"array_row", "array_col"} <= set(obs.columns)
            and not np.allclose(obs["array_row"].to_numpy(float), obs["imagerow"].to_numpy(float))):
        r = LinearRegression().fit(obs[["array_row"]].to_numpy(float), obs["imagerow"].to_numpy(float))
        c = LinearRegression().fit(obs[["array_col"]].to_numpy(float), obs["imagecol"].to_numpy(float))
        return float(np.sqrt(r.coef_[0] ** 2 + c.coef_[0] ** 2)), "array-coordinate regression (stLearn)"
    coor = obs[["imagecol", "imagerow"]].to_numpy(float)
    d, _ = cKDTree(coor).query(coor, k=2)
    return float(np.median(d[:, 1])), "median nearest-neighbour distance"


def sme_normalize_sparse(adata, weights, platform, rate=3.0, key_added="raw_SME_normalized"):
    """Memory-light re-implementation of stlearn SME_normalize(use_data='raw') for
    the histology-free weights 'physical_distance' and 'weights_matrix_pd_gd'.

    stLearn builds several dense n x n matrices (infeasible for Slide-seq /
    VisiumHD). Same rules as stLearn 0.4.x: neighbours = spots within 3 units;
    physical_distance -> equal weights over 6 neighbours (closest 6 here, stLearn
    takes an arbitrary 6 of the tied ones); pd_gd -> top-3 neighbours by Pearson
    correlation of X_pca. imputed = weighted mean of neighbours;
    SME = mean(X, imputed), entries where imputed == 0 keep X.
    """
    if weights not in ("physical_distance", "weights_matrix_pd_gd"):
        raise ValueError(f"SME weights '{weights}' are not supported by the sparse SME path; "
                         "use physical_distance or weights_matrix_pd_gd")
    coor = adata.obs[["imagecol", "imagerow"]].to_numpy(float)
    n = coor.shape[0]
    unit, how = _spot_unit(adata, platform)
    radius = rate * unit
    tree = cKDTree(coor)
    rows, cols, vals = [], [], []
    if weights == "physical_distance":
        d, idx = tree.query(coor, k=min(6, n), distance_upper_bound=radius)
        ok = np.isfinite(d)
        cnt = ok.sum(1)
        rows = np.repeat(np.arange(n), 6)[ok.ravel()]
        cols = idx[ok]
        vals = (1.0 / np.repeat(cnt, 6))[ok.ravel()]
    else:
        pcs = np.asarray(adata.obsm["X_pca"], dtype=float)
        pcs = pcs - pcs.mean(1, keepdims=True)
        pcs /= np.linalg.norm(pcs, axis=1, keepdims=True) + 1e-12
        rows, cols, vals = [], [], []
        for i, nb in enumerate(tree.query_ball_point(coor, radius)):
            nb = np.asarray(nb)
            w = np.clip(pcs[nb] @ pcs[i], 0, None)
            top = np.argsort(w)[-3:]
            nb, w = nb[top], w[top]
            if w.sum() > 0:
                rows.extend([i] * len(nb)); cols.extend(nb); vals.extend(w / w.sum())
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    X = X.astype(float)
    imputed = np.asarray(W @ X)
    adata.obsm[key_added] = np.where(imputed == 0, X, (X + imputed) / 2.0)
    print(f"SME (sparse): unit={unit:.3f} ({how}), radius={radius:.3f}, "
          f"mean neighbours used={W.getnnz() / n:.2f}", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-h5ad",       required=True,  dest="raw_h5ad")
    p.add_argument("--features",       required=True,  dest="feat_csv")
    p.add_argument("--norm-h5",        required=True,  dest="norm_h5",
                   help="Normalised expression of this workflow (02_<norm>_normalized.h5, "
                        "written by scripts/export_normalized.R)")
    p.add_argument("--out-h5ad",       required=True,  dest="out_h5ad")
    p.add_argument("--out-csv",        required=True,  dest="out_csv")
    p.add_argument("--out-rds",        required=True,  dest="out_rds")
    p.add_argument("--out-png",        required=True,  dest="out_png")
    p.add_argument("--use-morphology", action="store_true",       dest="use_morph")
    p.add_argument("--n-pcs",          type=int,    default=50,   dest="n_pcs")
    p.add_argument("--resolution",     type=float,  default=1.0,  dest="resolution")
    p.add_argument("--n-neighbors",    type=int,    default=10,   dest="n_neighbors")
    p.add_argument("--platform",       required=True, dest="platform")
    p.add_argument("--clustering-method", type=str, default="louvain",
                   choices=["louvain", "kmeans"], dest="clustering_method",
                   help="Clustering algorithm: louvain (default) or kmeans")
    p.add_argument("--n-clusters",    type=int,   default=None, dest="n_clusters",
                   help="Number of clusters (required when --clustering-method kmeans)")
    p.add_argument("--use-localization", action="store_true", dest="use_localization",
                   help="Apply DBSCAN spatial sub-clustering after initial clustering")
    p.add_argument("--localization-eps", type=float, default=20.0, dest="localization_eps",
                   help="DBSCAN epsilon radius for spatial localization (default: 20.0)")
    p.add_argument("--localization-min-samples", type=int, default=1,
                   dest="localization_min_samples",
                   help="DBSCAN min_samples for spatial localization (default: 1)")
    p.add_argument("--sme-weights", type=str, default="weights_matrix_all",
                   dest="sme_weights",
                   help="SME weight matrix type (default: weights_matrix_all; options: "
                        "weights_matrix_pd_gd, weights_matrix_pd_md, weights_matrix_gd_md, "
                        "physical_distance, morphological_distance, gene_expression_correlation)")
    p.add_argument("--sme-max-dense", type=int, default=15000, dest="sme_max_dense",
                   help="Use stLearn's own (dense n x n) SME_normalize for Visium samples up to "
                        "this many spots; larger samples / other platforms use the sparse "
                        "re-implementation (default: 15000)")
    args = p.parse_args()

    out_dir = Path(args.out_h5ad).parent
    tiling_dir = out_dir / "tiling"
    tiling_dir.mkdir(parents=True, exist_ok=True)
    print(f"stLearn tiling will write to {tiling_dir}", flush=True)

    # 1. Load raw AnnData (spot barcodes, obs coordinates, obsm['spatial'], uns)
    adata = sc.read(args.raw_h5ad)
    adata.var_names_make_unique()

    # 2. Replace X by THIS workflow's normalised values (Step 1) for ITS selected
    #    genes (Step 2). Spots not in the normalised object (e.g. extra Slide-seq
    #    beads) are dropped; obs/obsm/uns follow the reduced spot set.
    feats = read_feature_list(args.feat_csv)
    adata = apply_normalized(adata, args.norm_h5, feats)
   # 3. Convert for stLearn (adds spatial to .obsm) — ONLY for Visium / VisiumHD
    # If platform is Stereoseq or Slideseq we skip conversion and use adata as-is.
    if str(args.platform).lower() in ("visium", "visiumhd"):
        if "spatial" not in adata.uns:
            adata.uns["spatial"] = {
                "library": {
                    "scalefactors": {
                        "tissue_hires_scalef": 1.0,
                        "tissue_lowres_scalef": 1.0,
                    }
                }
            }
        adata = st.convert_scanpy(adata)
    else:
        # minimal change: skip convert_scanpy for Stereoseq/Slideseq
        print(f"Platform = {args.platform}; skipping st.convert_scanpy() and using input AnnData as-is.", flush=True)
        # Keep uns["spatial"] as a dict — assigning obsm["spatial"] (numpy array) would
        # break downstream calls that expect dict.keys() on uns["spatial"].
        if "spatial" not in adata.uns or not isinstance(adata.uns["spatial"], dict):
            adata.uns["spatial"] = {}

    # 4. Morphology (if requested)
    if args.use_morph:
        st.pp.tiling(adata, out_path=str(tiling_dir))
        st.pp.extract_feature(adata)

    # 5. Gene filtering only: X already holds the workflow's normalised values,
    #    so no normalize_total / log1p here.
    st.pp.filter_genes(adata, min_cells=1)

    # 6. PCA (needed by the gene-expression-correlation SME weights)
    n_pcs = min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    st.em.run_pca(adata, n_comps=n_pcs)

    # 7. SME normalisation (stLearn's spatial step) -- ALWAYS run.
    #    Without morphology, weights that need histology features are replaced
    #    by a histology-free option.
    weights = args.sme_weights
    if not args.use_morph and weights in MORPH_WEIGHTS:
        print(f"SME: morphology not used; weights '{weights}' need histology features -> "
              f"using '{NO_MORPH_DEFAULT_WEIGHTS}'", flush=True)
        weights = NO_MORPH_DEFAULT_WEIGHTS
    for k in ("imagerow", "imagecol"):
        if k not in adata.obs:
            adata.obs["imagecol"] = np.asarray(adata.obsm["spatial"])[:, 0]
            adata.obs["imagerow"] = np.asarray(adata.obsm["spatial"])[:, 1]
            break
    use_dense = (str(args.platform).lower() == "visium" and adata.n_obs <= args.sme_max_dense
                 and {"array_row", "array_col"} <= set(adata.obs.columns))
    print(f"SME: weights={weights}, use_morphology={args.use_morph}, "
          f"implementation={'stlearn.spatial.SME.SME_normalize' if use_dense else 'sparse re-implementation'}, "
          f"input = workflow-normalised X ({adata.n_obs} spots x {adata.n_vars} genes)", flush=True)
    if use_dense:
        dummy_morph = "X_morphology" not in adata.obsm
        if dummy_morph:
            # calculate_weight_matrix always computes a morphological distance;
            # a constant feature gives md == 1 and is unused by the chosen weights.
            adata.obsm["X_morphology"] = np.ones((adata.n_obs, 1))
        if not sp.issparse(adata.X):
            adata.X = np.asarray(adata.X)
        else:
            adata.X = sp.csr_matrix(adata.X)
        st.spatial.SME.SME_normalize(adata, use_data="raw", weights=weights, platform="Visium")
        if dummy_morph:
            del adata.obsm["X_morphology"]
        # drop the dense n x n matrices stLearn leaves in uns (huge, not needed)
        for k in ("gene_expression_correlation", "physical_distance", "morphological_distance",
                  "weights_matrix_all", "weights_matrix_pd_gd", "weights_matrix_pd_md",
                  "weights_matrix_gd_md"):
            adata.uns.pop(k, None)
        for k in ("imputed_data", "top_weights"):
            adata.obsm.pop(k, None)
    else:
        sme_normalize_sparse(adata, weights, args.platform)
    adata.uns["SME_weights"] = weights
    adata.X = adata.obsm.pop("raw_SME_normalized")
    st.pp.scale(adata)
    st.em.run_pca(adata, n_comps=n_pcs)
    print("SME done; PCA recomputed on SME-normalised data.", flush=True)

    # 8. Neighborhood graph and clustering
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=args.n_neighbors)

    if args.clustering_method == "kmeans":
        if args.n_clusters is None:
            raise ValueError("--n-clusters is required when --clustering-method kmeans")
        st.tl.clustering.kmeans(adata, n_clusters=args.n_clusters)
        cluster_key = "kmeans"
    else:
        sc.tl.louvain(adata, resolution=args.resolution)
        cluster_key = "louvain"

    # 9. Optional DBSCAN spatial sub-clustering (localization)
    if args.use_localization:
        st.spatial.clustering.localization(
            adata, use_label=cluster_key,
            eps=args.localization_eps,
            min_samples=args.localization_min_samples
        )
        loc_key = f"{cluster_key}_sub"
        if loc_key in adata.obs.columns:
            cluster_key = loc_key
        else:
            print(f"Warning: localization key '{loc_key}' not found; keeping '{cluster_key}' labels.", flush=True)

    # 10. Save results
    adata.obs['stlearn'] = adata.obs[cluster_key].astype(str)
    adata.write_h5ad(args.out_h5ad)
    adata.obs[['stlearn']].to_csv(args.out_csv)

    # 11. Plot spatial clusters
    _lib = list(adata.uns.get("spatial", {}).keys())
    _has_img = bool(_lib) and "images" in adata.uns["spatial"].get(_lib[0], {})
    if str(args.platform).lower() == "visium" and _has_img:
        fig = sc.pl.spatial(
            adata,
            color='stlearn',
            title='stLearn_louvain',
            show=False,
            return_fig=True,
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
