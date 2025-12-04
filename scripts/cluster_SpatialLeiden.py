#!/usr/bin/env python3
"""
scripts/cluster_SpatialLeiden.py

Run SpatialLeiden clustering and optionally perform resolution search (refine).
Outputs:
 - h5ad (annotated AnnData)
 - clusters csv (obs index, cluster label)
 - PNG plot
 - RDS placeholder (created from the csv via a small Rscript)
"""

import argparse
import os
import sys
import subprocess
import warnings
import traceback
import random

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import squidpy as sq
import spatialleiden as sl

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-h5ad", dest="raw_h5ad", required=True)
    p.add_argument("--features", dest="features_csv", required=False)
    p.add_argument("--out-h5ad", dest="out_h5ad", required=True)
    p.add_argument("--out-csv", dest="out_csv", required=True)
    p.add_argument("--out-png", dest="out_png", required=True)
    p.add_argument("--out-rds", dest="out_rds", required=True)
    p.add_argument("--n-neighbors", dest="n_neighbors", type=int, default=10)
    p.add_argument("--n-pcs", dest="n_pcs", type=int, default=30)
    p.add_argument("--resolution", dest="resolution", type=float, default=1.0)
    p.add_argument("--layer-ratio", dest="layer_ratio", type=float, default=1.8)
    p.add_argument("--coord-type", dest="coord_type", type=str, default="generic")
    p.add_argument("--spatial-n-neighs", dest="spatial_n_neighs", type=int, default=10)
    p.add_argument("--refine", dest="refine", action="store_true")
    p.add_argument("--search-n-clusters", dest="search_n_clusters", type=int, default=None)
    return p.parse_args()

def main():
    args = parse_args()
    
    # Set global seeds for reproducibility since spatialleiden argument failed
    random_state = 0
    np.random.seed(random_state)
    random.seed(random_state)
    
    sc.settings.autosave = False
    sc.settings.verbosity = 3

    # load anndata
    adata = ad.read_h5ad(args.raw_h5ad)
    print(f"Loaded AnnData: {adata.shape}", file=sys.stderr)

    # subset to features if provided (expect a CSV with a column of gene names)
    if args.features_csv and os.path.exists(args.features_csv):
        feats = pd.read_csv(args.features_csv, header=0)

        # try common column names
        if "gene" in feats.columns:
            feat_names = feats["gene"].astype(str).tolist()
        elif "feature" in feats.columns:
            feat_names = feats["feature"].astype(str).tolist()
        else:
            # take first column
            feat_names = feats.iloc[:, 0].astype(str).tolist()

        # remove duplicates in feat_names but preserve order (keep first occurrence)
        seen = set()
        feat_names_unique = []
        for g in feat_names:
            if g not in seen:
                feat_names_unique.append(g)
                seen.add(g)

        # Map gene names to the first matching position in adata.var_names
        positions = []
        # matched_names = [] # unused variable
        for g in feat_names_unique:
            # find indices where var_names equals g
            pos = np.where(adata.var_names == g)[0]
            if pos.size > 0:
                idx = int(pos[0])              # take the first match only
                if idx not in positions:       # avoid adding same index twice
                    positions.append(idx)
                    # matched_names.append(adata.var_names[idx])

        if len(positions) == 0:
            warnings.warn("No selected features found in adata.var_names; proceeding without subsetting.")
        else:
            # Subset by integer positions (safe even when var_names are not unique)
            adata = adata[:, positions].copy()
            print(f"Subset to {len(positions)} features (duplicates removed/first-match kept).", file=sys.stderr)

    # preprocess: log1p and PCA
    if "X_pca" not in adata.obsm:
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=args.n_pcs)
        print("Computed PCA.", file=sys.stderr)

    # compute latent neighbors (based on PCA)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X_pca")
    print("Computed latent neighbors.", file=sys.stderr)

    # build spatial neighbors
    # squidpy expects spatial coordinates in adata.obsm['spatial'] OR adata.obs with 'x','y'.
    try:
        sq.gr.spatial_neighbors(adata, coord_type=args.coord_type, n_neighs=args.spatial_n_neighs)
        print("Computed spatial neighbors via squidpy.", file=sys.stderr)
    except Exception as e:
        warnings.warn(f"Failed to compute spatial neighbors: {e}")

    # convert distances -> spatial_connectivities if spatial distances exist
    if "spatial_distances" in adata.obsp:
        try:
            adata.obsp["spatial_connectivities"] = sl.distance2connectivity(adata.obsp["spatial_distances"])
            print("Converted spatial distances -> spatial_connectivities.", file=sys.stderr)
        except Exception as e:
            warnings.warn(f"distance2connectivity failed: {e}")

    # optionally run resolution search to match desired n_clusters
    latent_resolution = args.resolution
    spatial_resolution = args.resolution
    
    if args.refine and (args.search_n_clusters is not None):
        n_clusters = int(args.search_n_clusters)
        print(f"Searching resolution to reach {n_clusters} clusters...", file=sys.stderr)

        try:
            latent_resolution, spatial_resolution = sl.search_resolution(
                adata,
                n_clusters,
                latent_kwargs={"random_state": random_state},
                spatial_kwargs={
                    "layer_ratio": args.layer_ratio,
                    # REMOVED "random_state": random_state because spatialleiden doesn't accept it
                    "directed": (False, True),
                },
            )
            print(f"Found latent_resolution={latent_resolution:.4f}, spatial_resolution={spatial_resolution:.4f}", file=sys.stderr)
        except Exception as e:
            print("search_resolution failed; printing traceback:", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            warnings.warn(f"search_resolution failed, falling back to provided resolution {args.resolution}: {e}")
            latent_resolution = args.resolution
            spatial_resolution = args.resolution


    # run standard leiden for comparison
    try:
        # Note: Added flavor='igraph' to silence future warning if your environment supports it
        sc.tl.leiden(adata, resolution=latent_resolution, key_added="leiden", directed=False)
    except Exception as e:
        warnings.warn(f"leiden failed: {e}")

    # run spatialleiden
    try:
        sl.spatialleiden(
            adata,
            layer_ratio=args.layer_ratio,
            resolution=(latent_resolution, spatial_resolution),
            directed=(False, True),
            key_added="spatialleiden"
        )
        print("spatialleiden finished.", file=sys.stderr)
    except Exception as e:
        print("spatialleiden failed:", e, file=sys.stderr)
        sys.exit(1)

    # write cluster labels csv (index, spatialleiden)
    if "spatialleiden" in adata.obs:
        df = adata.obs[["spatialleiden"]].copy()
        df["cell"] = df.index
        df = df.reset_index(drop=True)[["cell", "spatialleiden"]]
        df.to_csv(args.out_csv, index=False)
        print(f"Wrote clusters csv: {args.out_csv}", file=sys.stderr)
    else:
        warnings.warn("spatialleiden not present in adata.obs after run; writing empty csv.")
        pd.DataFrame(columns=["cell", "spatialleiden"]).to_csv(args.out_csv, index=False)

    # plot spatial embedding colored by clusters (use spatial coordinates)
    try:
        if "spatial" in adata.obsm_keys():
            sc.pl.embedding(adata, basis="spatial", color=["spatialleiden"], show=False)
            plt.savefig(args.out_png, bbox_inches="tight", dpi=150)
        else:
            # fallback: umap or pca embedding
            key = "X_pca" if "X_pca" in adata.obsm_keys() else "X_umap"
            sc.pl.embedding(adata, basis="pca" if key=="X_pca" else "umap", color=["spatialleiden"], show=False)
            plt.savefig(args.out_png, bbox_inches="tight", dpi=150)
        plt.close()
        print(f"Wrote cluster plot: {args.out_png}", file=sys.stderr)
    except Exception as e:
        warnings.warn(f"plotting failed: {e}")

    # write annotated anndata
    try:
        adata.write_h5ad(args.out_h5ad)
        print(f"Wrote annotated AnnData: {args.out_h5ad}", file=sys.stderr)
    except Exception as e:
        warnings.warn(f"Writing h5ad failed: {e}")

    # Call Rscript to save Seurat RDS
    try:
        subprocess.check_call([
            "Rscript", "-e",
            (
                "library(reticulate); "
                "library(Seurat); "
                "library(SeuratObject); "
                "library(schard); "
                f"seu_data <- h5ad2seurat('{args.out_h5ad}', assay='Spatial'); "
                f"SaveSeuratRds(seu_data, file = '{args.out_rds}')"
            )
        ])
        print(f"Seurat RDS saved to {args.out_rds}", flush=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Rscript failed with exit code {e.returncode}")

if __name__ == "__main__":
    main()
