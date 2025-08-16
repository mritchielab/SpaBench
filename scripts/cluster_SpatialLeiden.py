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
            feat_names = feats.iloc[:,0].astype(str).tolist()
        feat_present = [g for g in feat_names if g in adata.var_names]
        if len(feat_present) == 0:
            warnings.warn("No selected features found in adata.var_names; proceeding without subsetting.")
        else:
            adata = adata[:, feat_present].copy()
            print(f"Subset to {len(feat_present)} features.", file=sys.stderr)

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
        # continue; spatialleiden may need spatial connectivity; we try to proceed

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
                latent_kwargs={},  # no random state
                spatial_kwargs={
                    "layer_ratio": args.layer_ratio,
                    "directed": (False, True),
                    "spatial_resolution": None,  # search_resolution will vary this
                },
            )
            print(f"Found latent_resolution={latent_resolution:.3f}, spatial_resolution={spatial_resolution:.3f}", file=sys.stderr)
        except Exception as e:
            warnings.warn(f"search_resolution failed, falling back to provided resolution {args.resolution}: {e}")
            latent_resolution = args.resolution
            spatial_resolution = args.resolution

    # run standard leiden for comparison
    try:
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
        # prefer sq.pl.spatial if available coordinates; fallback to sc.pl.embedding
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

    # create a simple RDS placeholder by reading the csv in R and saving an RDS
    # This produces an R-readable .rds file so downstream R stages that expect .rds succeed.
    try:
        cmd = [
            "Rscript", "-e",
            (
                "d<-read.csv(file='{csv}', stringsAsFactors=FALSE); "
                "saveRDS(d, file='{rds}')"
            ).format(csv=args.out_csv.replace("'", "\\'"), rds=args.out_rds.replace("'", "\\'"))
        ]
        subprocess.run(cmd, check=True)
        print(f"Wrote RDS placeholder from CSV -> {args.out_rds}", file=sys.stderr)
    except Exception as e:
        warnings.warn(f"Creating RDS via Rscript failed: {e}")

if __name__ == "__main__":
    main()
