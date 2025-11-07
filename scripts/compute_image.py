#!/usr/bin/env python3
"""
Compute image features with Squidpy and store them in adata.obsm["img_feats"].

Usage example:
  python compute_image.py --adata /path/to/your.h5ad --image /path/to/your.tiff --out_dir /path/to/output/adata/with/imagefeature

Comments in this file are in English. The user-facing notes are in Chinese in the chat.
"""

import os
import argparse
import logging
from pathlib import Path

import scanpy as sc
import squidpy as sq
import anndata
import pandas as pd


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("compute_image_features")


def ensure_out_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)


def move_features_to_obsm(adata: anndata.AnnData, key: str):
    """Ensure that image features end up in adata.obsm[key].

    Squidpy may store results in adata.uns or adata.obsm depending on versions/options.
    This helper tries common locations and moves/sets the result into adata.obsm[key]
    as a pandas.DataFrame (preferred) or numpy array.
    """
    # If already in obsm, nothing to do
    if key in adata.obsm:
        logger.info(f"Features already present in adata.obsm['{key}']")
        return

    # Common fallback: some versions write to uns
    possible_uns_keys = [key, "img_features", "img_feats", "image_features"]
    for k in possible_uns_keys:
        if k in adata.uns:
            df = adata.uns.pop(k)
            # If it's a DataFrame, preserve index/columns
            if isinstance(df, pd.DataFrame):
                adata.obsm[key] = df
            else:
                # try to coerce to DataFrame using adata.obs_names
                try:
                    df2 = pd.DataFrame(df, index=adata.obs_names)
                    adata.obsm[key] = df2
                except Exception:
                    # last resort: place raw object into obsm
                    adata.obsm[key] = df
            logger.info(f"Moved features from adata.uns['{k}'] to adata.obsm['{key}']")
            return

    # If we reach here, we couldn't find results
    raise KeyError(f"Could not find image features in adata after calculate_image_features. Searched for '{key}' and common fallbacks.")


def main():
    parser = argparse.ArgumentParser(description="Compute image features with Squidpy and save them into adata.obsm['img_feats'].")
    parser.add_argument("--adata", required=True, help="Path to input AnnData (.h5ad)")
    parser.add_argument("--image", required=True, help="Path to high-resolution image (tiff)")
    parser.add_argument("--out_dir", required=True, help="Directory to save output AnnData with image features")
    parser.add_argument("--library_id", default=None, help="Optional library id or comma-separated list for multi-library AnnData")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of parallel jobs to use")
    parser.add_argument("--key", default="img_feats", help="Key to store features in adata.obsm (default: img_feats)")
    args = parser.parse_args()

    adata_path = args.adata
    image_path = args.image
    out_dir = args.out_dir
    lib = args.library_id
    n_jobs = args.n_jobs
    key = args.key

    ensure_out_dir(out_dir)

    logger.info(f"Reading AnnData from {adata_path}")
    adata = sc.read_h5ad(adata_path)

    logger.info(f"Loading image from {image_path} into Squidpy ImageContainer")
    img = sq.im.ImageContainer(image_path)

    # Define the feature set we want to compute. This attempts to cover the common
    # Squidpy feature extractors. If one fails (e.g., segmentation not available),
    # we try again without it.
    feature_sets = ["summary", "texture", "histogram"]

    try:
        logger.info(f"Calculating image features: {feature_sets}")
        # call Squidpy's high-level function
        sq.im.calculate_image_features(
            adata,
            img,
            library_id=(lib.split(",") if lib else None),
            features=feature_sets,
            key_added=key,
            n_jobs=n_jobs,
            show_progress_bar=True,
        )
    except ValueError as e:
        # common case: segmentation features require a segmentation mask which may not exist
        msg = str(e).lower()
        if "segmentation" in msg or "mask" in msg:
            logger.warning("Segmentation features failed, retrying without segmentation: %s", e)
            reduced = [f for f in feature_sets if f != "segmentation"]
            sq.im.calculate_image_features(
                adata,
                img,
                library_id=(lib.split(",") if lib else None),
                features=reduced,
                key_added=key,
                n_jobs=n_jobs,
                show_progress_bar=True,
            )
        else:
            logger.exception("calculate_image_features raised an unexpected ValueError")
            raise

    # Ensure results are in adata.obsm[key]
    try:
        move_features_to_obsm(adata, key)
    except KeyError as e:
        logger.exception("Failed to locate image features after calculation: %s", e)
        raise

    # Save the AnnData with a descriptive filename
    base = Path(adata_path).stem
    out_path = Path(out_dir) / f"{base}_with_image_features.h5ad"
    logger.info(f"Writing annotated AnnData with image features to {out_path}")
    adata.write(out_path)
    logger.info("Done.")


if __name__ == "__main__":
    main()

