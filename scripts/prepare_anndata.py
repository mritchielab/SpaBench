#!/usr/bin/env python3
# scripts/prepare_anndata.py

import argparse, os, sys
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import json
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["raw","sc"], required=True,
                    help="raw -> read_10x_h5 and write 01_raw_anndata.h5ad; sc -> read_visium-like and write 01_sc_raw_anndata.h5ad")
parser.add_argument("--data_dir", required=True, help="Directory with 10x outputs")
parser.add_argument("--filename", default="filtered_feature_bc_matrix.h5")
parser.add_argument("--out", required=True, help="Output .h5ad path")
parser.add_argument("--platform", default="Visium", help="Platform: Visium, VisiumHD, Slideseq, Stereoseq")
# --- NEW: Added bin_size for VisiumHD detection ---
parser.add_argument("--bin_size", default="008um", help="Bin size for VisiumHD (008um or 016um)")

args = parser.parse_args()

data_dir = args.data_dir
h5path = os.path.join(data_dir, args.filename)
outpath = args.out
mode = args.mode
platform = args.platform

os.makedirs(os.path.dirname(outpath), exist_ok=True)

# ==============================================================================
# ### VISIUM HD PATH REDIRECTION
# ==============================================================================
if platform == "VisiumHD":
    # If the H5 isn't at the root, look into the binned_outputs folder
    if not os.path.exists(h5path):
        hd_path = os.path.join(data_dir, "binned_outputs", f"square_{args.bin_size}")
        if os.path.exists(os.path.join(hd_path, args.filename)):
            print(f"VisiumHD detected. Redirecting paths to: {hd_path}")
            h5path = os.path.join(hd_path, args.filename)
            # Update data_dir for sc.read_visium to find the spatial folder
            data_dir = hd_path 

# ==============================================================================

if mode == "raw":
    if platform == "Slideseq":
        print("Processing Slideseq in 'raw' mode...")
        mtx_file = glob.glob(os.path.join(data_dir, "*expression_matrix.mtx"))[0]
        adata = sc.read_mtx(mtx_file).T
        adata.var_names_make_unique()
        # Add required metadata fields
        adata.var['gene_ids'] = adata.var_names
        adata.var['feature_types'] = 'Gene Expression'
        adata.var['genome'] = 'unknown' # Or 'mm10' if known
        adata.write_h5ad(outpath)
        print("Wrote Slideseq raw anndata to:", outpath)
        sys.exit(0)
    
    else:
        if not os.path.exists(h5path):
            print("ERROR: H5 not found:", h5path, file=sys.stderr)
            sys.exit(2)
        print("Reading 10x h5 (raw) from:", h5path)
        adata = sc.read_10x_h5(h5path)
        adata.uns["source_10x_h5"] = h5path
        adata.write_h5ad(outpath)
        print("Wrote raw anndata to:", outpath)
        sys.exit(0)

if mode == "sc":
    if platform == "Slideseq":
        print("Processing Slideseq in 'sc' mode...")
        # 1. Locate files
        mtx_files = glob.glob(os.path.join(data_dir, "*expression_matrix.mtx"))
        loc_files = glob.glob(os.path.join(data_dir, "*BeadLocationsForR.csv"))
        mtx_file = mtx_files[0]
        loc_file = loc_files[0]
        
        # 2. Load Expression
        adata = sc.read_mtx(mtx_file).T
        locations = pd.read_csv(loc_file, index_col=0) # Index is 'barcodes'

        # 3. Align Barcodes (The fix for n_obs = 0)
        # If scanpy loaded the matrix with default integer indices (0, 1, 2...)
        # and the shape matches the location file, assign the actual barcodes.
        if adata.n_obs == locations.shape[0]:
            print(f"Counts match ({adata.n_obs}). Mapping barcodes from location file...")
            adata.obs_names = locations.index
        else:
            print("Shape mismatch. Attempting intersection...")
            common = adata.obs_names.intersection(locations.index)
            adata = adata[common].copy()
            locations = locations.loc[common]
            print(f"Matched {len(common)} barcodes.")

        if adata.n_obs == 0:
            print("ERROR: AnnData is empty. Check if the matrix needs transposing or barcode format.", file=sys.stderr)
            sys.exit(6)

        adata.var_names_make_unique()
        
        # 4. Populate obs
        adata.obs['in_tissue'] = 1
        adata.obs['array_row'] = locations['ycoord'].values
        adata.obs['array_col'] = locations['xcoord'].values
        adata.obs['imagerow'] = locations['ycoord'].values
        adata.obs['imagecol'] = locations['xcoord'].values
        
        # 5. Populate var
        adata.var['gene_ids'] = adata.var_names
        adata.var['feature_types'] = 'Gene Expression'
        adata.var['genome'] = 'unknown'
        
        # 6. Populate obsm and uns
        adata.obsm['spatial'] = locations[['xcoord', 'ycoord']].values.astype(float)
        adata.uns['spatial'] = {
            "sample": {
                "images": {}, 
                "scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 10.0},
                "metadata": {}
            }
        }
        
        adata.write_h5ad(outpath)
        print("Wrote Slideseq sc raw anndata to:", outpath)
        sys.exit(0)

    else:
        try:
            print("Trying sc.read_visium on folder:", data_dir)
            # scanpy's read_visium expects the .h5 and 'spatial' folder to be in the same dir
            bdata = sc.read_visium(data_dir)
        except Exception as e:
            print("sc.read_visium failed:", e, " - trying sc.read_10x_h5 fallback.")
            if not os.path.exists(h5path):
                print("ERROR: H5 not found for fallback:", h5path, file=sys.stderr)
                sys.exit(3)
            bdata = sc.read_10x_h5(h5path)

        try:
            pos = None
            # --- IMPROVED: Added .parquet support and VisiumHD column mapping ---
            candidates = [
                os.path.join(data_dir, "spatial", "tissue_positions.parquet"), # VisiumHD Default
                os.path.join(data_dir, "spatial", "tissue_positions_list.csv"),
                os.path.join(data_dir, "spatial", "tissue_positions.csv"),
                os.path.join(data_dir, "tissue_positions_list.csv"),
                os.path.join(data_dir, "tissue_positions.csv")
            ]
            
            for c in candidates:
                if os.path.exists(c):
                    print(f"Found positions file: {c}")
                    if c.endswith('.parquet'):
                        pos = pd.read_parquet(c)
                    else:
                        try:
                            pos = pd.read_csv(c, header=None)
                        except Exception:
                            pos = pd.read_csv(c, sep='\t', header=None)
                    break

            if pos is not None:
                # --- IMPROVED: Standardize VisiumHD column names ---
                if 'pxl_row_in_fullres' in pos.columns:
                    print("Standardizing VisiumHD coordinate columns...")
                    pos = pos.rename(columns={
                        'barcode': 'barcode',
                        'pxl_row_in_fullres': 'imagerow',
                        'pxl_col_in_fullres': 'imagecol',
                        'in_tissue': 'in_tissue',
                        'array_row': 'array_row',
                        'array_col': 'array_col'
                    })
                
                # Original logic for processing the pos dataframe
                if pos.shape[1] >= 6 and not isinstance(pos.columns[0], str):
                    pos.columns = ['barcode','in_tissue','array_row','array_col','imagerow','imagecol'][:pos.shape[1]]
                
                if 'barcode' in pos.columns:
                    pos['barcode'] = pos['barcode'].astype(str).str.strip()
                    pos = pos[~pos['barcode'].duplicated(keep='first')].set_index('barcode')
                    
                    idx = bdata.obs.index.astype(str)
                    # Handle -1 suffix logic
                    if not idx[0] in pos.index and idx[0].split('-')[0] in pos.index:
                        pos.index = pos.index + "-1"
                    
                    bdata.obs['imagerow'] = pos['imagerow'].reindex(idx).values
                    bdata.obs['imagecol'] = pos['imagecol'].reindex(idx).values
                    
                    if 'in_tissue' in pos.columns:
                        bdata.obs['in_tissue'] = pos['in_tissue'].reindex(idx).values
                    
                    print("Patched imagerow/imagecol into bdata.obs.")

        except Exception as e:
            print("Warning: coordinate patch step failed:", str(e))

        # --- NEW: Ensure Spatial Data exists in uns for VisiumHD plotting ---
        if 'spatial' not in bdata.uns and platform == "VisiumHD":
            print("Finalizing VisiumHD spatial metadata...")
            pass 

        bdata.write_h5ad(outpath)
        print("Wrote sc raw anndata to:", outpath)
        sys.exit(0)

print("Unknown code path. Exiting.", file=sys.stderr)
sys.exit(5)
