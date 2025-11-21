#!/usr/bin/env python3
# scripts/prepare_anndata.py
# Usage examples:
# python scripts/prepare_anndata.py --mode raw --data_dir ./path/ --filename filtered_feature_bc_matrix.h5 --out data/SAMPLE/01_raw_anndata.h5ad --platform Visium
# python scripts/prepare_anndata.py --mode sc  --data_dir ./path/ --filename filtered_feature_bc_matrix.h5 --out data/SAMPLE/01_sc_raw_anndata.h5ad --platform Visium

import argparse, os, sys
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["raw","sc"], required=True,
                    help="raw -> read_10x_h5 and write 01_raw_anndata.h5ad; sc -> read_visium-like and write 01_sc_raw_anndata.h5ad")
parser.add_argument("--data_dir", required=True, help="Directory with 10x outputs (contains filtered_feature_bc_matrix.h5)")
parser.add_argument("--filename", default="filtered_feature_bc_matrix.h5")
parser.add_argument("--out", required=True, help="Output .h5ad path")
parser.add_argument("--platform", default="Visium", help="Platform: Visium, VisiumHD, Slideseq, Stereoseq")
args = parser.parse_args()

data_dir = args.data_dir
h5path = os.path.join(data_dir, args.filename)
outpath = args.out
mode = args.mode
platform = args.platform

os.makedirs(os.path.dirname(outpath), exist_ok=True)

if mode == "raw":
    if not os.path.exists(h5path):
        print("ERROR: H5 not found:", h5path, file=sys.stderr)
        sys.exit(2)
    print("Reading 10x h5 (raw) from:", h5path)
    adata = sc.read_10x_h5(h5path)
    adata.uns["source_10x_h5"] = h5path
    adata.write_h5ad(outpath)
    print("Wrote raw anndata to:", outpath)
    sys.exit(0)

# mode == "sc": read visium-like (scanpy) and patch coordinates (your snippet)
if mode == "sc":
    # prefer sc.read_visium if available
    try:
        print("Trying sc.read_visium on folder:", data_dir)
        # Depending on scanpy version signature - best effort:
        bdata = sc.read_visium(data_dir)
    except Exception as e:
        print("sc.read_visium failed:", e, " - trying sc.read_10x_h5 fallback.")
        if not os.path.exists(h5path):
            print("ERROR: H5 not found for fallback:", h5path, file=sys.stderr)
            sys.exit(3)
        bdata = sc.read_10x_h5(h5path)

    # Attempt the imagerow/imagecol fix you provided: look for pos table in data_dir if needed
    try:
        # Typical visium folder has 'spatial/tissue_positions_list.csv' or similar under data_dir
        # Try to locate positions file under <data_dir>/spatial or data_dir
        pos = None
        candidates = [
            os.path.join(data_dir, "spatial", "tissue_positions_list.csv"),
            os.path.join(data_dir, "spatial", "tissue_positions.csv"),
            os.path.join(data_dir, "spatial", "tissue_positions_list.tsv"),
            os.path.join(data_dir, "tissue_positions_list.csv"),
            os.path.join(data_dir, "tissue_positions.csv")
        ]
        for c in candidates:
            if os.path.exists(c):
                try:
                    pos = pd.read_csv(c, header=None)
                    # some files have header, try to read with header=0 as fallback
                    # but we will coerce to columns later
                except Exception:
                    pos = pd.read_csv(c, sep='\t', header=None)
                break

        # If pos not found via guessed files, try to extract from bdata.uns['spatial'] if present
        if pos is None:
            # some scanpy versions put positions in uns['spatial'][sample_id]['metadata'] or ['tissue_positions_list']
            if 'spatial' in bdata.uns and isinstance(bdata.uns['spatial'], dict):
                for k,v in bdata.uns['spatial'].items():
                    if isinstance(v, dict):
                        # try typical keys
                        for key in ("metadata", "tissue_positions_list", "tissue_positions"):
                            if key in v:
                                try:
                                    pos = pd.DataFrame(v[key])
                                    break
                                except Exception:
                                    pass
                    if pos is not None:
                        break

        # If still no pos, try bdata.obs for raw coords columns
        if pos is None and ('imagerow' in bdata.obs.columns or 'imagecol' in bdata.obs.columns):
            print("Found imagerow/imagecol in bdata.obs already.")
            bdata.obs['imagerow'] = bdata.obs.get('imagerow')
            bdata.obs['imagecol'] = bdata.obs.get('imagecol')
        elif pos is not None:
            # try to standardize pos to columns: barcode,in_tissue,array_row,array_col,imagerow,imagecol
            # your snippet expected first row to be header; attempt to detect
            if pos.shape[1] >= 6:
                # try common formatting: barcode,in_tissue,array_row,array_col,imagerow,imagecol
                if pos.shape[1] == 6:
                    pos.columns = ['barcode','in_tissue','array_row','array_col','imagerow','imagecol']
                else:
                    # if file includes header already, try to read with header=0
                    pos_try = pd.read_csv(c)
                    if set(['barcode','imagerow','imagecol']).issubset(pos_try.columns):
                        pos = pos_try
                    else:
                        # fallback name columns heuristically
                        pos.columns = ['barcode','in_tissue','array_row','array_col','imagerow','imagecol'][:pos.shape[1]]
                pos['barcode'] = pos['barcode'].astype(str).str.strip()
                pos['imagerow'] = pd.to_numeric(pos['imagerow'], errors='coerce')
                pos['imagecol'] = pd.to_numeric(pos['imagecol'], errors='coerce')
                pos = pos[~pos['barcode'].duplicated(keep='first')].set_index('barcode')
                idx = bdata.obs.index.astype(str)
                imagerow = pos['imagerow'].reindex(idx) if 'imagerow' in pos.columns else pd.Series(index=idx, dtype='float64')
                imagecol = pos['imagecol'].reindex(idx) if 'imagecol' in pos.columns else pd.Series(index=idx, dtype='float64')

                # handle -1 suffix cases as in your snippet
                if imagerow.isna().sum() > 0:
                    pos2 = pos.copy()
                    pos2.index = pos2.index.str.replace(r'-1$', '', regex=True)
                    idx2 = idx.str.replace(r'-1$', '', regex=True)
                    imagerow = imagerow.fillna(pd.Series(pos2['imagerow'].reindex(idx2).values, index=imagerow.index))
                    imagecol = imagecol.fillna(pd.Series(pos2['imagecol'].reindex(idx2).values, index=imagecol.index))

                bdata.obs['imagerow'] = pd.Series(imagerow.values, index=bdata.obs.index).astype('Int64')
                bdata.obs['imagecol'] = pd.Series(imagecol.values, index=bdata.obs.index).astype('Int64')
                print("Patched imagerow/imagecol into bdata.obs.")
            else:
                print("Position table present but unexpected shape; skipping patch.")
        else:
            print("No tissue_positions file found and no imagerow/imagecol in obs; continuing without patch.")

    except Exception as e:
        print("Warning: coordinate patch step failed:", str(e))

    bdata.write_h5ad(outpath)
    print("Wrote sc raw anndata (visium-like) to:", outpath)
    sys.exit(0)

# fallback
print("Unknown code path. Exiting.", file=sys.stderr)
sys.exit(5)

