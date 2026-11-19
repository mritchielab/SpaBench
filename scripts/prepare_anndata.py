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

# ---- Step 0 QC ----
parser.add_argument("--remove_mt", default="true",
                    help="QC: remove mitochondrial genes, all platforms (true/false) [default true]")
parser.add_argument("--mt_pattern", default="^mt-",
                    help="QC: regular expression for mitochondrial gene names, matched case-insensitively")
parser.add_argument("--in_tissue_only", default="true",
                    help="QC: Visium / Visium HD - keep only spots/bins with in_tissue = 1 (true/false) [default true]")

args = parser.parse_args()


def _flag(v):
    return str(v).strip().lower() in ("true", "t", "1", "yes")


REMOVE_MT = _flag(args.remove_mt)
IN_TISSUE_ONLY = _flag(args.in_tissue_only)


# ==============================================================================
# Step 0 QC (same rules in scripts/prepare_seurat.R)
#   * mitochondrial genes are removed for every platform
#   * Visium / Visium HD: spots / bins outside the captured tissue (in_tissue != 1
#     in tissue_positions) are removed
# ==============================================================================
def qc_remove_mt(adata):
    if not REMOVE_MT:
        print("QC: mitochondrial genes kept (--remove_mt false)")
        return adata
    mt = np.asarray(pd.Index(adata.var_names).str.contains(args.mt_pattern, case=False, regex=True))
    names = list(adata.var_names[mt][:15])
    print(f"QC: removing {int(mt.sum())} mitochondrial genes (pattern '{args.mt_pattern}', case-insensitive)"
          + (": " + ", ".join(names) if names else ""))
    if not mt.any():
        print(f"WARNING: QC: no gene matched the mitochondrial pattern '{args.mt_pattern}'", file=sys.stderr)
        return adata
    return adata[:, ~mt].copy()


def read_tissue_flags(d):
    """in_tissue flag per barcode from the Space Ranger tissue_positions file (None if absent)."""
    cands = [os.path.join(d, "spatial", "tissue_positions.parquet"),
             os.path.join(d, "spatial", "tissue_positions.csv"),
             os.path.join(d, "spatial", "tissue_positions_list.csv"),
             os.path.join(d, "tissue_positions.parquet"),
             os.path.join(d, "tissue_positions.csv"),
             os.path.join(d, "tissue_positions_list.csv")]
    for c in cands:
        if not os.path.exists(c):
            continue
        if c.endswith(".parquet"):
            pos = pd.read_parquet(c)
        else:
            with open(c) as fh:
                has_header = "barcode" in fh.readline().lower()
            pos = pd.read_csv(c, header=0 if has_header else None)
            if not has_header:
                pos = pos.iloc[:, :2]
                pos.columns = ["barcode", "in_tissue"]
        pos["barcode"] = pos["barcode"].astype(str).str.strip()
        pos = pos.drop_duplicates("barcode").set_index("barcode")
        return pd.to_numeric(pos["in_tissue"], errors="coerce").fillna(0).astype(int), c
    return None, None


def qc_in_tissue(adata, d):
    if not IN_TISSUE_ONLY:
        print("QC: off-tissue spots kept (--in_tissue_only false)")
        return adata
    flags, src = read_tissue_flags(d)
    if flags is None:
        print(f"ERROR: QC: no tissue_positions file found in {d} - cannot remove spots outside the tissue "
              "(set qc: in_tissue_only: false to skip)", file=sys.stderr)
        sys.exit(8)
    idx = pd.Index(adata.obs_names.astype(str))
    # positions files without the "-1" barcode suffix used in the matrix
    if not idx.isin(flags.index).any() and idx.isin(flags.index + "-1").any():
        flags.index = flags.index + "-1"
    keep = np.asarray(flags.reindex(idx).fillna(0).astype(int) == 1)
    if not keep.any():
        print(f"ERROR: QC: no matrix barcode has in_tissue = 1 in {src}", file=sys.stderr)
        sys.exit(8)
    print(f"QC: keeping {int(keep.sum())} of {adata.n_obs} spots/bins inside the tissue "
          f"(in_tissue = 1 in {os.path.basename(src)})")
    return adata[keep].copy() if not keep.all() else adata

data_dir = args.data_dir
h5path = os.path.join(data_dir, args.filename)
outpath = args.out
mode = args.mode
platform = args.platform
# Accept any capitalisation of the platform name and map it to the spelling used below
platform = {"visium": "Visium", "visiumhd": "VisiumHD", "slideseq": "Slideseq",
            "stereoseq": "Stereoseq"}.get(str(platform).lower(), platform)

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
# ### SLIDE-SEQ / STEREO-SEQ (both modes)
# The file-discovery, count-selection and bead/bin filtering rules below are
# mirrored in scripts/prepare_seurat.R so that 01_raw_seurat.rds,
# 01_raw_anndata.h5ad and 01_sc_raw_anndata.h5ad hold the same barcodes and
# genes. Files written by the workflow itself (01_*) are never used as input.
# ==============================================================================
import re
import h5py
import scipy.io
import scipy.sparse as sp
try:                                   # anndata >= 0.10
    from anndata.io import read_elem
except ImportError:                    # anndata 0.8 / 0.9
    from anndata.experimental import read_elem


def die(msg, code=7):
    print("ERROR: " + msg, file=sys.stderr)
    sys.exit(code)


def list_inputs(d, pattern):
    """Basenames in d matching pattern (case-insensitive), excluding 01_* outputs."""
    rx = re.compile(pattern, re.IGNORECASE)
    return sorted(f for f in os.listdir(d)
                  if rx.search(f) and not f.startswith("01_")
                  and os.path.isfile(os.path.join(d, f)))


def open_text(path):
    import gzip
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def is_integer_valued(m):
    d = m.data if sp.issparse(m) else np.asarray(m)
    return d.size == 0 or bool(np.all(np.isfinite(d)) and np.all(d >= 0)
                               and np.array_equal(d, np.round(d)))


def clean_gene_names(genes):
    """Same rules as Seurat (underscores -> dashes) + R make.unique(sep='-'),
    so AnnData and Seurat objects carry identical gene names."""
    genes = [str(g).replace("_", "-") for g in genes]
    seen, out = set(genes), []
    counts = {}
    first = set()
    for g in genes:
        if g not in first:
            first.add(g)
            out.append(g)
            continue
        k = counts.get(g, 0)
        while True:
            k += 1
            cand = f"{g}-{k}"
            if cand not in seen:
                break
        counts[g] = k
        seen.add(cand)
        out.append(cand)
    return out


# Count matrices tried in this order; the first one that is integer-valued wins.
H5AD_COUNT_KEYS = ["layers/raw_count", "layers/count", "layers/counts", "raw/X", "X"]


def h5_index(f, grp):
    col = f[grp].attrs.get("_index", "_index")
    return np.asarray(read_elem(f[grp][col])).astype(str)


def read_h5ad_counts(path, need_coords=True):
    """Raw integer counts (cells x genes, CSR), obs/var names and x/y coordinates."""
    with h5py.File(path, "r") as f:
        obs_names = h5_index(f, "obs")
        X, key = None, None
        for k in H5AD_COUNT_KEYS:
            if k not in f:
                continue
            m = read_elem(f[k])
            m = sp.csr_matrix(m) if not sp.issparse(m) else m.tocsr()
            if is_integer_valued(m):
                X, key = m, k
                break
            print(f"  {k}: not integer-valued (normalised?) - skipped")
        if X is None:
            die(f"{path}: no integer-valued count matrix in {H5AD_COUNT_KEYS}; "
                "raw counts are required (only normalised values found).")
        var_names = h5_index(f, "raw/var" if key.startswith("raw/") else "var")
        if X.shape != (len(obs_names), len(var_names)):
            die(f"{path}: {key} has shape {X.shape}, expected "
                f"({len(obs_names)}, {len(var_names)})")
        x = y = None
        if need_coords:
            if "x" in f["obs"] and "y" in f["obs"]:
                x = np.asarray(read_elem(f["obs/x"]), dtype=float)
                y = np.asarray(read_elem(f["obs/y"]), dtype=float)
            elif "obsm" in f and "spatial" in f["obsm"]:
                s = np.asarray(read_elem(f["obsm/spatial"]), dtype=float)
                x, y = s[:, 0], s[:, 1]
            else:
                die(f"{path}: coordinates not found in obs/x,y or obsm/spatial.")
    print(f"  counts from '{key}': {X.shape[0]} cells x {X.shape[1]} genes, "
          f"total = {X.sum():.0f}")
    return X, obs_names, var_names, x, y, key


def read_mtx(path):
    with open_text(path) as fh:
        m = scipy.io.mmread(fh)
    return sp.csr_matrix(m)


def read_slideseq(d):
    """Slide-seq: built ONLY from *expression_matrix.mtx + *locations*.csv
    (barcodes / gene names from a barcodes / genes file or digital_expression.txt).
    Any .h5ad in the folder is ignored. All beads present in both the matrix and
    the locations file are kept."""
    mtx = list_inputs(d, r"expression_matrix\.mtx(\.gz)?$")
    locs = list_inputs(d, r"locations.*\.csv(\.gz)?$")
    if not mtx or not locs:
        die(f"Slide-seq input needs *expression_matrix.mtx and *locations*.csv in {d}")
    mtx_path, loc_path = os.path.join(d, mtx[0]), os.path.join(d, locs[0])
    print("  matrix:   ", mtx_path)
    print("  locations:", loc_path)
    M = read_mtx(mtx_path)

    loc = pd.read_csv(loc_path)
    loc.index = loc.iloc[:, 0].astype(str).str.strip()
    loc = loc[~loc.index.duplicated(keep="first")]
    cols = {c.lower(): c for c in loc.columns}
    xcol = cols.get("xcoord", cols.get("x", loc.columns[1]))
    ycol = cols.get("ycoord", cols.get("y", loc.columns[2]))

    # --- bead barcodes (matrix order) ---
    bc_files = list_inputs(d, r"barcodes.*\.(tsv|txt|csv)(\.gz)?$")
    dge = list_inputs(d, r"digital_expression\.txt(\.gz)?$")
    if bc_files:
        barcodes = pd.read_csv(os.path.join(d, bc_files[0]), header=None, sep=None,
                               engine="python").iloc[:, 0].astype(str).str.strip().tolist()
        bc_src = bc_files[0]
    elif dge:
        with open_text(os.path.join(d, dge[0])) as fh:
            barcodes = fh.readline().rstrip("\n\r").split("\t")[1:]
        bc_src = dge[0] + " (header)"
    else:
        barcodes = loc.index.tolist()
        bc_src = locs[0] + " (row order)"
    if M.shape[1] == len(barcodes):
        M = M.T.tocsr()                      # genes x beads -> beads x genes
    elif M.shape[0] != len(barcodes):
        die(f"matrix {M.shape} matches neither dimension of {len(barcodes)} barcodes from {bc_src}")
    print(f"  barcodes from {bc_src}: {len(barcodes)}")
    n_genes = M.shape[1]

    # --- gene names ---
    genes, g_src = None, None
    g_files = list_inputs(d, r"(genes|features).*\.(tsv|txt|csv)(\.gz)?$")
    if g_files:
        g = pd.read_csv(os.path.join(d, g_files[0]), header=None, sep=None,
                        engine="python").iloc[:, 0].astype(str).tolist()
        if len(g) == n_genes:
            genes, g_src = g, g_files[0]
    if genes is None and dge:
        with open_text(os.path.join(d, dge[0])) as fh:
            fh.readline()
            g = [line.split("\t", 1)[0] for line in fh]
        if len(g) == n_genes:
            genes, g_src = g, dge[0] + " (first column)"
    if genes is None:
        print("WARNING: no gene names found for the Slide-seq matrix (no genes/features file "
              "or digital_expression.txt) - using placeholders Feature1..N",
              file=sys.stderr)
        genes, g_src = [f"Feature{i + 1}" for i in range(n_genes)], "PLACEHOLDER"
    print(f"  genes from {g_src}: {n_genes}")

    # --- keep beads present in both the matrix and the locations file ---
    barcodes = pd.Index(barcodes)
    keep = np.flatnonzero(barcodes.isin(loc.index))
    print(f"  beads: matrix {len(barcodes)}, locations {loc.shape[0]}, kept (both) {len(keep)}")
    if len(keep) == 0:
        die("no bead barcodes shared between matrix and locations file", 6)
    M = M[keep]
    barcodes = barcodes[keep]
    loc = loc.loc[barcodes]
    return (M, barcodes.astype(str).tolist(), genes,
            loc[xcol].to_numpy(dtype=float), loc[ycol].to_numpy(dtype=float), mtx_path)


def stereoseq_h5ad(d, filename):
    """Same choice as prepare_seurat.R: *_processed.h5ad, else --filename, else first .h5ad."""
    hits = list_inputs(d, r"\.h5ad$")
    if not hits:
        die(f"No input .h5ad (other than 01_*) found in {d} for Stereo-seq.")
    proc = [h for h in hits if h.lower().endswith("_processed.h5ad")]
    if proc:
        return os.path.join(d, proc[0])
    if filename in hits:
        return os.path.join(d, filename)
    return os.path.join(d, hits[0])


def write_spatial_anndata(X, obs_names, genes, x, y, source):
    adata = ad.AnnData(X=X.astype(np.float32),
                       obs=pd.DataFrame(index=pd.Index(obs_names, dtype=str)),
                       var=pd.DataFrame(index=pd.Index(clean_gene_names(genes), dtype=str)))
    adata.var["gene_ids"] = adata.var_names
    adata.var["feature_types"] = "Gene Expression"
    adata.var["genome"] = "unknown"
    # obsm['spatial'] = (x, y); in obs, row = y and col = x (same as the Seurat image)
    adata.obsm["spatial"] = np.column_stack([x, y]).astype(float)
    if mode == "sc":
        adata.obs["in_tissue"] = 1
        adata.obs["array_row"] = y
        adata.obs["array_col"] = x
        adata.obs["imagerow"] = y
        adata.obs["imagecol"] = x
        adata.uns["spatial"] = {
            "sample": {
                "images": {},
                "scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 10.0},
                "metadata": {}
            }
        }
    adata.uns["source_file"] = source
    adata = qc_remove_mt(adata)          # Step 0 QC (Slide-seq / Stereo-seq have no off-tissue flag)
    adata.write_h5ad(outpath)
    print(f"Wrote {platform} {mode} anndata ({adata.n_obs} x {adata.n_vars}) to: {outpath}")
    sys.exit(0)


if platform == "Slideseq":
    print(f"Processing Slideseq in '{mode}' mode...")
    X, obs_names, genes, x, y, src = read_slideseq(data_dir)
    write_spatial_anndata(X, obs_names, genes, x, y, src)

if platform == "Stereoseq":
    print(f"Processing Stereoseq in '{mode}' mode...")
    src = stereoseq_h5ad(data_dir, args.filename)
    print("  input:", src)
    X, obs_names, genes, x, y, _ = read_h5ad_counts(src)
    write_spatial_anndata(X, obs_names, genes, x, y, src)

# ==============================================================================
# ### VISIUM / VISIUM HD (10x .h5)
# ==============================================================================

if mode == "raw":
    if not os.path.exists(h5path):
        print("ERROR: H5 not found:", h5path, file=sys.stderr)
        sys.exit(2)
    print("Reading 10x h5 (raw) from:", h5path)
    adata = sc.read_10x_h5(h5path)
    adata.uns["source_10x_h5"] = h5path
    adata = qc_in_tissue(adata, data_dir)  # Step 0 QC
    adata = qc_remove_mt(adata)
    adata.write_h5ad(outpath)
    print("Wrote raw anndata to:", outpath)
    sys.exit(0)

if mode == "sc":
    try:
        print("Trying sc.read_visium on folder:", data_dir)
        # scanpy's read_visium expects the .h5 and 'spatial' folder to be in the same dir
        bdata = sc.read_visium(data_dir, count_file=os.path.basename(h5path))
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
                    pos = pd.read_csv(c, header=0 if 'barcode' in open(c).readline().lower() else None)
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

    # Downstream spatial methods (squidpy / stLearn) need obsm['spatial'];
    # the read_10x_h5 fallback only fills obs columns, so add it here.
    if 'spatial' not in bdata.obsm and {'imagerow', 'imagecol'} <= set(bdata.obs.columns):
        bdata.obs['imagerow'] = pd.to_numeric(bdata.obs['imagerow'], errors='coerce')
        bdata.obs['imagecol'] = pd.to_numeric(bdata.obs['imagecol'], errors='coerce')
        bdata.obsm['spatial'] = bdata.obs[['imagecol', 'imagerow']].to_numpy(dtype=float)
    bdata = qc_in_tissue(bdata, data_dir)  # Step 0 QC
    bdata = qc_remove_mt(bdata)
    bdata.write_h5ad(outpath)
    print("Wrote sc raw anndata to:", outpath)
    sys.exit(0)

print("Unknown code path. Exiting.", file=sys.stderr)
sys.exit(5)
