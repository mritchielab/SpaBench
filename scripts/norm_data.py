"""Load the normalised expression written by scripts/export_normalized.R.

Used by cluster_SpaGCN.py, cluster_stlearn.py and cluster_SpatialLeiden.py so that
each workflow clusters the values produced by ITS normalisation method (Step 1),
restricted to the genes chosen by ITS feature-selection method (Step 2).
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _str(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def read_normalized(path):
    """Return (matrix spots x genes as CSR, gene names, spot barcodes, assay name)."""
    with h5py.File(path, "r") as f:
        n_genes, n_spots = (int(v) for v in f["shape"][:])
        m = sp.csc_matrix((f["data"][:], f["indices"][:], f["indptr"][:]), shape=(n_genes, n_spots))
        genes, barcodes = _str(f["genes"][:]), _str(f["barcodes"][:])
        assay = _str(np.atleast_1d(f["assay"][()]))[0] if "assay" in f else "unknown"
    return m.T.tocsr(), genes, barcodes, assay


def read_feature_list(path):
    """Gene list written by the R feature-selection scripts (one column with a header)."""
    return pd.read_csv(path, header=0, dtype=str).iloc[:, 0].dropna().tolist()


def apply_normalized(adata, norm_h5, features):
    """Replace adata.X by the workflow's normalised values for the selected genes.

    * spots: those present in BOTH adata and the normalised object, in adata order
      (so every method of a workflow clusters the same spots as the R methods)
    * genes: the selected features found in the normalised object, in list order
    obs / obsm (coordinates) of adata are kept. Raises if nothing matches.
    """
    import anndata as ad

    X, genes, barcodes, assay = read_normalized(norm_h5)
    spot_pos, gene_pos = {}, {}
    for i, b in enumerate(barcodes):
        spot_pos.setdefault(b, i)
    for i, g in enumerate(genes):
        gene_pos.setdefault(g, i)
    keep_spots = [b for b in adata.obs_names if b in spot_pos]
    if not keep_spots:
        raise ValueError("No spot barcodes shared between the AnnData and the normalised object")
    feats = [g for g in dict.fromkeys(features) if g in gene_pos]
    if not feats:
        raise ValueError("None of the selected features are in the normalised object")
    missing = len(set(features)) - len(feats)

    rows = np.array([spot_pos[b] for b in keep_spots])
    cols = np.array([gene_pos[g] for g in feats])
    sub = adata[keep_spots, :].copy()
    out = ad.AnnData(X=X[rows][:, cols].tocsr().astype(np.float32),
                     obs=sub.obs.copy(), var=pd.DataFrame(index=pd.Index(feats)))
    for k, v in sub.obsm.items():
        out.obsm[k] = v
    out.uns = dict(sub.uns)
    print(f"Normalised data ({assay}): {out.n_obs} spots x {out.n_vars} selected genes "
          f"({adata.n_obs - out.n_obs} spots not in the normalised object dropped; "
          f"{missing} selected genes not found)", flush=True)
    return out
