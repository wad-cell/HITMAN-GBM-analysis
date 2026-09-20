# -*- coding: utf-8 -*-
"""Shared scRNA scoring utilities for C2 (GSE324656 project).
Scoring definition (frozen C2.6):
  AUCell/UCell-style rank-based score on DETECTED genes per cell:
    for cell c with n_detected genes, ranks over its non-zero expression (average ranks).
    k = # signature genes detected in cell; sig_rank = mean rank among them.
    score = (sig_rank - (k+1)/2) / (n_detected - k)   if n_detected>k else NaN
  This equals the Mann-Whitney U statistic scaled to [0,1] comparing signature
  expression ranks vs all other detected genes within the cell (UCell definition
  restricted to expressed genes; documented as such).
Alternative scoring (sensitivity):
  alt_score = mean( log1p(expr) over signature genes detected )  (per-cell; NaN if k=0)
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

def load_signatures(gene_tsv, removal_tsv=None, drop_dup=True):
    df = pd.read_csv(gene_tsv, sep="\t")
    sig = {}
    for name, g in df.groupby("signature"):
        syms = list(g["symbol"].dropna().unique())
        if drop_dup:
            syms = list(dict.fromkeys(syms))
        sig[name] = syms
    rem = None
    if removal_tsv is not None:
        rem = pd.read_csv(removal_tsv, sep="\t")
        rem_col = "symbol" if "symbol" in rem.columns else rem.columns[1]
        rem_set = set(rem[rem_col].dropna().astype(str))
        rem = rem_set
    return sig, rem

def make_variants(sig, rem=None):
    """Return dict of score-ready gene lists: original + CC-depleted (if rem given)."""
    out = {}
    for name, syms in sig.items():
        out[name] = syms
        if rem is not None:
            dep = [s for s in syms if s not in rem]
            out[name + "_CCDEP"] = dep
    return out

def cell_auc_for_sigs(csr_genes_cells, genes, variants, cell_idx=None, max_cells=None, verbose=True, already_csc=False):
    """csr_genes_cells: genes x cells sparse. genes: list matching rows.
    Returns DataFrame cells x variants with rank-based AUC scores.
    Operates per-cell on non-zero expression values (detected genes).
    If already_csc=True, the input is expected to be CSC."""
    G = csr_genes_cells.shape[0]
    gi = {g: i for i, g in enumerate(genes)}
    var_names = list(variants.keys())
    memb = np.zeros((G, len(var_names)), dtype=bool)
    for vi, name in enumerate(var_names):
        for s in variants[name]:
            ii = gi.get(s)
            if ii is not None:
                memb[ii, vi] = True
    if cell_idx is None:
        cell_idx = np.arange(csr_genes_cells.shape[1])
    if max_cells is not None and len(cell_idx) > max_cells:
        cell_idx = cell_idx[:max_cells]
    # work on CSC slice so that per-column .indices are gene row indices
    sl = csr_genes_cells if already_csc else csr_genes_cells.tocsc()
    if cell_idx is not None and (not already_csc or not np.array_equal(cell_idx, np.arange(csr_genes_cells.shape[1]))):
        sl = sl[:, cell_idx]
    V = len(var_names)
    out = np.full((sl.shape[1], V), np.nan)
    for cj in range(sl.shape[1]):
        col = sl[:, cj]
        d = col.data
        r = col.indices
        nd = d.size
        if nd == 0:
            continue
        uniq, inv, cnts = np.unique(d, return_inverse=True, return_counts=True)
        sums = np.zeros(uniq.size)
        # average-rank via within-unique mean after initial sorted rank
        order = np.argsort(d, kind="mergesort")
        tmp = np.empty(nd, dtype=np.float64)
        tmp[order] = np.arange(1, nd + 1, dtype=np.float64)
        np.add.at(sums, inv, tmp)
        ranks = sums[inv] / cnts[inv]
        gmem = memb[r, :]  # nd x V
        k = gmem.sum(axis=0).astype(np.float64)
        sr = np.where(k > 0, (ranks[:, None] * gmem).sum(axis=0) / np.where(k > 0, k, 1.0), np.nan)
        ok = (k > 0) & (nd > k)
        out[cj, ok] = (sr[ok] - (k[ok] + 1) / 2.0) / (nd - k[ok])
    return pd.DataFrame(out, index=cell_idx, columns=var_names)

def cell_altmean_for_sigs(csr_genes_cells, genes, variants, cell_idx=None, verbose=True):
    gi = {g: i for i, g in enumerate(genes)}
    var_idx = {}
    for name, syms in variants.items():
        var_idx[name] = np.array([gi[s] for s in syms if s in gi], dtype=np.int64)
    if cell_idx is None:
        cell_idx = np.arange(csr_genes_cells.shape[1])
    X = csr_genes_cells.tocsc()
    rows = {}
    for name, idx in var_idx.items():
        sub = X[idx, :][:, cell_idx]
        # log1p of values; mean over rows with any detection
        vals = sub.toarray()
        vals = np.log1p(vals)
        det = (vals > 0).sum(axis=0)
        s = vals.sum(axis=0)
        rows[name] = np.where(det > 0, s / np.maximum(det, 1), np.nan)
    return pd.DataFrame(rows, index=cell_idx)
