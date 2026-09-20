# -*- coding: utf-8 -*-
"""C3 robustness part 2: gene-level perturbation (frozen method, gseapy ssgsea rank).
(a) remove top-10% highest-expressed genes per module (high-expression dominance check)
(b) 80% random gene dropout x5 seeds
Report Q1 expected-pair direction-match and key pair rhos per cohort.
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats
import gseapy as gp

C3D = os.path.dirname(os.path.abspath(__file__))
frozen = json.load(open(os.path.join(C3D, "frozen_sets.json"), encoding="utf-8"))
sets = frozen["sets"]
A8 = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
      "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]
EXPECTED = {
 ("M_cholesterol_SREBP", "M_TNF_NFkB"): 1, ("M_cholesterol_SREBP", "M_p53_apoptosis"): 1,
 ("M_TNF_NFkB", "M_p53_apoptosis"): 1, ("C4_UP", "M_cholesterol_SREBP"): 1,
 ("C4_UP", "M_TNF_NFkB"): 1, ("C4_UP", "M_p53_apoptosis"): 1,
 ("M_shared_antiproliferative", "M_shared_ER_UPR"): 1,
 ("C4_DOWN", "M_shared_antiproliferative"): 1, ("C4_DOWN", "M_shared_ER_UPR"): 1,
 ("C4_UP", "C4_DOWN"): -1, ("M_shared_antiproliferative", "M_cholesterol_SREBP"): -1,
 ("M_shared_antiproliferative", "M_TNF_NFkB"): -1, ("M_shared_antiproliferative", "M_p53_apoptosis"): -1,
 ("M_shared_ER_UPR", "M_cholesterol_SREBP"): -1, ("M_shared_ER_UPR", "M_TNF_NFkB"): -1,
 ("M_shared_ER_UPR", "M_p53_apoptosis"): -1, ("C4_DOWN", "M_cholesterol_SREBP"): -1,
 ("C4_DOWN", "M_TNF_NFkB"): -1, ("C4_DOWN", "M_p53_apoptosis"): -1,
 ("C4_UP", "M_shared_antiproliferative"): -1, ("C4_UP", "M_shared_ER_UPR"): -1,
}
PAIRS = list(EXPECTED.keys())

def load(tag):
    z = np.load(os.path.join(C3D, f"{tag}_matrix.npz"), allow_pickle=True)
    mat, samples = z["mat"], list(z["samples"])
    genes = json.load(open(os.path.join(C3D, "tcga_genes.json")))["gene"] if tag == "tcga" else list(z["genes"])
    X = np.log2(mat.astype(np.float64) + 1.0)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    df = pd.DataFrame(X.T, index=genes, columns=samples)
    df = df[(df.sum(axis=1) > 0)]
    return df

def score(df, gs):
    res = gp.ssgsea(data=df, gene_sets=gs, outdir=None, sample_norm_method="rank",
                    min_size=2, max_size=20000, permutation_num=0, verbose=False)
    wide = res.res2d.pivot_table(index="Name", columns="Term", values="ES")
    return wide

def pair_stats(wide):
    rows = []
    for a, b in PAIRS:
        rho, p = stats.spearmanr(wide[a], wide[b], nan_policy="omit")
        rows.append((a, b, rho, p))
    return rows

def summarize(rows):
    match = 0
    out = {}
    for a, b, rho, p in rows:
        exp = EXPECTED[(a, b)]
        hit = bool(np.sign(rho) == exp)
        match += hit
        out[f"{a}|{b}"] = rho
    return match, out

def universe_present(df):
    return set(df.index)

def run_all(tag):
    df = load(tag)
    uni = universe_present(df)
    base_gs = {k: [g for g in sets[k] if g in uni] for k in A8}
    record = []
    # baseline
    wide0 = score(df, base_gs)
    m0, rhos0 = summarize(pair_stats(wide0))
    record.append({"tag": tag, "run": "baseline", "dir_match": m0, **{f"rho_{k.replace(chr(124),'__')}": v for k, v in rhos0.items()}})
    # (a) drop top-10% expressed genes per module
    trim_gs = {}
    for k in A8:
        gs = base_gs[k]
        if len(gs) < 20:
            trim_gs[k] = gs
            continue
        mean_expr = df.loc[gs].mean(axis=1)
        keep = mean_expr.index[mean_expr.values < np.percentile(mean_expr.values, 90)]
        trim_gs[k] = list(keep)
    wide_t = score(df, trim_gs)
    mt, rhost = summarize(pair_stats(wide_t))
    record.append({"tag": tag, "run": "drop_top10expr", "dir_match": mt, **{f"rho_{k.replace(chr(124),'__')}": v for k, v in rhost.items()}})
    # (b) 80% gene dropout x5
    rng = np.random.default_rng(20260909)
    for i in range(5):
        seed = rng.integers(0, 2**31)
        rg = np.random.default_rng(seed)
        gs_i = {}
        for k in A8:
            gs = np.asarray(base_gs[k], dtype=object)
            n = max(2, int(np.ceil(len(gs) * 0.8)))
            gs_i[k] = list(rg.choice(gs, size=n, replace=False))
        wi = score(df, gs_i)
        mi, rhosi = summarize(pair_stats(wi))
        record.append({"tag": tag, "run": f"dropout80_{i+1}", "dir_match": mi, **{f"rho_{k.replace(chr(124),'__')}": v for k, v in rhosi.items()}})
    return pd.DataFrame(record)

recs = []
for tag in ["tcga", "cgga"]:
    recs.append(run_all(tag))
out = pd.concat(recs, ignore_index=True)
out.to_csv(os.path.join(C3D, "_robust_geneperturb_intermediate.tsv"), sep="\t", index=False)
print(out[["tag", "run", "dir_match"]].to_string(index=False))
pd.set_option("display.width", 250)
print(out.to_string(index=False))
print("done")
