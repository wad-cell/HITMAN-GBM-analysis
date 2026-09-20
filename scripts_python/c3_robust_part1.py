# -*- coding: utf-8 -*-
"""C3 robustness part 1 (no gene re-scoring):
- cross-cohort consistency summary
- partial Spearman controlling per-sample global expression (median log2)
- CGGA IDH (Wildtype vs Mutant) and prs_type (Primary vs Recurrent) stratified Q1 + axis
- TCGA main results after excluding 4 recurrent samples
- signature per-sample coverage audit (min per sample)
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C3D = os.path.join(BASE, "temp", "c3_data")
OUT = os.path.join(BASE, "output")

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

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    adj = np.minimum.accumulate(ranked * n / np.arange(1, n + 1))
    out = np.empty(n); out[order] = np.minimum(adj, 1.0)
    return out

def load_scores(tag):
    df = pd.read_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t").set_index("Name")
    return df

def q1_pairs(df, sample_mask=None):
    if sample_mask is not None:
        df = df.loc[sample_mask]
    rows = []
    for i, a in enumerate(A8):
        for b in A8[i + 1:]:
            rho, p = stats.spearmanr(df[a], df[b], nan_policy="omit")
            rows.append({"module_a": a, "module_b": b, "rho": rho, "p": p,
                         "expected": EXPECTED.get((a, b), EXPECTED.get((b, a), np.nan))})
    r = pd.DataFrame(rows)
    return r

def global_expr(tag):
    z = np.load(os.path.join(C3D, f"{tag}_matrix.npz"), allow_pickle=True)
    mat = z["mat"]; samples = list(z["samples"])
    genes = json.load(open(os.path.join(C3D, "tcga_genes.json")))["gene"] if tag == "tcga" else list(z["genes"])
    X = np.log2(mat.astype(np.float64) + 1.0)
    # drop all-zero rows then mean per sample (proxy of global transcript level)
    med = np.nanmedian(X, axis=1)
    mean = np.nanmean(X, axis=1)
    return pd.Series(mean, index=samples, name="global_mean"), pd.Series(med, index=samples, name="global_median")

def partial_spearman(x, y, c):
    def rankr(s):
        return stats.rankdata(s)
    rx, ry, rc = rankr(x), rankr(y), rankr(c)
    # residuals
    def resid(r, rc):
        A = np.vstack([rc, np.ones(len(rc))]).T
        coef, *_ = np.linalg.lstsq(A, r, rcond=None)
        return r - A @ coef
    ex, ey = resid(rx, rc), resid(ry, rc)
    return stats.spearmanr(ex, ey, nan_policy="omit")

def run():
    t = load_scores("tcga"); g = load_scores("cgga")
    cc = pd.read_csv(os.path.join(C3D, "cgga_clinical_iv.tsv"), sep="\t").set_index("id")
    sel = json.load(open(os.path.join(C3D, "tcga_gbm_selected.json")))
    rec_cases = [x["case"] for x in sel if x["sample_type"] == "Recurrent Tumor"]
    print("TCGA recurrent cases:", rec_cases)

    gm_t, gmed_t = global_expr("tcga")
    gm_g, gmed_g = global_expr("cgga")
    print("global mean available:", gm_t.index.isin(t.index).sum(), gm_g.index.isin(g.index).sum())

    rows_cross = []
    for a in A8:
        for b in A8[A8.index(a) + 1:]:
            rt = q1_pairs(t, None)
            rt = rt[(rt.module_a == a) & (rt.module_b == b)].iloc[0]
            rg = q1_pairs(g, None)
            rg = rg[(rg.module_a == a) & (rg.module_b == b)].iloc[0]
            rows_cross.append({"module_a": a, "module_b": b, "expected": rt.expected,
                               "rho_TCGA": rt.rho, "rho_CGGA": rg.rho,
                               "same_direction": np.sign(rt.rho) == np.sign(rg.rho)})
    cross = pd.DataFrame(rows_cross)
    cross["rho_TCGA_sign"] = np.sign(cross["rho_TCGA"])
    cross["rho_CGGA_sign"] = np.sign(cross["rho_CGGA"])
    print("cross-cohort sign-consistent:", cross["same_direction"].sum(), "/", len(cross))

    # partial spearman controlling global mean expression
    partial_rows = []
    for tag, df, gmean in [("TCGA", t, gm_t), ("CGGA", g, gm_g)]:
        idx = df.index.intersection(gmean.index)
        cc_vals = gmean.loc[idx]
        for a in A8:
            for b in A8[A8.index(a) + 1:]:
                rho_raw, _ = stats.spearmanr(df.loc[idx, a], df.loc[idx, b])
                rho_p, p_p = partial_spearman(df.loc[idx, a].values, df.loc[idx, b].values, cc_vals.values)
                partial_rows.append({"cohort": tag, "module_a": a, "module_b": b,
                                     "rho_raw": rho_raw, "rho_partial_global_mean": rho_p, "p_partial": p_p})
    part = pd.DataFrame(partial_rows)
    # how many expected-sign pairs keep sign after partial
    for coh in ["TCGA", "CGGA"]:
        pp = part[part.cohort == coh].copy()
        exp = pp.apply(lambda r: EXPECTED.get((r.module_a, r.module_b), EXPECTED.get((r.module_b, r.module_a), np.nan)), axis=1)
        pp["expected"] = exp.values
        pp["raw_match"] = np.sign(pp["rho_raw"]) == exp.values
        pp["partial_match"] = np.sign(pp["rho_partial_global_mean"]) == exp.values
        pp = pp[exp.notna()]
        print(coh, "raw dir-match:", pp["raw_match"].sum(), "/", len(pp),
              " partial dir-match:", pp["partial_match"].sum(), "/", len(pp))
    part.to_csv(os.path.join(C3D, "_robust_partial_intermediate.tsv"), sep="\t", index=False)

    # CGGA IDH strata & prs_type strata: Q1 architecture summary + key MES-axis pairs
    strata = {}
    strata["CGGA_IDH_WT"] = cc.index[cc["idh"] == "Wildtype"]
    strata["CGGA_IDH_Mut"] = cc.index[cc["idh"] == "Mutant"]
    strata["CGGA_Primary"] = cc.index[cc["prs_type"] == "Primary"]
    strata["CGGA_Recurrent"] = cc.index[cc["prs_type"] == "Recurrent"]
    strat_rows = []
    for name, idx in strata.items():
        m = g.index.intersection(idx)
        q = q1_pairs(g, m)
        qe = q[q.expected.notna()]
        strat_rows.append({"stratum": name, "n": int(len(m)),
                           "dir_match": int((np.sign(qe.rho) == qe.expected).sum()),
                           "total_pairs": int(len(qe)),
                           "rho_C4UP_C4DOWN": float(q[(q.module_a == "C4_UP") & (q.module_b == "C4_DOWN")]["rho"].iloc[0]),
                           "rho_antiProlif_TNF": float(q[(q.module_a == "M_shared_antiproliferative") & (q.module_b == "M_TNF_NFkB")]["rho"].iloc[0])})
    strat = pd.DataFrame(strat_rows)
    print(strat.to_string())

    # TCGA remove recurrent
    t_prim = t.loc[~t.index.isin(rec_cases)]
    print("TCGA primary-only n:", len(t_prim))
    qp = q1_pairs(t_prim)
    qe = qp[qp.expected.notna()]
    print("TCGA primary-only dir-match:", int((np.sign(qe.rho) == qe.expected).sum()), "/", len(qe))

    # C4_UP-axis association per stratum (exploratory context for C5/state)
    axis_rows = []
    for tag, df, gmean in [("TCGA", t, gm_t), ("CGGA", g, gm_g)]:
        for a in A8 + ["C5_UP", "C5_STABLE_UP", "C5_DOWN", "C5_STABLE_DOWN"]:
            rho, p = stats.spearmanr(df["C4_UP"], df[a], nan_policy="omit")
            axis_rows.append({"cohort": tag, "score": a, "C4UP_rho": rho, "p": p})
    axis = pd.DataFrame(axis_rows)
    axis["FDR_BH"] = bh(axis["p"].values)
    axis.to_csv(os.path.join(C3D, "_robust_axis_intermediate.tsv"), sep="\t", index=False)

    # coverage audit: min per-sample coverage is not stored; we stored per-set coverage json.
    cov = {}
    for tag in ["tcga", "cgga"]:
        cov[tag] = json.load(open(os.path.join(C3D, f"{tag}_set_coverage.json")))
    covdf = pd.DataFrame(cov).T
    covdf.to_csv(os.path.join(C3D, "_robust_coverage_intermediate.tsv"), sep="\t")
    print("coverage min per set:", covdf.min().min(), "below 0.80? ", (covdf < 0.8).any().any())

    # save intermediates
    strat.to_csv(os.path.join(C3D, "_robust_strata_intermediate.tsv"), sep="\t", index=False)
    cross.to_csv(os.path.join(C3D, "_robust_cross_intermediate.tsv"), sep="\t", index=False)
    print("done")

if __name__ == "__main__":
    run()
