# -*- coding: utf-8 -*-
"""Rebuild C3_PRIMARY_RESULTS.tsv with unified schema: Q1 architecture + Q2 MES context blocks."""
import os
import numpy as np
import pandas as pd

CONV = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C3D = os.path.join(CONV, "temp", "c3_data")
OUT = os.path.join(CONV, "output")

A8 = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
      "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]
EXPECTED = {("M_cholesterol_SREBP", "M_TNF_NFkB"): 1, ("M_cholesterol_SREBP", "M_p53_apoptosis"): 1,
 ("M_TNF_NFkB", "M_p53_apoptosis"): 1, ("C4_UP", "M_cholesterol_SREBP"): 1,
 ("C4_UP", "M_TNF_NFkB"): 1, ("C4_UP", "M_p53_apoptosis"): 1,
 ("M_shared_antiproliferative", "M_shared_ER_UPR"): 1, ("C4_DOWN", "M_shared_antiproliferative"): 1,
 ("C4_DOWN", "M_shared_ER_UPR"): 1, ("C4_UP", "C4_DOWN"): -1,
 ("M_shared_antiproliferative", "M_cholesterol_SREBP"): -1, ("M_shared_antiproliferative", "M_TNF_NFkB"): -1,
 ("M_shared_antiproliferative", "M_p53_apoptosis"): -1, ("M_shared_ER_UPR", "M_cholesterol_SREBP"): -1,
 ("M_shared_ER_UPR", "M_TNF_NFkB"): -1, ("M_shared_ER_UPR", "M_p53_apoptosis"): -1,
 ("C4_DOWN", "M_cholesterol_SREBP"): -1, ("C4_DOWN", "M_TNF_NFkB"): -1,
 ("C4_DOWN", "M_p53_apoptosis"): -1, ("C4_UP", "M_shared_antiproliferative"): -1,
 ("C4_UP", "M_shared_ER_UPR"): -1}

def load_scores(tag):
    return pd.read_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t").set_index("Name")

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); rp = p[o]
    raw = rp * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(raw[::-1])[::-1]   # suffix-min = standard BH monotone
    out = np.empty(n); out[o] = np.minimum(adj, 1.0); return out

def fisher_ci(rho, n):
    z = np.arctanh(np.clip(rho, -0.9999, 0.9999)); se = 1 / np.sqrt(n - 3); q = 1.96
    return np.tanh(z - q * se), np.tanh(z + q * se)

from scipy import stats
# coverage json is flat {set: fraction}; use plain json load
import json as _json
covj = {"TCGA": _json.load(open(os.path.join(C3D, "tcga_set_coverage.json"))),
        "CGGA": _json.load(open(os.path.join(C3D, "cgga_set_coverage.json")))}

t, g = load_scores("tcga"), load_scores("cgga")
rows = []
for coh, df, nn in [("TCGA", t, len(t)), ("CGGA", g, len(g))]:
    for i, a in enumerate(A8):
        for b in A8[i + 1:]:
            rho, p = stats.spearmanr(df[a], df[b], nan_policy="omit")
            lo, hi = fisher_ci(rho, nn)
            exp = EXPECTED.get((a, b), EXPECTED.get((b, a)))
            rows.append({"block": "Q1_architecture", "cohort": coh, "module_a": a, "module_b": b,
                         "effect_rho": rho, "CI95_low": lo, "CI95_high": hi, "P": p, "FDR_BH": np.nan,
                         "expected_sign": np.nan if exp is None else float(exp),
                         "direction_match": np.nan if exp is None else bool(np.sign(rho) == exp),
                         "min_pair_set_coverage": round(min(covj[coh][a], covj[coh][b]), 3),
                         "n": nn, "grade": "A", "note": ""})
q1 = pd.DataFrame(rows)
for coh in ["TCGA", "CGGA"]:
    m = q1["cohort"] == coh
    q1.loc[m, "FDR_BH"] = bh(q1.loc[m, "P"].values)

q2 = pd.read_csv(os.path.join(C3D, "_q2_mes_intermediate.tsv"), sep="\t")
q2fdr = bh(q2["mw_p"].values)
q2b = []
for (_, r), f in zip(q2.iterrows(), q2fdr):
    q2b.append({"block": "Q2_MES_vs_nonMES_context", "cohort": "TCGA_official_Verhaak_n44",
                "module_a": r["score"], "module_b": "MES(n=15)_vs_nonMES(n=29)",
                "effect_rho": np.nan, "CI95_low": r["CI95_low"], "CI95_high": r["CI95_high"],
                "P": r["mw_p"], "FDR_BH": f,
                "expected_sign": np.nan, "direction_match": np.nan, "min_pair_set_coverage": np.nan,
                "n": 44, "grade": "A",
                "note": f"effect=median ssGSEA score difference (MES minus non-MES)={r['median_diff_1minus0']:.3f}"})
q2f = pd.DataFrame(q2b)
final = pd.concat([q1, q2f], ignore_index=True, sort=False)
final.to_csv(os.path.join(OUT, "C3_PRIMARY_RESULTS.tsv"), sep="\t", index=False)
print("PRIMARY rebuilt rows", final.shape)
print(q2f[["module_a", "effect_rho", "P", "FDR_BH", "note"]].to_string(index=False))
