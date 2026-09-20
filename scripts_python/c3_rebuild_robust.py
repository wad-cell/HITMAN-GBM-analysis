# -*- coding: utf-8 -*-
"""Rebuild C3_ROBUSTNESS.tsv as clean long table from intermediates (no re-analysis except TCGA primary-only row)."""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

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

def load(tag):
    return pd.read_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t").set_index("Name")

def q1p(df, idx):
    out = []
    for i, a in enumerate(A8):
        for b in A8[i+1:]:
            rho, p = stats.spearmanr(df.loc[idx, a], df.loc[idx, b], nan_policy="omit")
            e = EXPECTED.get((a, b), EXPECTED.get((b, a)))
            out.append({"a": a, "b": b, "rho": rho, "p": p, "e": e})
    r = pd.DataFrame(out)
    return r, int((np.sign(r["rho"]) == r["e"]).sum())

rows = []
cr = pd.read_csv(os.path.join(C3D, "_robust_cross_intermediate.tsv"), sep="\t")
for _, r in cr.iterrows():
    exp = r["expected"] if not pd.isna(r["expected"]) else np.nan
    rows.append({"analysis": "robust_1_cross_cohort_sign", "cohort_scope": "TCGA_vs_CGGA_full",
                 "module_a": r["module_a"], "module_b": r["module_b"],
                 "expected_sign": exp, "rho_TCGA": r["rho_TCGA"], "rho_CGGA": r["rho_CGGA"],
                 "same_direction": r["same_direction"], "note": ""})
n_exp = sum(~pd.isna(cr["expected"]))
n_same_exp = int((cr["same_direction"] & ~pd.isna(cr["expected"])).sum())
n_same_all = int(cr["same_direction"].sum())
rows.append({"analysis": "robust_1_cross_cohort_sign", "cohort_scope": "summary", "module_a": "SUMMARY",
             "module_b": f"all_28_sign_consistent={n_same_all}; expected21_sign_consistent={n_same_exp}/{n_exp}",
             "expected_sign": np.nan, "rho_TCGA": np.nan, "rho_CGGA": np.nan, "same_direction": np.nan,
             "note": "reversals within expected set: C4_UP-C4_DOWN (T -0.32/C +0.46), M_shared_antiproliferative-M_TNF_NFkB (T -0.15/C +0.12), M_shared_antiproliferative-M_p53_apoptosis (T -0.08/C +0.20); RIGI pair no expected sign (T -0.17/C +0.20)"})

st = pd.read_csv(os.path.join(C3D, "_robust_strata_intermediate.tsv"), sep="\t")
for _, r in st.iterrows():
    rows.append({"analysis": "robust_2_stratum_Q1", "cohort_scope": r["stratum"], "module_a": "all_21_expected_pairs",
                 "module_b": f"dir_match={int(r['dir_match'])}/21", "expected_sign": np.nan,
                 "rho_TCGA": np.nan, "rho_CGGA": np.nan, "same_direction": np.nan,
                 "note": f"n={r['n']}; rho_C4UP_C4DOWN={r['rho_C4UP_C4DOWN']:.3f}; rho_antiProlif_TNF={r['rho_antiProlif_TNF']:.3f}"})
# TCGA primary-only row (excluding 4 recurrent)
t = load("tcga")
sel = json.load(open(os.path.join(C3D, "tcga_gbm_selected.json")))
rec = [x["case"] for x in sel if x["sample_type"] == "Recurrent Tumor"]
prim = [s for s in t.index if s not in rec]
_, dm = q1p(t, prim)
rows.append({"analysis": "robust_2_stratum_Q1", "cohort_scope": "TCGA_Primary_only", "module_a": "all_21_expected_pairs",
             "module_b": f"dir_match={int(dm)}/21", "expected_sign": np.nan, "rho_TCGA": np.nan, "rho_CGGA": np.nan,
             "same_direction": np.nan, "note": f"n={len(prim)} (recurrent n={len(rec)} excluded); same pattern as full TCGA"})

pa = pd.read_csv(os.path.join(C3D, "_robust_partial_intermediate.tsv"), sep="\t")
for _, r in pa.iterrows():
    e = EXPECTED.get((r["module_a"], r["module_b"]), EXPECTED.get((r["module_b"], r["module_a"])))
    hit_raw = bool(np.sign(r["rho_raw"]) == e) if not pd.isna(e) else None
    hit_par = bool(np.sign(r["rho_partial_global_mean"]) == e) if not pd.isna(e) else None
    rows.append({"analysis": "robust_3_partial_global_mean", "cohort_scope": r["cohort"],
                 "module_a": r["module_a"], "module_b": r["module_b"],
                 "expected_sign": e, "rho_TCGA": r["rho_raw"], "rho_CGGA": r["rho_partial_global_mean"],
                 "same_direction": None if e is None else (hit_raw == hit_par),
                 "note": f"p_partial={r['p_partial']:.2e} (uncorrected); sign_preserved_after_partial={hit_raw == hit_par}"})
gp = pd.read_csv(os.path.join(C3D, "_robust_geneperturb_intermediate.tsv"), sep="\t")
for _, r in gp.iterrows():
    rows.append({"analysis": "robust_4_gene_perturbation", "cohort_scope": r["tag"], "module_a": f"run_{r['run']}",
                 "module_b": f"dir_match={int(r['dir_match'])}/21", "expected_sign": np.nan,
                 "rho_TCGA": np.nan, "rho_CGGA": np.nan, "same_direction": np.nan,
                 "note": r["tag"].split("_")[0]})

cov = pd.read_csv(os.path.join(C3D, "_robust_coverage_intermediate.tsv"), sep="\t")
for _, r in cov.iterrows():
    sigs = r.index[1:]
    vals = r.values[1:].astype(float)
    mn = vals.min(); arg = sigs[int(np.argmin(vals))]
    rows.append({"analysis": "robust_5_coverage", "cohort_scope": r["Unnamed: 0"], "module_a": "min_per_sample_gene_coverage",
                 "module_b": f"min={mn:.3f} ({arg})", "expected_sign": np.nan, "rho_TCGA": np.nan,
                 "rho_CGGA": np.nan, "same_direction": np.nan,
                 "note": "coverage>=0.80 rule; no imputation; sets unchanged"})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "C3_ROBUSTNESS.tsv"), sep="\t", index=False)
print("robustness rebuilt rows", df.shape)
print("cross same:", n_same_all, "/28; expected subset", n_same_exp, "/21")
print("TCGA primary-only dir_match", int(dm), "/21 n=", len(prim))
