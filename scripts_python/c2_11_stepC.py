# C2.11 Step C: patient-level effects + LOPO + mapping diagnostic + LB state tests placeholder
import os, numpy as np, pandas as pd
from scipy.stats import wilcoxon

T = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(T), "output")
pool_auc = pd.read_csv(os.path.join(T, "c2_11_pool_auc.tsv"), sep="\t", index_col=0)
info = pd.read_csv(os.path.join(T, "c2_11_pool_info_mapped.tsv"), sep="\t")
STATE5 = ["MES", "AC", "OPC", "NPC", "Cycling"]
la = info[info["is_LA"]].copy()
coreA = la[la["core5"]].copy()
la5 = la[la["mapped"].isin(STATE5)].copy()
focus = ["C5_UP", "C5_STABLE_UP", "C5_DOWN", "C5_STABLE_DOWN", "M_HSF1_heat_shock"]
focus_cc = [c + "_CCDEP" for c in focus]

def per_patient(df, sigcols, states):
    rows = []
    for c in sigcols:
        for s in states:
            for p, g in df.groupby("patient", dropna=True):
                gs = g[g["mapped"] == s]
                go = g[g["mapped"] != s]
                if len(gs) >= 5 and len(go) >= 5:
                    rows.append(dict(sig=c, state=s, patient=p, n_state=len(gs), n_other=len(go),
                                     delta=float(pool_auc.loc[gs.index, c].mean() - pool_auc.loc[go.index, c].mean())))
    return pd.DataFrame(rows)

def lopo(df, sigcols, states, out_key):
    rows = []
    for c in sigcols:
        for s in states:
            for rem in sorted(df["patient"].dropna().unique()):
                sub = df[df["patient"] != rem]
                vals = []
                for p, g in sub.groupby("patient", dropna=True):
                    gs = g[g["mapped"] == s]; go = g[g["mapped"] != s]
                    if len(gs) >= 5 and len(go) >= 5:
                        vals.append(float(pool_auc.loc[gs.index, c].mean() - pool_auc.loc[go.index, c].mean()))
                if len(vals) >= 6:
                    v = np.array(vals)
                    try: pv = wilcoxon(v, zero_method="wilcox")[1]
                    except Exception: pv = np.nan
                    rows.append(dict(cohort=out_key, sig=c, state=s, removed=rem, n_pat=len(v),
                                     delta=float(v.mean()), frac_pos=float((v > 0).mean()), p=pv))
    return pd.DataFrame(rows)

for name, df in [("A_core", coreA), ("B_mapped_all", la5)]:
    pp = per_patient(df, focus + focus_cc, ["MES", "NPC", "Cycling"])
    pp.to_csv(os.path.join(OUT, f"C2_11_PATIENT_EFFECTS_{name}.tsv"), sep="\t", index=False)
    lo = lopo(df, focus, ["MES", "NPC"], name)
    lo.to_csv(os.path.join(OUT, f"C2_11_LOPO_{name}.tsv"), sep="\t", index=False)
    print(name, "patient rows", len(pp), "lopo rows", len(lo))

# mapping diagnostic: original expanded category -> mapped state
mm = la[~la["core5"]][["mp_top", "mapped", "patient"]].copy()
mm["mp_top"] = mm["mp_top"].fillna("NA")
tab = pd.crosstab(mm["mp_top"], mm["mapped"])
tab.insert(0, "cells", mm.groupby("mp_top").size())
tab.insert(1, "patients", mm.groupby("mp_top")["patient"].nunique())
tab.to_csv(os.path.join(OUT, "C2_11_EXPANDED_MAPPING_DIAG.tsv"), sep="\t")
print("mapping diag saved", tab.shape)

# LB state tests file (not testable rows documented)
lb = info[info["is_LB"]].copy()
lb5 = lb[lb["mapped"].isin(STATE5)].copy()
rows = []
for s in STATE5:
    sub = lb5[lb5["mapped"] == s]
    rows.append(dict(state=s, cells=int(len(sub)), patients=int(sub["patient"].nunique()),
                     reason="per-patient state test not applicable (cells<3 or patients<2)",
                     n_pat_min_required=2, n_cell_min_required=3))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "C2_LAYERB_STATE_TESTS.tsv"), sep="\t", index=False)
print("LB state tests rows:")
print(pd.DataFrame(rows).to_string(index=False))
print("STEP C DONE")
