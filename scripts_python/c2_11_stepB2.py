# C2.11 Step B2: tests with correct row alignment (load cached pool AUC/info/mapped)
import os, numpy as np, pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

T = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(T), "output")

pool_auc = pd.read_csv(os.path.join(T, "c2_11_pool_auc.tsv"), sep="\t", index_col=0)
info = pd.read_csv(os.path.join(T, "c2_11_pool_info_mapped.tsv"), sep="\t")
core_sigs = [c for c in pool_auc.columns if not c.endswith("_CCDEP")]
STATE5 = ["MES", "AC", "OPC", "NPC", "Cycling"]

def patient_test(df, auc, sigcols, min_cell=5, min_pat=6, states=STATE5, key="mapped"):
    """df.index must equal auc row index (pool positions). Patient-level delta (state vs other-state), cells-only within cohort."""
    rows = []
    for c in sigcols:
        for s in states:
            vals, pats = [], []
            for p, g in df.groupby("patient", dropna=True):
                gs = g[g[key] == s]
                go = g[g[key] != s]
                if len(gs) >= min_cell and len(go) >= min_cell:
                    vals.append(float(auc.loc[gs.index, c].mean() - auc.loc[go.index, c].mean()))
                    pats.append(p)
            if len(vals) >= min_pat:
                v = np.array(vals)
                try:
                    w, pv = wilcoxon(v, zero_method="wilcox")
                except Exception:
                    pv = np.nan
                rows.append(dict(sig=c, state=s, n_pat=len(v), delta=float(v.mean()),
                                 frac_pos=float((v > 0).mean()), p=float(pv), patients=",".join(map(str, pats))))
    r = pd.DataFrame(rows)
    if len(r):
        r["fdr"] = multipletests(r["p"].fillna(1), method="fdr_bh")[1]
    return r

def pooled_test(df, auc, sigcols, label, min_cell=3, min_pat=6):
    rows = []
    for c in sigcols:
        vals, pats = [], []
        for p, g in df.groupby("patient", dropna=True):
            gn = g[g["core5"] == False]
            gc = g[g["core5"] == True]
            if len(gn) >= min_cell and len(gc) >= min_cell:
                vals.append(float(auc.loc[gn.index, c].mean() - auc.loc[gc.index, c].mean()))
                pats.append(p)
        if len(vals) >= min_pat:
            v = np.array(vals)
            try:
                w, pv = wilcoxon(v, zero_method="wilcox")
            except Exception:
                pv = np.nan
            rows.append(dict(label=label, sig=c, n_pat=len(v), delta=float(v.mean()),
                             frac_pos=float((v > 0).mean()), p=float(pv),
                             patients=",".join(map(str, pats))))
    r = pd.DataFrame(rows)
    if len(r):
        r["fdr"] = multipletests(r["p"].fillna(1), method="fdr_bh")[1]
    return r

la = info[info["is_LA"]].copy()
coreA = la[la["core5"]].copy()
la5 = la[la["mapped"].isin(STATE5)].copy()
print("LA core", len(coreA), "patients", coreA.patient.nunique())
print("LA mapped5", len(la5), "patients", la5.patient.nunique())

cc = [c for c in pool_auc.columns if c.endswith("_CCDEP")]
rA = patient_test(coreA, pool_auc, core_sigs)
rAcc = patient_test(coreA, pool_auc, cc)
rB = patient_test(la5, pool_auc, core_sigs)
rBcc = patient_test(la5, pool_auc, cc)

focus = ["C5_UP", "C5_STABLE_UP", "C5_DOWN", "C5_STABLE_DOWN", "M_HSF1_heat_shock"]
print("== A core base ==")
print(rA[rA.sig.isin(focus) & rA.state.isin(["MES", "NPC"])].to_string(index=False))
print("== A core CCDEP ==")
print(rAcc[rAcc.sig.str.replace("_CCDEP","").isin(focus) & rAcc.state.isin(["MES", "NPC"])].to_string(index=False))
print("== B mapped all ==")
print(rB[rB.sig.isin(focus) & rB.state.isin(["MES", "NPC"])].to_string(index=False))
print("== B CCDEP ==")
print(rBcc[rBcc.sig.str.replace("_CCDEP","").isin(focus) & rBcc.state.isin(["MES", "NPC"])].to_string(index=False))

# expanded pooled noncore vs core
rP = pooled_test(la, pool_auc, core_sigs, "LA_noncore_vs_core")
print("== LA pooled noncore vs core ==")
print(rP[rP.sig.isin(focus)].to_string(index=False))

# LB
lb = info[info["is_LB"]].copy()
lb5 = lb[lb["mapped"].isin(STATE5)].copy()
print("LB5 mapped counts", lb5["mapped"].value_counts().to_dict(), "patients", lb5.patient.nunique())
rLB = patient_test(lb5, pool_auc, core_sigs, min_cell=3, min_pat=2)
rLBcc = patient_test(lb5, pool_auc, cc, min_cell=3, min_pat=2)
print("== LB base =="); print(rLB.to_string(index=False))
print("== LB ccdep =="); print(rLBcc.to_string(index=False))

# persistence with full patient vectors (for patient-level tables)
def clean(r): return r.drop(columns=["patients"])
for name, r in [("A_core", rA), ("A_ccdep", rAcc), ("B_mapped", rB), ("B_ccdep", rBcc), ("LB", rLB), ("LB_ccdep", rLBcc), ("pooled", rP)]:
    clean(r).assign(cohort=name).to_csv(os.path.join(T, f"c2_11_test_{name}.tsv"), sep="\t", index=False)
print("STEP B2 DONE")
