# C2.11 Step B: mapping, sensitivity tests, patient-level diagnostics
import os, numpy as np, pandas as pd, scipy.sparse as sp, sys, time
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

T = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(T), "output")

pool_auc = pd.read_csv(os.path.join(T, "c2_11_pool_auc.tsv"), sep="\t", index_col=0)
info = pd.read_csv(os.path.join(T, "c2_11_pool_info.tsv"), sep="\t")
pos = np.load(os.path.join(T, "c2_11_pool_pos.npy"))
lb_pos = np.load(os.path.join(T, "c2_11_pool_lb.npy"))
is_lb = np.isin(pos, lb_pos)
core_sigs = [c for c in pool_auc.columns if not c.endswith("_CCDEP")]
STATE5 = ["MES", "AC", "OPC", "NPC", "Cycling"]
OUT_STATE5 = ["MES", "AC", "OPC", "NPC"]
info["is_LA"] = (info["pool"] == "LA_tissue") | (info["pool"] == "both")
info["is_LB"] = info["pool"].isin(["LB", "both"])
info["core5"] = info["state6"].isin(STATE5)

# ---------- 1) log1p normalized pool matrix (for mapping only) ----------
t0 = time.time()
Ml = sp.load_npz(os.path.join(T, "c2_11_pool_csc.npz")).tocsc()
cs = np.asarray(Ml.sum(axis=0)).ravel()
cs[cs == 0] = 1.0
scale = 1e4 / cs
n = Ml.shape[1]
counts = np.diff(Ml.indptr).astype(np.int64)
col_of = np.repeat(np.arange(n, dtype=np.int64), counts)
Xlog = Ml.copy()
Xlog.data = Xlog.data * scale[col_of]
Xlog.data = np.log1p(Xlog.data)
print("log matrix", Xlog.shape, f"{t0:.0f}->{time.time()-t0:.0f}s")

# centroids from LA core four non-cycling states
core4 = info["is_LA"] & info["state6"].isin(OUT_STATE4 := ["MES", "AC", "OPC", "NPC"])
sub = Xlog[:, np.where(core4)[0]]
csum = np.asarray(sub.sum(axis=1)).ravel() / max(1, int(core4.sum()))
cent_all = {s: csum for s in ["NA"]}  # placeholder replaced below by per-state means
cent = {}
for s in OUT_STATE4:
    mask = (info["state6"] == s) & info["is_LA"]
    sm = Xlog[:, np.where(mask)[0]]
    cent[s] = (np.asarray(sm.sum(axis=1)).ravel() / max(1, int(mask.sum()))).astype(np.float64)
Cmat = np.stack([cent[s] for s in OUT_STATE4], axis=1)   # genes x 4
Cm = Cmat - Cmat.mean(axis=0)
Cn = np.linalg.norm(Cm, axis=0)
Cn[Cn == 0] = 1.0
Cm /= Cn

def map_chunk(cols, batch=600):
    """corr-based mapping: returns assigned state + max corr + margin"""
    ass = np.full(len(cols), "", dtype=object)
    mcorr = np.full(len(cols), np.nan)
    margin = np.full(len(cols), np.nan)
    for b in range(0, len(cols), batch):
        ci = cols[b:b+batch]
        Xb = Xlog[:, ci].T.toarray().astype(np.float64)
        Xb = Xb - Xb.mean(axis=1, keepdims=True)
        nrm = np.linalg.norm(Xb, axis=1)
        nrm[nrm == 0] = 1.0
        Xb /= nrm[:, None]
        corr = Xb @ Cm                       # cells x 4
        order = np.argsort(-corr, axis=1)
        top = corr[np.arange(len(ci)), order[:, 0]]
        sec = corr[np.arange(len(ci)), order[:, 1]]
        mc = corr[np.arange(len(ci)), order[:, 0]]
        for k in range(len(ci)):
            if mc[k] <= 0.05 or (mc[k] - sec[k]) <= 0.02:
                ass[b+k] = "OTHER"
            else:
                ass[b+k] = OUT_STATE4[order[k, 0]]
            mcorr[b+k] = mc[k]
            margin[b+k] = mc[k] - sec[k]
    return ass, mcorr, margin

# ---------- 2) mapping for LA non-core and LB ----------
info["phase_cycle"] = info["cell_cycle_phase"].fillna("").astype(str).isin(["G1/S", "G2/M"])
def assign_mapped(df):
    st = df["state6"].astype(str).copy()
    mp = df["mp_top"].fillna("").astype(str)
    out = pd.Series(np.where(df["core5"], st, "unset"), index=df.index)
    # residual non-core cells
    cand = df.index[(df["core5"] == False)].tolist()
    if cand:
        ass, mcorr, marg = map_chunk(np.array(cand))
        aS = pd.Series(ass, index=cand)
        out.loc[cand] = np.where(aS.values == "OTHER", "OTHER", aS.values)
        # cycle-priority: phase cycle or mp_top Cell Cycle => Cycling
        cyc = [i for i in cand if (df.loc[i, "phase_cycle"]) or str(mp.loc[i]).startswith("Cell Cycle")]
        for i in cyc:
            out.loc[i] = "Cycling"
    out[df.index[df["is_LB"] & ~df["core5"] & df["phase_cycle"].values]] = "Cycling"
    return out, None, None

info["mapped"] = "unset"
info["mapped"], _, _ = assign_mapped(info)
print("mapped counts:\n", info["mapped"].value_counts(dropna=False).to_dict())
info.to_csv(os.path.join(T, "c2_11_pool_info_mapped.tsv"), sep="\t", index=False)

# ---------- 3) helpers ----------
def patient_test(df, auc, sigcols, min_cell=5, min_pat=6, states=STATE5):
    """per-patient within-state vs other-states deltas; returns DataFrame rows"""
    rows = []
    dd = df.reset_index(drop=True)
    aa = auc.iloc[dd.index]
    for c in sigcols:
        for s in states:
            vals = []
            for p, g in dd.groupby("patient", dropna=True):
                gs = g[g["mapped"] == s]
                go = g[g["mapped"] != s]
                if len(gs) >= min_cell and len(go) >= min_cell:
                    vals.append(float(aa.loc[gs.index, c].mean() - aa.loc[go.index, c].mean()))
            if len(vals) >= min_pat:
                v = np.array(vals)
                try:
                    w, pv = wilcoxon(v, zero_method="wilcox")
                except Exception:
                    pv = np.nan
                rows.append(dict(sig=c, state=s, n_pat=len(v), delta=float(v.mean()),
                                 frac_pos=float((v > 0).mean()), p=pv, deltas=v))
    r = pd.DataFrame([{k: v for k, v in x.items() if k != "deltas"} for x in rows])
    if len(r):
        r["fdr"] = multipletests(r["p"].fillna(1), method="fdr_bh")[1]
    return r

# ---------- 4) cohort tests ----------
la = info[info["is_LA"]].copy()
la5 = la[la["mapped"].isin(STATE5)].copy()
coreA = la[la["core5"]].copy()

rA = patient_test(coreA, pool_auc, core_sigs)
print("== A core LA ==")
print(rA[rA.sig.isin(["C5_UP", "C5_STABLE_UP", "C5_DOWN", "M_HSF1_heat_shock"])].to_string(index=False))
rB = patient_test(la5, pool_auc, core_sigs)
print("== B all mapped LA ==")
print(rB[rB.sig.isin(["C5_UP", "C5_STABLE_UP", "C5_DOWN", "M_HSF1_heat_shock"])].to_string(index=False))

# CCDEP
cc_sigs = [c for c in pool_auc.columns if c.endswith("_CCDEP")]
rA_cc = patient_test(coreA, pool_auc, cc_sigs)
rB_cc = patient_test(la5, pool_auc, cc_sigs)
print("== A CCDEP MES ==")
print(rA_cc[rA_cc.sig.str.startswith(("C5_UP", "C5_STABLE_UP", "M_HSF1")) & (rA_cc.state=="MES")].to_string(index=False))
print("== B CCDEP MES ==")
print(rB_cc[rB_cc.sig.str.startswith(("C5_UP", "C5_STABLE_UP", "M_HSF1")) & (rB_cc.state=="MES")].to_string(index=False))

# LB tests
lb = info[info["is_LB"]].copy()
lb5 = lb[lb["mapped"].isin(STATE5)].copy()
print("LB5 table:", lb5["mapped"].value_counts().to_dict(), "patients:", lb5["patient"].nunique())
rLB = patient_test(lb5, pool_auc, core_sigs, min_cell=3, min_pat=2)
print("== LB ==")
print(rLB.head(40).to_string(index=False))

# pooled non-core vs core within patient (expanded-state sensitivity D)
def pooled_test(df, auc, sigcols, label, min_cell=3):
    rows = []
    dd = df.reset_index(drop=True)
    aa = auc.iloc[dd.index]
    for c in sigcols:
        vals = []
        for p, g in dd.groupby("patient", dropna=True):
            gn = g[g["core5"] == False]
            gc = g[g["core5"] == True]
            if len(gn) >= min_cell and len(gc) >= min_cell:
                vals.append(float(aa.loc[gn.index, c].mean() - aa.loc[gc.index, c].mean()))
        if len(vals) >= 6:
            v = np.array(vals)
            try:
                w, pv = wilcoxon(v, zero_method="wilcox")
            except Exception:
                pv = np.nan
            rows.append(dict(label=label, sig=c, n_pat=len(v), delta=float(v.mean()),
                             frac_pos=float((v > 0).mean()), p=pv))
    r = pd.DataFrame(rows)
    if len(r):
        r["fdr"] = multipletests(r["p"].fillna(1), method="fdr_bh")[1]
    return r

rP = pooled_test(la, pool_auc, core_sigs, "noncore_vs_core")
print("== pooled noncore vs core (per patient) ==")
print(rP[rP.sig.isin(["C5_UP", "C5_STABLE_UP", "M_HSF1_heat_shock"])].to_string(index=False))

# category summary (expanded-state) descriptive
la2 = la.copy()
cat = la2[la2["core5"] == False]["mp_top"].fillna("OTHER").value_counts()
cat_df = []
for lbl, n in cat.items():
    sub = la2[(la2["mp_top"].fillna("OTHER") == lbl)]
    cat_df.append(dict(expanded_label=lbl, cells=int(n), patients=sub["patient"].nunique()))
cat_out = pd.DataFrame(cat_df).sort_values("cells", ascending=False)
print("categories:", len(cat_out))
cat_out.to_csv(os.path.join(OUT, "C2_11_EXPANDED_CATEGORY_COUNTS.tsv"), sep="\t", index=False)

# ---------- 5) persist ----------
def clean(r):
    return r.drop(columns=["deltas"]) if "deltas" in r else r

for name, r in [("A_core", rA), ("B_mapped", rB), ("A_ccdep", rA_cc), ("B_ccdep", rB_cc), ("LB", rLB), ("pooled", rP)]:
    clean(r).assign(cohort=name).to_csv(os.path.join(T, f"c2_11_test_{name}.tsv"), sep="\t", index=False)
print("ALL STEP B DONE", time.time())
