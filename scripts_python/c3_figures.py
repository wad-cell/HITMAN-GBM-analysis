# -*- coding: utf-8 -*-
"""C3 Figure-6 candidate & QC plots (honest figures, English labels, no beautification)."""
import os, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

CONV = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C3D = os.path.join(CONV, "temp", "c3_data")
OUT = os.path.join(CONV, "output")
FIG = os.path.join(OUT, "C3_FIG6_CANDIDATE")
os.makedirs(FIG, exist_ok=True)

A8 = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
      "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]
LAB = {k: k for k in A8}
LAB.update({"M_shared_ER_UPR": "shared_ER-UPR", "M_cholesterol_SREBP": "SREBP/cholesterol",
            "M_p53_apoptosis": "p53/apoptosis", "M_TNF_NFkB": "TNF/NF-kB",
            "M_RIGI_typeI_IFN": "RIG-I/IFN", "M_shared_antiproliferative": "shared_antiProlif"})
ORDER = ["C4_UP", "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN",
         "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR"]
EXPECTED = {("C4_UP", "M_cholesterol_SREBP"): 1, ("C4_UP", "M_TNF_NFkB"): 1, ("C4_UP", "M_p53_apoptosis"): 1,
 ("M_cholesterol_SREBP", "M_TNF_NFkB"): 1, ("M_cholesterol_SREBP", "M_p53_apoptosis"): 1,
 ("M_TNF_NFkB", "M_p53_apoptosis"): 1, ("M_shared_antiproliferative", "M_shared_ER_UPR"): 1,
 ("C4_DOWN", "M_shared_antiproliferative"): 1, ("C4_DOWN", "M_shared_ER_UPR"): 1,
 ("C4_UP", "C4_DOWN"): -1, ("C4_UP", "M_shared_antiproliferative"): -1, ("C4_UP", "M_shared_ER_UPR"): -1,
 ("C4_DOWN", "M_cholesterol_SREBP"): -1, ("C4_DOWN", "M_TNF_NFkB"): -1, ("C4_DOWN", "M_p53_apoptosis"): -1,
 ("M_shared_antiproliferative", "M_cholesterol_SREBP"): -1, ("M_shared_antiproliferative", "M_TNF_NFkB"): -1,
 ("M_shared_antiproliferative", "M_p53_apoptosis"): -1, ("M_shared_ER_UPR", "M_cholesterol_SREBP"): -1,
 ("M_shared_ER_UPR", "M_TNF_NFkB"): -1, ("M_shared_ER_UPR", "M_p53_apoptosis"): -1}

def load_scores(tag):
    return pd.read_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t").set_index("Name")
T, G = load_scores("tcga"), load_scores("cgga")

def rho_mat(df, sets=ORDER):
    m = np.eye(len(sets))
    for i, a in enumerate(sets):
        for j, b in enumerate(sets):
            if j <= i:
                continue
            r, p = stats.spearmanr(df[a], df[b], nan_policy="omit")
            m[i, j] = m[j, i] = r
    return m

def plot_heatmap(ax, mat, title, lab):
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax.imshow(mat, cmap="RdBu_r", norm=norm)
    ax.set_xticks(range(len(lab))); ax.set_yticks(range(len(lab)))
    ax.set_xticklabels([LAB[k] for k in lab], rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels([LAB[k] for k in lab], fontsize=8)
    ax.set_title(title, fontsize=11)
    for i in range(len(lab)):
        for j in range(len(lab)):
            if i != j:
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(mat[i, j]) > 0.55 else "black")
    return im

def fisher_ci(rho, n):
    z = np.arctanh(np.clip(rho, -0.9999, 0.9999)); se = 1 / np.sqrt(n - 3); q = 1.96
    return np.tanh(z - q * se), np.tanh(z + q * se)

# ---------- Fig A/B: heatmaps ----------
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
im1 = plot_heatmap(axes[0], rho_mat(T), f"TCGA-GBM (n={len(T)})  Grade-A module correlation", ORDER)
im2 = plot_heatmap(axes[1], rho_mat(G), f"CGGA WHO IV (n={len(G)})  Grade-A module correlation", ORDER)
fig.colorbar(im1, ax=axes[0], shrink=0.8); fig.colorbar(im2, ax=axes[1], shrink=0.8)
fig.suptitle("Fig6 candidate A/B: Grade-A transcriptional-state architecture", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(FIG, "fig6_AB_gradeA_correlation_heatmap.png"), dpi=200)
plt.close(fig)

# ---------- Fig C: cross-cohort rho forest with 95% CI ----------
pairs = list(EXPECTED.keys())
rows = []
for a, b in pairs:
    rt, _ = stats.spearmanr(T[a], T[b], nan_policy="omit")
    rg, _ = stats.spearmanr(G[a], G[b], nan_policy="omit")
    lo_t, hi_t = fisher_ci(rt, len(T)); lo_g, hi_g = fisher_ci(rg, len(G))
    rows.append({"a": a, "b": b, "exp": EXPECTED[(a, b)], "rt": rt, "lo_t": lo_t, "hi_t": hi_t,
                 "rg": rg, "lo_g": lo_g, "hi_g": hi_g})
R = pd.DataFrame(rows)
R["pair"] = R.apply(lambda r: f"{LAB[r['a']]} vs {LAB[r['b']]}", axis=1)
R["match_t"] = np.sign(R.rt) == R.exp
R["match_g"] = np.sign(R.rg) == R.exp
R["both_match"] = R.match_t & R.match_g
R = R.sort_values(["exp", "both_match", "rt"], ascending=[False, False, True])
fig, ax = plt.subplots(figsize=(8.5, 11))
y = np.arange(len(R))
for i, (_, r) in enumerate(R.iterrows()):
    ax.errorbar([r.rt], [i + 0.16], xerr=[[r.rt - r.lo_t], [r.hi_t - r.rt]], fmt="o", ms=5,
                color=("#c0392b" if r.match_t else "#f5b7b1"), ecolor=("#c0392b" if r.match_t else "#f5b7b1"),
                elinewidth=1, capsize=2)
    ax.errorbar([r.rg], [i - 0.16], xerr=[[r.rg - r.lo_g], [r.hi_g - r.rg]], fmt="o", ms=5,
                color=("#1f618d" if r.match_g else "#aed6f1"), ecolor=("#1f618d" if r.match_g else "#aed6f1"),
                elinewidth=1, capsize=2)
    ax.text(1.02, i + 0.16, "P" if r.match_t else "X", va="center", fontsize=8, color="#7b241c")
    ax.text(1.08, i - 0.16, "P" if r.match_g else "X", va="center", fontsize=8, color="#1b4f72")
    ax.text(1.14, i, f"(exp {r.exp:+d})", va="center", fontsize=7, color="grey")
ax.plot([], [], "o", color="#c0392b", label="TCGA")
ax.plot([], [], "o", color="#1f618d", label="CGGA")
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(R["pair"], fontsize=7.5)
ax.invert_yaxis(); ax.set_xlim(-0.85, 1.35)
ax.set_xlabel("Spearman rho (95% CI)")
ax.set_title("Fig6 candidate C: cross-cohort effect-size forest plot\n(red=TCGA, blue=CGGA; P=expected direction hit, X=mismatch; exp +1/-1 = prereg expected sign)", fontsize=10)
ax.legend(loc="lower right"); fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig6_C_cross_cohort_rho_forest.png"), dpi=200)
plt.close(fig)

# ---------- Fig D: robustness ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
par = pd.read_csv(os.path.join(C3D, "_robust_partial_intermediate.tsv"), sep="\t")
ax = axes[0]
for coh, c in [("TCGA", "#c0392b"), ("CGGA", "#1f618d")]:
    s = par[par.cohort == coh]
    ax.scatter(s.rho_raw, s.rho_partial_global_mean, s=10, color=c, alpha=0.65, label=coh)
lim = [-1, 1]
ax.plot(lim, lim, ls="--", color="grey", lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel("raw rho"); ax.set_ylabel("partial rho (global expression)")
ax.set_title("D1: raw vs partial Spearman (Q1 pairs)", fontsize=10); ax.legend()
gp = pd.read_csv(os.path.join(C3D, "_robust_geneperturb_intermediate.tsv"), sep="\t")
ax = axes[1]
xt = np.arange(len(gp))
cols = [("#c0392b" if r.tag == "TCGA" else "#1f618d") for _, r in gp.iterrows()]
ax.bar(xt, gp.dir_match, color=cols, alpha=0.8)
ax.set_xticks(xt); ax.set_xticklabels([f"{r.tag[0].upper()}{r.run.split('_')[-1]}" if r.run.startswith("dropout") else r.run for _, r in gp.iterrows()], rotation=45, fontsize=7)
ax.axhline(10.5, color="grey", ls="--", lw=0.8)
ax.set_ylabel("direction-matched pairs (of 21 expected)")
ax.set_title("D2: gene-level perturbation (top-10% expr drop; 80% dropout x5)", fontsize=10)
fig.suptitle("Fig6 candidate D: robustness", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(FIG, "fig6_D_robustness.png"), dpi=200)
plt.close(fig)

# ---------- Fig E: conceptual (expected vs observed signed structure) ----------
exp_mat = np.zeros((len(ORDER), len(ORDER)))
obs_mat = np.zeros_like(exp_mat)
for i, a in enumerate(ORDER):
    for j, b in enumerate(ORDER):
        if i == j:
            continue
        if (a, b) in EXPECTED:
            exp_mat[i, j] = EXPECTED[(a, b)]
        elif (b, a) in EXPECTED:
            exp_mat[i, j] = EXPECTED[(b, a)]
        rt, _ = stats.spearmanr(T[a], T[b], nan_policy="omit")
        rg, _ = stats.spearmanr(G[a], G[b], nan_policy="omit")
        obs_mat[i, j] = (rt + rg) / 2
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
for ax, mat, tt in [(axes[0], exp_mat, "Pre-registered expectation"), (axes[1], obs_mat, "Observed mean rho (TCGA+CGGA)")]:
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax.imshow(mat, cmap="RdBu_r", norm=norm)
    ax.set_xticks(range(len(ORDER))); ax.set_yticks(range(len(ORDER)))
    ax.set_xticklabels([LAB[k] for k in ORDER], rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels([LAB[k] for k in ORDER], fontsize=8)
    ax.set_title(tt, fontsize=11)
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            if i != j and mat[i, j] != 0:
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(mat[i, j]) > 0.55 else "black")
fig.colorbar(im, ax=axes, shrink=0.8)
fig.suptitle("Fig6 candidate E: conceptual summary - expected vs observed signed architecture", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(FIG, "fig6_E_conceptual_expected_vs_observed.png"), dpi=200)
plt.close(fig)

# ---------- QC overview: coverage & C5 axis ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
covt = json.load(open(os.path.join(C3D, "tcga_set_coverage.json")))
covg = json.load(open(os.path.join(C3D, "cgga_set_coverage.json")))
ax = axes[0]
keys = sorted(set(covt) | set(covg))
xx = np.arange(len(keys)); w = 0.4
ax.bar(xx - w / 2, [covt.get(k, 0) for k in keys], w, label="TCGA", color="#c0392b", alpha=0.8)
ax.bar(xx + w / 2, [covg.get(k, 0) for k in keys], w, label="CGGA", color="#1f618d", alpha=0.8)
ax.axhline(0.80, color="grey", ls="--", lw=0.8)
ax.set_xticks(xx); ax.set_xticklabels(keys, rotation=60, fontsize=7)
ax.set_ylabel("fraction of frozen set present in matrix")
ax.set_title("QC: frozen-set gene coverage by cohort", fontsize=10)
ax.legend()
ax = axes[1]
c5axd = pd.read_csv(os.path.join(C3D, "_c5_axis_intermediate.tsv"), sep="\t")
for coh, c, n in [("TCGA", "#c0392b", len(T)), ("CGGA", "#1f618d", len(G))]:
    s = c5axd[c5axd.cohort == coh]
    lo = [np.tanh(np.arctanh(np.clip(x, -0.9999, 0.9999)) - 1.96 / np.sqrt(n - 3)) for x in s.C4UP_rho]
    hi = [np.tanh(np.arctanh(np.clip(x, -0.9999, 0.9999)) + 1.96 / np.sqrt(n - 3)) for x in s.C4UP_rho]
    ax.errorbar(s.C4UP_rho, s.score, xerr=[np.array(s.C4UP_rho) - np.array(lo), np.array(hi) - np.array(s.C4UP_rho)],
                fmt="o", ms=6, color=c, label=coh)
ax.set_xlabel("Spearman rho with C4_UP (HITMAN-enriched axis)"); ax.set_ylabel("C5 signature")
ax.set_title("QC: C5 signatures vs C4_UP axis (exploratory)", fontsize=10)
ax.legend()
fig.suptitle("Fig6/QC: coverage & exploratory C5 contextualization", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(FIG, "qc_coverage_C5axis.png"), dpi=200)
plt.close(fig)

# figlist md
with open(os.path.join(FIG, "figlist.md"), "w", encoding="utf-8") as f:
    f.write("# C3_FIG6_CANDIDATE\n\n候选图（只读真实结果，不做美化；失败的负向预期在图中如实显示）。\n\n")
    f.write("- fig6_AB_gradeA_correlation_heatmap.png  Q1 architecture（TCGA / CGGA）\n")
    f.write("- fig6_C_cross_cohort_rho_forest.png      跨队列 effect-size forest plot（P/X 标预期方向命中）\n")
    f.write("- fig6_D_robustness.png                   partial correlation + gene-perturbation\n")
    f.write("- fig6_E_conceptual_expected_vs_observed.png  预期 vs 观察 signed structure\n")
    f.write("- qc_coverage_C5axis.png                  QC：frozen-set coverage 与 C5 exploratory axis\n")
print("figures written")
for f in sorted(os.listdir(FIG)):
    print(f, round(os.path.getsize(os.path.join(FIG, f)) / 1024, 1), "KB")
