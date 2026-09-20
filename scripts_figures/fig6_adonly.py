# -*- coding: utf-8 -*-
"""Fig6 (FINAL-LEGEND COMPLIANT: panels A-D only) - TCGA/CGGA human GBM contextualization.

Layout-only variant of fig6.py: panel E is NOT rendered, because the FINAL legend
defines Fig 6 as A-D.  The content that the legend assigns to panel D (coverage,
stratification, gene-level inconsistency) is kept inside panel D as two stacked
sub-axes with a single "D" label.  Frozen tables only; no statistics recomputed.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from fig_common import *

pr = rd(os.path.join(PT, "C3_PRIMARY_RESULTS.tsv"))
rb = rd(os.path.join(PT, "C3_ROBUSTNESS.tsv"))
sec = rd(os.path.join(PT, "C3_SECONDARY_RESULTS.tsv"))
q2 = rd(os.path.join(C3D, "_q2_for_output.tsv"))

q1 = [r for r in pr if r["block"] == "Q1_architecture"]
SHORT = {"C4_UP": "C4_UP", "C4_DOWN": "C4_DOWN", "M_shared_antiproliferative": "antiprolif.",
         "M_shared_ER_UPR": "ER-UPR", "M_cholesterol_SREBP": "SREBP", "M_TNF_NFkB": "TNF/NF-kB",
         "M_p53_apoptosis": "p53", "M_RIGI_typeI_IFN": "RIG-I"}
mods = []
for r in q1:
    for m in (r["module_a"], r["module_b"]):
        if m not in mods:
            mods.append(m)
mods = [m for m in ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
                    "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"] if m in mods]
idx = {m: i for i, m in enumerate(mods)}


def matrix(cohort):
    M = np.full((len(mods), len(mods)), np.nan)
    for r in q1:
        if r["cohort"] != cohort:
            continue
        i, j = idx[r["module_a"]], idx[r["module_b"]]
        v = num(r["effect_rho"])
        M[i, j] = M[j, i] = v
    return M


RV = {("C4_UP", "C4_DOWN"), ("M_shared_antiproliferative", "M_TNF_NFkB"),
      ("M_shared_antiproliferative", "M_p53_apoptosis")}

fig = plt.figure(figsize=(7.4, 13.2))
gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.05, 2.6, 1.30], hspace=0.72)

# ---- A / B heatmaps ----
for k, cohort in enumerate(["TCGA", "CGGA"]):
    ax = fig.add_subplot(gs[k]); panel_label(ax, "AB"[k], dx=-0.24, dy=1.03)
    M = matrix(cohort)
    im = ax.imshow(np.ma.masked_invalid(M), cmap="coolwarm", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1))
    for i in range(len(mods)):
        for j in range(len(mods)):
            if i == j or np.isnan(M[i, j]):
                continue
            rev = (mods[i], mods[j]) in RV or (mods[j], mods[i]) in RV
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.6,
                    fontweight="bold" if rev else "normal", color="black")
    ax.set_xticks(range(len(mods))); ax.set_xticklabels([SHORT[m] for m in mods], rotation=45, ha="right", fontsize=6.2)
    ax.set_yticks(range(len(mods))); ax.set_yticklabels([SHORT[m] for m in mods], fontsize=6.2)
    n = 288 if cohort == "TCGA" else 249
    ax.set_title(f"({'AB'[k]}) {cohort} relationship structure (n = {n}); Spearman rho between module scores")
    fig.colorbar(im, ax=ax, shrink=0.75, label="rho")
    if k == 0:
        NOTE_A = ("Contextualisation only: these are cross-sectional correlations in bulk tumours, not causal or "
                  "replication evidence. Bold values mark the three sign-reversed pairs.")

# inset for B: TCGA MES vs non-MES
axB = fig.axes[1]
axin = axB.inset_axes([0.58, 0.26, 0.40, 0.60])
ins = [(r["module"], num(r["median_diff_MES_minus_nonMES"]), num(r["P"])) for r in q2]
aut = [r for r in sec if r["module"] == "M_autophagy_lysosome" and r["cohort"].startswith("TCGA")]
if aut:
    ins.append(("M_autophagy_lysosome", num(aut[0]["median_diff_MES_minus_nonMES"]), num(aut[0]["P"])))
ins = sorted(ins, key=lambda t: t[1])
yy = np.arange(len(ins))
cols = ["#7f8c8d" if (t[2] or 1) >= 0.05 else "#1e8449" for t in ins]
axin.barh(yy, [t[1] for t in ins], color=cols)
axin.set_yticks(yy); axin.set_yticklabels([SHORT.get(t[0], t[0]) for t in ins], fontsize=4.6)
axin.set_ylim(-0.6, len(ins) - 0.4)
axin.axvline(0, color="black", lw=0.6)
for yi, t in zip(yy, ins):
    axin.text(t[1] + 60, yi, ("n.s." if (t[2] or 1) >= 0.05 else f"P={t[2]:.1e}"), fontsize=4.6, va="center")
axin.set_xticks([]); axin.tick_params(axis="x", length=0)
axin.set_title("TCGA MES (n = 15) vs non-MES (n = 29)", fontsize=5.6)

# ---- C forest plot ----
axC = fig.add_subplot(gs[2]); panel_label(axC, "C", dx=-0.20, dy=1.02)
f = [r for r in rb if r["analysis"] == "robust_1_cross_cohort_sign"]
cnt = [r for r in f if r["module_a"] == "SUMMARY"]
f = [r for r in f if r["module_a"] != "SUMMARY"]
yy = np.arange(len(f))[::-1]
for yi, r in zip(yy, f):
    a, b = num(r["rho_TCGA"]), num(r["rho_CGGA"])
    rev = str(r["same_direction"]).strip() == "False"
    col = "#c0392b" if rev else "#2c3e50"
    axC.plot([a, b], [yi, yi], color=col, lw=0.9, alpha=0.85)
    axC.plot([a], [yi], "o", ms=3.4, color="#e67e22")
    axC.plot([b], [yi], "s", ms=3.2, color="#2980b9")
    if str(r["expected_sign"]).strip():
        axC.plot([num(r["expected_sign"])], [yi], "x", ms=3.6, color="black")
axC.axvline(0, color="black", lw=0.8)
axC.set_yticks(yy)
axC.set_yticklabels([f"{SHORT.get(r['module_a'], r['module_a'])}-{SHORT.get(r['module_b'], r['module_b'])}" for r in f],
                    fontsize=5.0)
axC.set_xlabel("Spearman rho", fontsize=7)
axC.set_xlim(-0.75, 1.0)
axC.set_title("Cross-cohort forest plot of the 28 module relationships (circles: TCGA; squares: CGGA; x: expected sign; "
              "red: 3 sign-reversed pairs)")
axC.legend(handles=[plt.Line2D([], [], color="#e67e22", marker="o", ls="", ms=4, label="TCGA"),
                    plt.Line2D([], [], color="#2980b9", marker="s", ls="", ms=4, label="CGGA"),
                    plt.Line2D([], [], color="black", marker="x", ls="", ms=4, label="expected sign"),
                    plt.Line2D([], [], color="#c0392b", lw=1.2, label="sign reversal")],
           frameon=False, fontsize=5.8, ncol=4, loc="lower left")
note = cnt[0]["note"] if cnt else ""
NOTE_C = f"Panel C summary (robust_1_cross_cohort_sign): {note}"

# ---- D robustness (single letter "D"; two stacked sub-axes, no panel E) ----
gsD = gs[3].subgridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.62)
st = [r for r in rb if r["analysis"] == "robust_2_stratum_Q1"]
gp = [r for r in rb if r["analysis"] == "robust_4_gene_perturbation"]

axD = fig.add_subplot(gsD[0]); panel_label(axD, "D", dx=-0.30, dy=1.03)
labs, vals = [], []
for r in st:
    labs.append(r["cohort_scope"].replace("CGGA_", "CGGA ").replace("TCGA_", "TCGA ").replace("_", " "))
    vals.append(num(str(r["expected_sign"]).split("=")[-1]) or 0)
yy = np.arange(len(labs))[::-1]
axD.barh(yy, vals, color="#5dade2", height=0.66)
for yi, v in zip(yy, vals):
    axD.text(v + 0.15, yi, f"{int(v)}/21", va="center", fontsize=5.8)
axD.set_yticks(yy); axD.set_yticklabels(labs, fontsize=5.8)
axD.set_xlim(0, 14); axD.set_xlabel("expected-direction relationships recovered (/21)", fontsize=6.4)
axD.set_title("Robustness and negative results: stratification (IDH, primary/recurrent)")
axD.set_ylabel("")

axD2 = fig.add_subplot(gsD[1])
for k, coh in enumerate(["tcga", "cgga"]):
    rs = [r for r in gp if r["cohort_scope"] == coh]
    v = [num(str(r["expected_sign"]).split("=")[-1]) or 0 for r in rs]
    axD2.plot(range(len(v)), v, "o-", ms=3.6, lw=1.0, color=["#e67e22", "#2980b9"][k], label=coh.upper())
axD2.set_xticks(range(7))
axD2.set_xticklabels(["base", "drop10", "d80-1", "d80-2", "d80-3", "d80-4", "d80-5"], fontsize=5.2, rotation=45)
axD2.set_ylabel("recovered (/21)", fontsize=6.4); axD2.set_ylim(6, 15)
cov = [r for r in rb if r["analysis"] == "robust_5_coverage"]
axD2.set_title("gene-perturbation sensitivity (no imputation); coverage rule >= 0.80", fontsize=6.6, loc="left")
axD2.legend(frameon=False, fontsize=5.8, loc="upper right")
axD2.text(0.0, -0.55, "Coverage: " + "; ".join(f"{r['cohort_scope'].upper()} {r['expected_sign']}" for r in cov),
          transform=axD2.transAxes, fontsize=6.0, va="top")

add_caption(fig, ["Panel A/B note: " + NOTE_A, NOTE_C], fontsize=6.0, width=150)

save_fig(fig, "Fig6")
# plotted values for Fig6 are already recorded in provenance/plotted_values/pre_release/Fig6_values.tsv
# (A-D content is identical to the pre-fix export); not re-appended here to avoid duplicate rows.
print("Fig6 (A-D, FINAL legend compliant) done")
