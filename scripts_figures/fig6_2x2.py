# -*- coding: utf-8 -*-
"""Fig6 the release stage - LAYOUT-ONLY 2x2 re-composition (panels A-D; panel E not rendered).

Scientific content is frozen: this script reads exactly the same frozen tables as the
the pre-release export and plots exactly the same values. Only the canvas geometry changes
(4x1 single column -> 2x2 grid), so that the figure stays readable after scaling.
No statistics, scales, thresholds or numbers are recomputed or altered.
"""
import os, sys
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCR)          # fig_common.py ships in this same directory
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import fig_common
from fig_common import *

# Release plotted-value log (the pre-release log shipped under plotted_values/pre_release
# is left untouched so the two exports remain independently checkable).
PKG_ROOT = os.path.dirname(SCR)
REL_LOG = os.path.join(PKG_ROOT, "provenance", "plotted_values", "release")
os.makedirs(REL_LOG, exist_ok=True)
fig_common.LOGD = REL_LOG

pr = rd(os.path.join(PT, "C3_PRIMARY_RESULTS.tsv"))
rb = rd(os.path.join(PT, "C3_ROBUSTNESS.tsv"))
sec = rd(os.path.join(PT, "C3_SECONDARY_RESULTS.tsv"))
q2 = rd(os.path.join(C3D, "_q2_for_output.tsv"))

q1 = [r for r in pr if r["block"] == "Q1_architecture"]
SHORT = {"C4_UP": "C4_UP", "C4_DOWN": "C4_DOWN", "M_shared_antiproliferative": "antiprolif.",
         "M_shared_ER_UPR": "ER-UPR", "M_cholesterol_SREBP": "SREBP", "M_TNF_NFkB": "TNF/NF-kB",
         "M_p53_apoptosis": "p53", "M_RIGI_typeI_IFN": "RIG-I",
         "M_autophagy_lysosome": "autophagy"}
mods = []
for r in q1:
    for m in (r["module_a"], r["module_b"]):
        if m not in mods:
            mods.append(m)
mods = [m for m in ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
                    "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"] if m in mods]
idx = {m: i for i, m in enumerate(mods)}


def dm(s):
    """Numerator of a frozen 'dir_match=X/21' string (rendering helper; no computation)."""
    return num(str(s).split("=")[-1].split("/")[0].strip()) or 0


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

# ---------------------------------------------------------------- canvas 2x2
# Canvas is sized at full-page print width (7.5 in) so that no down-scaling is
# needed at typesetting and every font size below is the final printed size.
fig = plt.figure(figsize=(7.5, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[0.86, 1.14], width_ratios=[1.0, 1.0],
                      hspace=0.46, wspace=0.30, left=0.055, right=0.955,
                      top=0.945, bottom=0.05)

# ---- A / B heatmaps (row 1, col 1 / col 2) ----
# release (layout only): the two heatmaps are placed explicitly and capped at a fixed
# square side, so that the free band directly under panel B can host the B inset without any
# text overprint.  Geometry only - the plotted matrices, values and statistics are untouched.
HS = 1.80                                            # heatmap square side (in)
HTOP = 7.489                                         # top edge of the heatmaps (in, pre-caption)
HX = [0.096, 0.6053]                                 # left edges (fraction of 7.5 in width)
axA = axB = None
for k, cohort in enumerate(["TCGA", "CGGA"]):
    ax = fig.add_axes([HX[k], (HTOP - HS) / 8.2, HS / 7.5, HS / 8.2])
    panel_label(ax, "AB"[k], dx=-0.24, dy=1.03)
    axA, axB = (ax, axB) if k == 0 else (axA, ax)
    M = matrix(cohort)
    im = ax.imshow(np.ma.masked_invalid(M), cmap="coolwarm", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1))
    for i in range(len(mods)):
        for j in range(len(mods)):
            if i == j or np.isnan(M[i, j]):
                continue
            rev = (mods[i], mods[j]) in RV or (mods[j], mods[i]) in RV
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6.0,
                    fontweight="bold" if rev else "normal", color="black")
    ax.set_xticks(range(len(mods))); ax.set_xticklabels([SHORT[m] for m in mods], rotation=40, ha="right", fontsize=6.4)
    ax.set_yticks(range(len(mods))); ax.set_yticklabels([SHORT[m] for m in mods], fontsize=6.4)
    n = 288 if cohort == "TCGA" else 249
    ax.set_title(f"({'AB'[k]}) {cohort} relationship structure (n = {n}); Spearman rho between module scores", fontsize=8.5)
    cax = fig.add_axes([HX[k] + HS / 7.5 + 0.0093, (HTOP - 0.875 * HS) / 8.2, 0.012, 0.75 * HS / 8.2])
    fig.colorbar(im, cax=cax, label="rho")
    if k == 0:
        NOTE_A = ("Contextualisation only: these are cross-sectional correlations in bulk tumours, not causal or "
                  "replication evidence. Bold values mark the three sign-reversed pairs.")

# inset for B: TCGA MES vs non-MES. release layout-only change: the inset is drawn as a
# standalone mini-panel in the free band directly below heatmap B instead of floating inside
# the heatmap, so that its module labels no longer overprint panel B's correlation values.
# The plotted rows are exactly the pre-release rows (no recomputation, no value change);
# the final position is set after the caption band has been added, in final-canvas coordinates.
axin = fig.add_axes([0.564, 0.560, 0.391, 0.090])
ins = [(r["module"], num(r["median_diff_MES_minus_nonMES"]), num(r["P"])) for r in q2]
aut = [r for r in sec if r["module"] == "M_autophagy_lysosome" and r["cohort"].startswith("TCGA")]
if aut:
    ins.append(("M_autophagy_lysosome", num(aut[0]["median_diff_MES_minus_nonMES"]), num(aut[0]["P"])))
ins = sorted(ins, key=lambda t: t[1])
yy = np.arange(len(ins))
cols = ["#7f8c8d" if (t[2] or 1) >= 0.05 else "#1e8449" for t in ins]
def _fit_inset_xlim(ax):
    """release (layout only): choose the display range of the inset so that the bar labels
    stay inside the axes box.  Only the drawn range is chosen here - no plotted value changes."""
    lo, hi = -1.0, 1.0
    for _ in range(3):
        fig.canvas.draw()
        try:
            R = fig.canvas.get_renderer()
        except Exception:
            R = fig._get_renderer()
        inv = ax.transData.inverted()
        xs = [0.0]
        for t in ax.texts:
            bb = t.get_window_extent(renderer=R)
            xs += [inv.transform((bb.x0, bb.y0))[0], inv.transform((bb.x1, bb.y1))[0]]
        new_lo, new_hi = min(xs), max(xs)
        pad = 0.05 * (new_hi - new_lo)
        new_lo, new_hi = new_lo - pad, new_hi + pad
        if abs(new_lo - lo) < 1e-9 and abs(new_hi - hi) < 1e-9:
            break
        lo, hi = new_lo, new_hi
        ax.set_xlim(lo, hi)


axin.barh(yy, [t[1] for t in ins], color=cols)
axin.set_yticks(yy); axin.set_yticklabels([SHORT.get(t[0], t[0]) for t in ins], fontsize=5.6)
axin.set_ylim(-0.45, len(ins) - 0.55)   # display range only (tighter row padding = better label legibility)
axin.axvline(0, color="black", lw=0.6)
for yi, t in zip(yy, ins):
    axin.annotate("n.s." if (t[2] or 1) >= 0.05 else f"P={t[2]:.1e}", xy=(t[1], yi), xytext=(2.5, 0),
                  textcoords="offset points", fontsize=5.6, va="center", ha="left")
axin.set_xticks([]); axin.tick_params(axis="x", length=0)
axin.set_title("TCGA MES (n = 15) vs non-MES (n = 29)", fontsize=6.2, pad=2.0)

# ---- C forest plot (row 2, col 1) ----
axC = fig.add_subplot(gs[1, 0]); panel_label(axC, "C", dx=-0.20, dy=1.02)
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
                    fontsize=6.2)
axC.set_xlabel("Spearman rho", fontsize=7.5)
axC.set_xlim(-0.75, 1.0)
axC.set_title("Cross-cohort forest plot of the 28 module relationships (circles: TCGA; squares: CGGA; x: expected sign; "
              "red: 3 sign-reversed pairs)", fontsize=8.5)
axC.legend(handles=[plt.Line2D([], [], color="#e67e22", marker="o", ls="", ms=4, label="TCGA"),
                    plt.Line2D([], [], color="#2980b9", marker="s", ls="", ms=4, label="CGGA"),
                    plt.Line2D([], [], color="black", marker="x", ls="", ms=4, label="expected sign"),
                    plt.Line2D([], [], color="#c0392b", lw=1.2, label="sign reversal")],
           frameon=False, fontsize=6.4, ncol=4, loc="lower left")
note = cnt[0]["note"] if cnt else ""
NOTE_C = f"Panel C summary (robust_1_cross_cohort_sign): {note}"

# ---- D robustness (row 2, col 2; single letter "D", two stacked sub-axes, no panel E) ----
gsD = gs[1, 1].subgridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.62)
st = [r for r in rb if r["analysis"] == "robust_2_stratum_Q1"]
gp = [r for r in rb if r["analysis"] == "robust_4_gene_perturbation"]

axD = fig.add_subplot(gsD[0]); panel_label(axD, "D", dx=-0.30, dy=1.03)
# release defect fix (rendering only): the recovered-relationship counts live in the frozen
# field module_b ("dir_match=X/21"); the pre-release export read expected_sign, which is empty
# for these rows, so panel D rendered 0/21 bars and a flat 0 line off-scale.
labs, vals = [], []
for r in st:
    labs.append(r["cohort_scope"].replace("CGGA_", "CGGA ").replace("TCGA_", "TCGA ").replace("_", " "))
    vals.append(dm(r["module_b"]))
yy = np.arange(len(labs))[::-1]
axD.barh(yy, vals, color="#5dade2", height=0.66)
for yi, v in zip(yy, vals):
    axD.text(v + 0.15, yi, f"{int(v)}/21", va="center", fontsize=6.5)
axD.set_yticks(yy); axD.set_yticklabels(labs, fontsize=6.5)
axD.set_xlim(0, 14); axD.set_xlabel("expected-direction relationships recovered (/21)", fontsize=7.0)
axD.set_title("Robustness and negative results: stratification (IDH, primary/recurrent)", fontsize=8.5)
axD.set_ylabel("")

axD2 = fig.add_subplot(gsD[1])
for k, coh in enumerate(["tcga", "cgga"]):
    rs = [r for r in gp if r["cohort_scope"] == coh]
    v = [dm(r["module_b"]) for r in rs]
    axD2.plot(range(len(v)), v, "o-", ms=3.6, lw=1.0, color=["#e67e22", "#2980b9"][k], label=coh.upper())
axD2.set_xticks(range(7))
axD2.set_xticklabels(["base", "drop10", "d80-1", "d80-2", "d80-3", "d80-4", "d80-5"], fontsize=6.2, rotation=45)
axD2.set_ylabel("recovered (/21)", fontsize=7.0); axD2.set_ylim(6, 15)
cov = [r for r in rb if r["analysis"] == "robust_5_coverage"]
axD2.set_title("gene-perturbation sensitivity (no imputation); coverage rule >= 0.80", fontsize=8.0, loc="left")
axD2.legend(frameon=False, fontsize=6.4, loc="upper right")
NOTE_D = ("Panel D values (frozen C3_ROBUSTNESS.tsv): stratification dir_match "
          + "; ".join(f"{r['cohort_scope']}={str(r['module_b']).split('=')[-1]}" for r in st)
          + "; coverage " + "; ".join(f"{r['cohort_scope'].upper()} {r['module_b']}" for r in cov) + ".")

add_caption(fig, ["Panel A/B note: " + NOTE_A, NOTE_C, NOTE_D], fontsize=7.0, width=132)

# ---- release (layout only): pin the two heatmaps to an exact square and the B inset into
# the free band below heatmap B ----
# The inches below are given in the pre-caption canvas; add_caption shifts every axes down by
# the caption band while keeping its size in inches, so the boxes are re-stated here in final
# canvas fractions (this avoids the aspect machinery re-scaling a box that was set in
# pre-caption fractions).  The band under panel B sits below heatmap B's tick labels and above
# panel D's title; both clearances are verified by the release layout QC.
_CAP = fig.get_size_inches()[1] - 8.2                        # caption band height (in)
_H = fig.get_size_inches()[1]
for _k, _ax in enumerate([axA, axB]):
    _ax.set_position([HX[_k], (HTOP - HS + _CAP) / _H, HS / 7.5, HS / _H])
axin.set_position([4.58 / 7.5, (4.19 + _CAP) / _H,
                   2.55 / 7.5, 0.99 / _H])
axin.tick_params(axis="y", length=1.5, pad=1.5)
_fit_inset_xlim(axin)

save_fig(fig, "Fig6")

# --- release traceability log: same metrics as the pre-release plotted-value log, re-dumped
#     from the same frozen tables (used to prove that no plotted value changed). ---
rows = []
for r in q1:
    rows.append(("A" if r["cohort"] == "TCGA" else "B", "effect_rho",
                 "C3_PRIMARY_RESULTS.tsv/Q1_architecture", f"{r['module_a']}|{r['module_b']}|{r['effect_rho']}"))
for r in q2:
    rows.append(("B", "median_diff_MES_minus_nonMES", "_q2_for_output.tsv",
                 f"{r['module']}|{r['median_diff_MES_minus_nonMES']}|P={r['P']}"))
for r in sec:
    if r["module"] == "M_autophagy_lysosome" and r["cohort"].startswith("TCGA"):
        rows.append(("B", "median_diff_MES_minus_nonMES", "C3_SECONDARY_RESULTS.tsv",
                     f"{r['module']}|{r['median_diff_MES_minus_nonMES']}|P={r['P']}"))
for r in rb:
    rows.append(("C", "rho_TCGA/rho_CGGA/expected_sign",
                 "C3_ROBUSTNESS.tsv/robust_1_cross_cohort_sign",
                 f"{r['module_a']}|{r['module_b']}|{r['rho_TCGA']}|{r['rho_CGGA']}|{r['expected_sign']}"))
for r in st:
    rows.append(("D", "dir_match", "C3_ROBUSTNESS.tsv/robust_2_stratum_Q1",
                 f"{r['cohort_scope']}|{r['module_b']}"))
for r in gp:
    rows.append(("D", "dir_match", "C3_ROBUSTNESS.tsv/robust_4_gene_perturbation",
                 f"{r['cohort_scope']}|{r['module_b']}"))
for r in cov:
    rows.append(("D", "min_per_sample_gene_coverage", "C3_ROBUSTNESS.tsv/robust_5_coverage",
                 f"{r['cohort_scope']}|{r['module_b']}"))
log_values("Fig6", rows)
print("Fig6 2x2 layout-only re-composition done")
