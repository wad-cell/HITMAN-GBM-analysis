# -*- coding: utf-8 -*-
"""Fig5 - replication, robustness and retained negatives. Frozen tables only."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from fig_common import *

ev = rd(os.path.join(PT, "C2_11_EVIDENCE_GRADES.tsv"))
lopo_d = rd(os.path.join(TP, "ss2_lopo.tsv"))
lopo_la = rd(os.path.join(PT, "C2_couturier_layerA_lopo.tsv"))
null_cell = rd(os.path.join(TP, "ss2_cellnull_ALL.tsv"))
null_pb = rd(os.path.join(TP, "ss2_pseudobulk_null.tsv"))
ccdep = {r["signature"]: r for r in rd(os.path.join(PT, "C2_3_CC_DEPLETION_SUMMARY.tsv"))}

GRADE_COL = {"A": "#1e8449", "B": "#b9770e", "C": "#7f8c8d", "D": "#a6a6a6", "E": "#c0392b"}
fig = plt.figure(figsize=(7.2, 15.0))
gs = fig.add_gridspec(4, 1, height_ratios=[1.35, 2.30, 2.30, 0.85], hspace=0.30)

# ---- A: Layer A replication ----
axA = fig.add_subplot(gs[0]); panel_label(axA, "A", dx=-0.19, dy=1.02)
labs, deltas, fdrv, cols = [], [], [], []
for r in ev:
    d = (r["layerA_delta"] or "").split("/")[0]
    f = (r["layerA_fdr"] or "").split("/")[0]
    lbl = f"{r['signature']}:{r['anchor_state']}"
    labs.append(lbl); deltas.append((num(d) or 0.0) * 1000); fdrv.append(num(f))
    cols.append(GRADE_COL.get(r["final_evidence_grade"].strip(), "#7f8c8d"))
y = np.arange(len(labs))[::-1]
axA.barh(y, deltas, color=cols, height=0.72)
for yi, v, f in zip(y, deltas, fdrv):
    axA.text(v + (0.4 if v >= 0 else -0.4), yi, f"FDR {f:.3g}" + star(f), va="center",
             ha="left" if v >= 0 else "right", fontsize=5.6)
axA.set_yticks(y); axA.set_yticklabels(labs, fontsize=5.4)
axA.set_xlabel("Layer A mean patient-level AUC difference (x 10^-3)", fontsize=7)
axA.axvline(0, color="black", lw=0.7)
axA.set_xlim(-4, 46)
axA.set_title("Independent Layer A replication of the anchor programmes (bars coloured by frozen evidence grade)")
axA.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=GRADE_COL[g]) for g in "ABCDE"],
           labels=[f"grade {g}" for g in "ABCDE"], frameon=False, fontsize=5.8, ncol=5, loc="lower right")

# ---- B: LOPO stability ----
NOTE_B = []
gsB = gs[1].subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.55)
for k, (tab, key, ttl, note) in enumerate([
        (lopo_d, "loo_flips_top", "Discovery (20 patients, leave-one-out)",
         "one top-state change (CC-depleted C4_DOWN axis)"),
        (lopo_la, "flip", "Layer A (12 leave-one-out runs)", "all anchors stable")]):
    ax = fig.add_subplot(gsB[k])
    if k == 0:
        panel_label(ax, "B", dx=-0.42, dy=1.02)
    vals = [num(r[key]) or 0 for r in tab]
    yy = np.arange(len(tab))[::-1]
    flag = [str(r.get("flag_patient_sensitive", "")).strip() for r in tab]
    cols1 = ["#c0392b" if (v or 0) > 0 else "#1e8449" for v in vals]
    ax.barh(yy, vals, color=cols1, height=0.72)
    for yi, v in zip(yy, vals):
        ax.text(v + 0.06, yi, str(int(v)), va="center", fontsize=5.6)
    ax.set_yticks(yy); ax.set_yticklabels([r["sig"] for r in tab], fontsize=5.0)
    ax.set_xlabel("leave-one-patient-out top-state changes", fontsize=6.4)
    ax.set_xlim(0, max(2, max(vals) + 0.6))
    ax.set_title(ttl, fontsize=7)
    NOTE_B.append("%s (n = %d): %s" % (ttl.split(" (")[0], len(tab), note))

# ---- C: random gene-set null controls ----
axC = fig.add_subplot(gs[2]); panel_label(axC, "C", dx=-0.19, dy=1.02)
sigs = []
for r in null_cell:
    if r["sig"] not in sigs:
        sigs.append(r["sig"])
STATES = ["MES", "AC", "OPC", "NPC", "Cycling"]
Z = {(r["sig"], r["state"]): num(r["z"]) for r in null_cell}
P = {(r["sig"], r["state"]): num(r["emp_p"]) for r in null_cell}
M = np.array([[Z.get((s, st)) or 0.0 for st in STATES] for s in sigs])
im = axC.imshow(M, cmap="viridis", aspect="auto")
for i, s in enumerate(sigs):
    for j, st in enumerate(STATES):
        p = P.get((s, st))
        if p is not None:
            axC.text(j, i, f"{M[i, j]:.1f}" + ("*" if p < 0.05 else ""), ha="center", va="center",
                     fontsize=5.6, color="white" if M[i, j] < (M.max() * 0.6) else "black")
axC.set_xticks(range(len(STATES))); axC.set_xticklabels(STATES, fontsize=6.2)
axC.set_yticks(range(len(sigs))); axC.set_yticklabels(sigs, fontsize=5.0)
axC.set_title("Random gene-set null controls: observed vs size-matched random sets (cell-level, z; * empirical P < 0.05)")
fig.colorbar(im, ax=axC, shrink=0.7, label="z")
NOTE_C5 = ("Pseudobulk-level null controls (ss2_pseudobulk_null.tsv) are reported for the same signatures; "
           "cell-cycle-depleted variants of the C5 signatures remove %s%%-%s%% of genes (C5_UP/C5_DOWN)."
           % (ccdep["C5_UP"]["removed_fraction_pct"], ccdep["C5_DOWN"]["removed_fraction_pct"]))

# ---- D: retained negatives ----
axD = fig.add_subplot(gs[3]); panel_label(axD, "D", dx=-0.19, dy=1.05)
axD.set_axis_off()
lines = [
 "Retained negative and sensitivity findings carried into the submitted figures:",
 "  C4_DOWN and shared anti-proliferative mappings attenuate under expanded-state mapping (FDR 0.594 / 0.249); primary anchors retained with this caveat.",
 "  C5 base Layer A MES contrast not significant (FDR 0.0563) -> Grade C, incompletely replicated.",
 "  Autophagy-lysosome (MES) not significant in TCGA bulk (n = 44 MES, P = 0.334) -> Grade B sensitivity only.",
 "  HSF1 heat-shock reversible between cohorts (Discovery MES FDR 0.002 vs Layer A MES FDR 1.000) -> Grade E, reported as heterogeneity.",
 "  CGGA stratification by IDH status and primary/recurrent status does not improve expected-direction recovery (9/21 in every stratum).",
 "  Gene-level cross-platform concordance with the comparators is absent (|rho| <= 0.14); shared programmes are a module-level concept only.",
 "  C5 does not reproduce on conventional electric-field platforms (gene-level and module-level checks) and remains HITMAN-field-specific.",
]
_fsD5, _nD5 = fit_rows_in_ax(axD, lines, fs_max=6.4, fs_min=5.0, top=0.98, lead=1.28, row_gap=0.40)
print("Fig5 panel D fontsize/rows: %.1f / %d" % (_fsD5, _nD5))
axD.set_title("Retained negative findings (no element is removed or smoothed away)", fontsize=8)

from fig_common import add_caption
add_caption(fig, ["Panel B (leave-one-out top-state changes): " + "; ".join(NOTE_B) + ".",
                  "Panel C: " + NOTE_C5], fontsize=6.0, width=150)

save_fig(fig, "Fig5")
log_values("Fig5", [("A", "layerA_delta", "C2_11_EVIDENCE_GRADES.tsv", f"{r['signature']}|{r['layerA_delta']}|{r['layerA_fdr']}") for r in ev]
          + [("B", "loo_flips_top", "ss2_lopo.tsv", f"{r['sig']}={r['loo_flips_top']}") for r in lopo_d]
          + [("B", "flip", "C2_couturier_layerA_lopo.tsv", f"{r['sig']}={r['flip']}") for r in lopo_la])
print("Fig5 done")
