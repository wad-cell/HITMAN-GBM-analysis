# -*- coding: utf-8 -*-
"""Fig4 - patient-level state localisation, grades and retained negatives. Frozen tables only."""
import os, sys, collections, textwrap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from fig_common import *

ev = rd(os.path.join(PT, "C2_11_EVIDENCE_GRADES.tsv"))
disc = {(r["sig"], r["state"]): r for r in rd(os.path.join(PT, "C2_ss2_state_tests_unified.tsv"))}
la = {(r["sig"], r["state"]): r for r in rd(os.path.join(PT, "C2_couturier_layerA_state_tests_all.tsv"))}
qc2 = rd(os.path.join(PT, "C2_4_QC_NEFTEL_SS2.tsv"))
qc10 = rd(os.path.join(PT, "C2_4_QC_NEFTEL_10X.tsv"))
roles = rd(os.path.join(PT, "C2_4_COUTURIER_SAMPLE_ROLES.tsv"))
ccdep = {r["signature"]: r for r in rd(os.path.join(PT, "C2_3_CC_DEPLETION_SUMMARY.tsv"))}
lb = rd(os.path.join(PT, "C2_LAYERB_STATE_TESTS.tsv"))
lbcnv = rd(os.path.join(PT, "C2_LAYERB_CNV_QC.tsv"))
c5x = rd(os.path.join(PT, "C3_C5_EXPLORATORY.tsv"))
neg = {r["negative_id"]: r for r in rd(os.path.join(PT, "C4_NEGATIVE_RESULTS.tsv"))}

STATES = ["MES", "AC", "OPC", "NPC", "Cycling"]
order = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
         "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN",
         "M_autophagy_lysosome", "C5_UP", "C5_STABLE_UP", "C5_DOWN", "C5_STABLE_DOWN",
         "M_HSF1_heat_shock"]
grade = {r["signature"]: r["final_evidence_grade"] for r in ev}

fig = plt.figure(figsize=(7.2, 11.2))
gs = fig.add_gridspec(4, 1, height_ratios=[1.2, 1.9, 0.95, 2.5], hspace=0.58)

# ---- A: layers ----
axA = fig.add_subplot(gs[0]); panel_label(axA, "A")
axA.set_axis_off()
ad = [r for r in qc2 if r.get("age_group") == "Adult"]
ss2_n = len(ad); ss2_mal = sum(int(r["malignant_cells"]) for r in ad)
lb_cells = sum(int(r["cnv_malignant_n"]) for r in lbcnv)
# Layer A values as frozen in the FINAL Fig 4 legend and C2_3_10 report (12 samples / 11 patients /
# 18,252 malignant cells; state tests on 8,505 core-state-annotated cells).
LA_SAMPLES, LA_PATIENTS, LA_MALIG, LA_CORE = 12, 11, 18252, 8505
txt = [
 "Layer A (Discovery): Neftel SS2 adult GBM, %d patients / %s malignant cells / 5 malignant states (8 pediatric patients, sensitivity only)" % (ss2_n, format(ss2_mal, ",")),
 "Layer A (Replication): Couturier GBM tissue, %d samples / %d independent patients / %s malignant cells; state tests on %s core-state-annotated cells" % (LA_SAMPLES, LA_PATIENTS, format(LA_MALIG, ","), format(LA_CORE, ",")),
 "Layer B (CNV-defined candidate cells): %d patients; infercnvpy; per-state patient-level test required by preregistration" % (len(set(r["patient"] for r in lbcnv))),
 "Layer B per-state testability: " + "; ".join(f"{r['state']} {r['cells']} cells/{r['patients']} patients vs min {r['n_cell_min_required']}/{r['n_pat_min_required']}" for r in lb) + " -> NOT TESTABLE",
]
_fsA, _nA = fit_rows_in_ax(axA, txt, fs_max=6.7, fs_min=5.4, top=0.96, lead=1.28, row_gap=0.45)
print("Fig4 panel A fontsize/rows: %.1f / %d" % (_fsA, _nA))
axA.set_title("Two independent single-cell layers and the preregistered Layer B coverage gate", fontsize=8)

# ---- B: state-localisation heatmaps ----
def make_mat(tab, cols=STATES):
    M = np.full((len(order), len(cols)), np.nan); S = np.full((len(order), len(cols)), False)
    for i, s in enumerate(order):
        for j, st in enumerate(cols):
            r = tab.get((s, st))
            if r and num(r["delta"]) is not None:
                M[i, j] = num(r["delta"]); S[i, j] = (num(r["fdr"]) or 1) < 0.05
    return M, S

gsL = gs[1].subgridspec(1, 2, wspace=0.08)
for k, (tab, ttl) in enumerate([(disc, "Discovery (Neftel SS2, 20 patients)"),
                                (la, "Layer A (Couturier, 11 patients)")]):
    ax = fig.add_subplot(gsL[k])
    if k == 0:
        panel_label(ax, "B", dx=-0.30, dy=1.02)
    M, S = make_mat(tab)
    im = ax.imshow(np.nan_to_num(M), cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.02, vcenter=0, vmax=0.02), aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]*1000:.1f}" + ("*" if S[i, j] else ""), ha="center", va="center",
                        fontsize=5.6, color="white" if abs(M[i, j]) > 0.014 else "black")
    ax.set_xticks(range(len(STATES))); ax.set_xticklabels(STATES, fontsize=6.2, rotation=45, ha="right")
    if k == 0:
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([f"{s}  [grade {grade.get(s,'-')}]" for s in order], fontsize=6.2)
    else:
        ax.set_yticks(range(len(order))); ax.set_yticklabels([])
    ax.set_title(ttl, fontsize=7.5)
    ax.set_aspect("auto")
cb = fig.colorbar(im, ax=[fig.axes[-2], fig.axes[-1]], shrink=0.7, pad=0.02,
                  ticks=[-0.02, -0.01, 0, 0.01, 0.02])
cb.ax.set_yticklabels(["-20", "-10", "0", "+10", "+20"], fontsize=6)
cb.set_label("AUC difference (x 10^-3)", fontsize=6.5)
# (the panel-B technical note is placed in the figure caption band at the
#  bottom of the canvas, see add_caption below - keeps it inside the canvas)

# ---- C: C5 diagnostic ----
axC = fig.add_subplot(gs[2]); panel_label(axC, "C")
c5sigs = ["C5_UP", "C5_STABLE_UP", "C5_DOWN", "C5_STABLE_DOWN"]
e5 = {r["signature"]: r for r in ev}
axes_ = [("discovery_delta", "discovery_fdr", "Discovery base", "#1f618d"),
         ("layerA_delta", "layerA_fdr", "Layer A base", "#c0392b"),
         ("ccdep_la_delta", "ccdep_la_fdr", "Layer A CC-depleted", "#7fb3d5")]
x = np.arange(len(c5sigs)); w = 0.26
for k, (dc, fc, lab, col) in enumerate(axes_):
    vals = [num(e5[s][dc]) * 1000 if num(e5[s][dc]) is not None else 0 for s in c5sigs]
    b = axC.bar(x + (k - 1) * w, vals, w, color=col, label=lab)
    for xi, s, v in zip(x + (k - 1) * w, c5sigs, vals):
        f = num(e5[s][fc])
        if f is not None:
            axC.text(xi, v + (0.15 if v >= 0 else -0.45), star(f), ha="center", fontsize=7)
axC.set_xticks(x); axC.set_xticklabels([f"{s}\n(CC removed {ccdep[s]['removed_fraction_pct']}%)" for s in c5sigs], fontsize=6.2)
axC.set_ylabel("mean patient-level\nAUC difference (x 10^-3)", fontsize=7)
axC.axhline(0, color="black", lw=0.7)
axC.legend(frameon=False, fontsize=6, ncol=3, loc="upper left")
axC.set_title("C5 is directionally concordant but incompletely replicated (Grade C; base Layer A MES n.s., FDR 0.0563)")
NOTE_C = ("Layer B and CGGA WHO-IV subtype contrasts are NOT TESTABLE (no molecular subtype annotation; insufficient cells/patients). "
          "TCGA-Verhaak (n=44): C5_UP MES vs non-MES median difference +476.2, FDR %.1e (Grade C, exploratory)."
          % num([r for r in c5x if r["signature"] == "C5_UP" and r["analysis"] == "C5_exploratory_MES_context"][0]["FDR_BH"]))

# ---- D: retained negatives ----
axD = fig.add_subplot(gs[3]); panel_label(axD, "D")
axD.set_axis_off()
hsf = [("Discovery / MES", disc.get(("M_HSF1_heat_shock", "MES"))), ("Layer A / MES", la.get(("M_HSF1_heat_shock", "MES"))),
       ("Layer A / NPC", la.get(("M_HSF1_heat_shock", "NPC")))]
lines = ["Retained negative / inconsistent findings (not removed from the figures):"]
for lab, r in hsf:
    if r:
        lines.append(f"  HSF1 heat-shock {lab}: delta {num(r['delta'])*1000:+.1f}e-3, FDR {num(r['fdr']):.3g}"
                     + ("  (reversed direction)" if lab.endswith("MES") and "Layer" in lab else ""))
lines.append("  Grade E: HSF1 heat-shock module is cohort-dependent and not replicated; it is retained as a heterogeneity finding.")
en = {
 "N6": "Layer B (CNV-defined cells) fails the preregistered coverage gate: per-state patient-level tests NOT TESTABLE.",
 "N7": "C5 base Layer A MES contrast not significant (FDR 0.0563); reported as directionally concordant, not replicated.",
 "N9": "C4_DOWN expanded-state mapping is attenuated (FDR 0.594); the Cycling anchor is retained with this sensitivity.",
 "N10": "Shared anti-proliferative expanded Cycling mapping is not significant (FDR 0.249); primary anchor mapping unchanged.",
 "N3": "Autophagy-lysosome (MES) is not significant in human bulk (TCGA n=44 MES P = 0.334); kept as Grade-B sensitivity only.",
}
for nid in ["N6", "N7", "N9", "N10", "N3"]:
    lines.append(f"  {nid}: " + en[nid])
lines.append("  Layer B: no figure element is derived from Layer B; coverage insufficiency is shown explicitly.")
_fsD, _nD = fit_rows_in_ax(axD, lines, fs_max=6.3, fs_min=5.4, top=0.98, lead=1.28, row_gap=0.40)
print("Fig4 panel D fontsize/rows: %.1f / %d" % (_fsD, _nD))
axD.set_title("Kept-in negatives and sensitivity caveats inherited from the frozen layer-B mapping", fontsize=8)

from fig_common import add_caption
add_caption(fig, ["Panel B: mean patient-level AUC difference; * FDR < 0.05 by per-patient paired Wilcoxon test. "
                  "RNA velocity / CNV panels are absent: Layer B is not testable.",
                  "Panel C: " + NOTE_C], fontsize=6.0, width=150)

save_fig(fig, "Fig4")
log_values("Fig4", [("A", "Layer A discovery patients/malignant cells", "C2_4_QC_NEFTEL_SS2.tsv (Adult rows)",
                     f"{ss2_n}|{ss2_mal}"),
                    ("A", "Layer A replication samples/patients/malignant/core-state", "FINAL Fig4 legend / C2_3_10 report",
                     f"{LA_SAMPLES}|{LA_PATIENTS}|{LA_MALIG}|{LA_CORE}"),
                    ("A", "Layer B CNV pool cells", "C2_LAYERB_CNV_QC.tsv", lb_cells)]
                   + [("B", "mean patient-level AUC delta", "C2_ss2_state_tests_unified.tsv / C2_couturier_layerA_state_tests_all.tsv",
                       f"{s}|Discovery={disc.get((s,st),{}).get('delta','NA')};LayerA={la.get((s,st),{}).get('delta','NA')}")
                      for s in order for st in STATES])
print("Fig4 done")
