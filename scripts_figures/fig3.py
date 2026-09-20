# -*- coding: utf-8 -*-
"""Fig3 - cross-platform / cross-model architecture. Frozen tables only, no recomputation."""
import os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from fig_common import *

REACT = "Contextualisation of HITMAN architecture across models"
meta = rd(os.path.join(PT, "B1_EMTAB9043_METADATA.tsv"))
focus = rd(os.path.join(TP, "C1_focus_pathways.tsv"))
xplat = rd(os.path.join(PT, "C1_C4_CROSS_PLATFORM_COMPARISON.tsv"))
c5x = rd(os.path.join(PT, "C1_C5_CROSS_PLATFORM_COMPARISON.tsv"))
neg = {r["negative_id"]: r for r in rd(os.path.join(PT, "C4_NEGATIVE_RESULTS.tsv"))}

STIM = {"none": "untreated", "Tumor Treating Fields (200KHz)": "TTFields 200 kHz",
        "Deep Brain Stimulation electric field (10V)": "10-V electric field"}
lines = ["GIN-28", "KNS-42", "U87-MG"]
ct = collections.Counter((r["cell_line_factor"], r["stimulus_factor"]) for r in meta)

fig = plt.figure(figsize=(7.2, 8.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.85], hspace=0.42)

# ---- A: comparator design ----
axA = fig.add_subplot(gs[0]); panel_label(axA, "A")
xs = np.arange(len(lines)); w = 0.26
colors = {"none": "#9e9e9e", "Tumor Treating Fields (200KHz)": "#1f618d",
          "Deep Brain Stimulation electric field (10V)": "#c0392b"}
for k, (key, lab) in enumerate(STIM.items()):
    vals = [ct[(l, key)] for l in lines]
    b = axA.bar(xs + (k - 1) * w, vals, w, label=lab, color=colors[key])
    axA.bar_label(b, fontsize=7)
axA.set_xticks(xs); axA.set_xticklabels(lines)
axA.set_ylabel("microarray samples (n)")
axA.set_ylim(0, 7)
axA.set_title("E-MTAB-9043: n = 40 arrays, 3 GBM lines, 6 field comparators (no matrix merging)")
axA.legend(frameon=False, loc="upper right")

# ---- B: pathway-level signed NES ----
axB = fig.add_subplot(gs[1]); panel_label(axB, "B", dx=-0.30, dy=1.02)
cols = ["GSE_C4", "GSE_C5", "EMTAB_U87MG_TTFields", "EMTAB_U87MG_DBS",
        "EMTAB_KNS42_TTFields", "EMTAB_KNS42_DBS", "EMTAB_GIN28_TTFields", "EMTAB_GIN28_DBS"]
clab = ["HITMAN\nC4", "HITMAN\nC5", "U87-MG\nTTFields", "U87-MG\n10 V",
        "KNS-42\nTTFields", "KNS-42\n10 V", "GIN-28\nTTFields", "GIN-28\n10 V"]
frows = [r for r in focus if r["db"] == "HALLMARK"]
frows.sort(key=lambda r: (num(r["GSE_C4"]) or 0))
M = np.array([[num(r[c]) or 0.0 for c in cols] for r in frows])
S = np.array([[1.0 if str(r[c]).endswith("*") else 0.0 for c in cols] for r in frows])
im = axB.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3), aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        mark = "*" if S[i, j] else ""
        axB.text(j, i, f"{M[i, j]:.2f}{mark}", ha="center", va="center", fontsize=5.8,
                 color="white" if abs(M[i, j]) > 1.9 else "black")
axB.set_xticks(range(len(cols))); axB.set_xticklabels(clab, fontsize=6.4)
axB.set_yticks(range(len(frows))); axB.set_yticklabels([r["pathway"] for r in frows], fontsize=6.4)
axB.set_title("Pathway-level signed NES (module level; * = FDR/adj. P < 0.05 as frozen in C1)")
fig.colorbar(im, ax=axB, shrink=0.55, label="signed NES")
t = []
t.append("Shared anti-proliferative / ER-UPR programmes are more consistently down across conventional")
t.append("electric-field comparators than across the HITMAN fields (module-level descriptive summary).")
t.append("U87-MG 10-V column is the outlier comparator that broadly activates most pathways (N11).")
gc = ", ".join(f"{r['cell_line']}-{r['comparator']} rho={float(r['spearman_rho']):.2f}"
               for r in xplat[:1] + xplat[3:4])
t.append(f"Gene-level signed concordance with comparators is absent ({gc}; N15).")
t.append(f"C5 cross-platform: {c5x[0]['cell_line']} {c5x[0]['comparator']} rho={float(c5x[0]['spearman_rho']):.3f} (n.s.); "
         f"C5 stable set dir_concord={c5x[0]['stable_dir']} (N12).")
t.append("C5 is therefore reported as a HITMAN-field-specific interaction programme, not a field-general signature.")
t.append("(A) U87-MG untreated n = 2 (comparator-outlier caveat, N11); SREBP/cholesterol and interferon "
         "programmes are the HITMAN-enriched sets; C5 is not reproduced on conventional electric-field platforms (N12).")
add_caption(fig, t, fontsize=6.2, width=140)

save_fig(fig, "Fig3")
log_values("Fig3", [("A", "arrays per cell line x stimulus", "B1_EMTAB9043_METADATA.tsv", f"{k}={v}") for k, v in sorted(ct.items())]
          + [("B", "signed NES", "C1_focus_pathways.tsv", f"{r['pathway']}|{c}={r[c]}") for r in frows for c in cols])
print("Fig3 done")
