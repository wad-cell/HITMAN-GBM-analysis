# -*- coding: utf-8 -*-
"""Fig2 (FINAL-LEGEND COMPLIANT: panels A-C) - differential response.

Layout-only variant of fig2.py: the FINAL legend defines Fig 2 as A-C, and the
manuscript body never cites a panel 2D, so the evidence-grade bar panel is not
rendered here (the superseded A-D export is retained outside this package).  Panel C spans the
bottom row so the figure has no empty slot.  Frozen tables only; no recomputation.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig_common import PT, TP, C3D, rd, num, panel_label, save_fig, add_caption

fp = rd(os.path.join(TP, "C1_focus_pathways.tsv"))
arms = ["GSE_C4", "GSE_C5", "EMTAB_U87MG_TTFields", "EMTAB_KNS42_TTFields", "EMTAB_GIN28_TTFields",
        "EMTAB_U87MG_DBS", "EMTAB_KNS42_DBS", "EMTAB_GIN28_DBS"]

shared_down = ["HALLMARK\tE2F Targets", "HALLMARK\tMyc Targets V1",
               "HALLMARK\tUnfolded Protein Response", "REACTOME\tCell Cycle Checkpoints R-HSA-69620",
               "REACTOME\tApoptosis R-HSA-109581"]
enriched_up = ["HALLMARK\tCholesterol Homeostasis", "REACTOME\tCholesterol Biosynthesis R-HSA-191273",
               "REACTOME\tActivation Of Gene Expression By SREBF (SREBP) R-HSA-2426168",
               "HALLMARK\tTNF-alpha Signaling via NF-kB", "HALLMARK\tp53 Pathway",
               "REACTOME\tMacroautophagy R-HSA-1632852"]
idx = {(r["db"], r["pathway"]): r for r in fp}

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, wspace=0.30, hspace=0.45, left=0.06, right=0.97, top=0.90, bottom=0.09)


def dotplot(ax, keys, title):
    arm_styles = [("GSE_C4", "#b2182b", "o", 46), ("GSE_C5", "#8c510a", "D", 40),
                  ("TT", "#2f4b7c", "^", 34), ("DBS", "#7f9dc0", "s", 30)]
    for j, key in enumerate(keys):
        r = idx[tuple(key.split("\t"))]
        y = len(keys) - 1 - j
        for k, (tag, c, mk, s) in enumerate(arm_styles):
            if tag in ("GSE_C4", "GSE_C5"):
                vals = [num(r[tag])]
            else:
                vals = [num(r[a]) for a in arms if tag in a]
            for v in vals:
                ax.scatter(v, y + (k - 1.5) * 0.15, c=c, marker=mk, s=s, zorder=3,
                           edgecolor="white", linewidth=0.4)
        from fig_common import wrap_text
        lab = key.split("\t")[1]
        _p = ax.get_position()
        _x0 = (_p.x0 + 0.02 * _p.width) * fig.get_size_inches()[0]
        _avail = (_p.x1 - 0.03 * _p.width) * fig.get_size_inches()[0] - _x0
        _fs = 7.0
        _parts = wrap_text(lab, _avail, _fs)
        if len(_parts) > 2:
            _fs = 6.2
            _parts = wrap_text(lab, _avail, _fs)
        if len(_parts) > 3:
            _fs = 5.8
            _parts = wrap_text(lab, _avail, _fs)
        ax.text(0.02, y, "\n".join(_parts), transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=_fs,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6))
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_yticks([]); ax.set_ylim(-0.7, len(keys) - 0.3)
    ax.set_xlabel("Signed NES (frozen GSEA/ssGSEA tables)")
    ax.set_title(title, fontsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)


# ---- A shared down-regulated programs ----
axA = fig.add_subplot(gs[0, 0])
dotplot(axA, shared_down, "Shared down-regulated programs (HITMAN arms vs comparators)")
cvals = [num(idx[tuple(k.split("\t"))][a]) for k in shared_down
         for a in arms if "TTFields" in a or "DBS" in a]
axA.text(0.99, 0.02, "comparator signed NES range\n%.2f to %.2f" % (min(cvals), max(cvals)),
         transform=axA.transAxes, ha="right", va="bottom", fontsize=6.6, color="#2f4b7c")
panel_label(axA, "A")

# ---- B HITMAN-enriched programs ----
axB = fig.add_subplot(gs[0, 1])
dotplot(axB, enriched_up, "HITMAN-enriched programs (SREBP/cholesterol, TNF/NF-kB, p53, autophagy)")
panel_label(axB, "B")
axB.legend(handles=[Line2D([], [], marker=m, ls="", color=c, label=l, markersize=5)
                    for l, c, m in [("GSE C4", "#b2182b", "o"), ("GSE C5", "#8c510a", "D"),
                                    ("E-MTAB-9043 TTFields", "#2f4b7c", "^"),
                                    ("E-MTAB-9043 10 V (DBS)", "#7f9dc0", "s")]],
           loc="lower right", frameon=False, fontsize=5.6, handlelength=1.1,
           handletextpad=0.45, labelspacing=0.32, borderpad=0.2)

# ---- C C5 exploratory contextualization (spans the bottom row) ----
axC = fig.add_subplot(gs[1, :])
_expall = rd(os.path.join(PT, "C3_C5_EXPLORATORY.tsv"))
exp = [r for r in _expall if r["analysis"] == "C5_exploratory_MES_context"]
labs, eff, lo, hi, stars = [], [], [], [], []
for r in exp:
    labs.append(r["signature"].replace("C5_", "").replace("_", " "))
    eff.append(num(r["effect"])); lo.append(num(r["CI95_low"])); hi.append(num(r["CI95_high"]))
    p = num(r["P"])
    stars.append("***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s.")
y = np.arange(len(labs))[::-1]
axC.errorbar(eff, y, xerr=[np.array(eff) - np.array(lo), np.array(hi) - np.array(eff)],
             fmt="o", color="#2f4b7c", ecolor="#7f9dc0", capsize=2.5, ms=5)
for i, (yy, s) in enumerate(zip(y, stars)):
    axC.text(hi[i] + 30, yy, s, fontsize=7, va="center")
axC.axvline(0, color="#888", lw=0.8, ls="--")
axC.set_yticks(y); axC.set_yticklabels(labs, fontsize=7)
axC.set_xlabel("MES - non-MES median score difference (TCGA n = 44)")
axC.set_title("C5 exploratory overlay in TCGA MES context (Grade C)", fontsize=8)
axC.spines[["top", "right"]].set_visible(False)
axC.set_xlim(-200, max(hi) * 1.30)
axC.text(0.99, 0.03, "C5-up axis correlation with frozen C4-up axis:\n"
                     "TCGA rho = 0.56, CGGA rho = 0.96 (frozen table)",
         transform=axC.transAxes, ha="right", va="bottom", fontsize=6.4, color="#555")
panel_label(axC, "C", dx=-0.10)

add_caption(fig, ["Evidence grades are unchanged and carried over verbatim from the frozen table (C5 = C, HSF1 = E); "
                  "no grade was upgraded."],
            fontsize=6.4, width=150)

save_fig(fig, "Fig2")
print("Fig2 (A-C, FINAL legend compliant) done")
