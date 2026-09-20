"""Fig1 - Study design, estimands, and the frozen module architecture.
Frozen sources: B2_FROZEN_SAMPLE_TABLE.tsv, B2_CONTRAST_DICTIONARY.tsv,
B35_DESEQ2_SUMMARY.tsv, C2_2_FROZEN_SIGNATURES_COVERAGE.tsv,
C2_12_MASTER_EVIDENCE.tsv, C2_LAYERB_STATE_TESTS.tsv
"""
import os, collections
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from fig_common import PT, rd, num, panel_label, save_fig, log_values, add_caption

fc = rd(os.path.join(PT, "B2_FROZEN_SAMPLE_TABLE.tsv"))
cont = rd(os.path.join(PT, "B2_CONTRAST_DICTIONARY.tsv"))
b35 = rd(os.path.join(PT, "B35_DESEQ2_SUMMARY.tsv"))
cov = rd(os.path.join(PT, "C2_2_FROZEN_SIGNATURES_COVERAGE.tsv"))
lb = rd(os.path.join(PT, "C2_LAYERB_STATE_TESTS.tsv"))

b1 = rd(os.path.join(PT, "B1_MASTER_SAMPLE_SHEET.tsv"))
d_cnt = collections.Counter(r["group_short"] for r in b1)
m_cnt = collections.Counter(r["group_short"] for r in fc if r.get("in_matrix", "Y") == "Y")
n_design = {"MAN-F": d_cnt.get("MAN-F", 0), "ND-F": d_cnt.get("ND-F", 0),
            "MAN-NF": d_cnt.get("MAN-NF", 0), "ND-NF": d_cnt.get("ND-NF", 0)}
n_matrix = {"MAN-F": m_cnt.get("MAN-F", 0), "ND-F": m_cnt.get("ND-F", 0),
            "MAN-NF": m_cnt.get("MAN-NF", 0), "ND-NF": m_cnt.get("ND-NF", 0)}
n_total = sum(n_design.values())
n_matrix_total = sum(n_matrix.values())
sens_excl = [r["matrix_col"] for r in fc if r.get("sensitivity_exclude", "").strip()]

d = {r["contrast"]: r for r in b35}
c4_r, c5_r = num(d["C4"]["r_n_FDR005"]), num(d["C5"]["r_n_FDR005"])
c4_py, c5_py = num(d["C4"]["py_n_FDR005"]), num(d["C5"]["py_n_FDR005"])
stable_n = sum(num(r["original_n_ENSG"]) for r in cov if r["signature"].startswith("C5_STABLE"))

fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, wspace=0.28, hspace=0.42, left=0.09, right=0.97, top=0.93, bottom=0.08)

# ---------------- Panel A: factorial design ----------------
axA = fig.add_subplot(gs[0, 0])
axA.set_xlim(0, 2); axA.set_ylim(0, 2); axA.set_axis_off()
cells = {(1, 1): ("MAN-F", n_design["MAN-F"], n_matrix["MAN-F"]),
         (0, 1): ("ND-F", n_design["ND-F"], n_matrix["ND-F"]),
         (1, 0): ("MAN-NF", n_design["MAN-NF"], n_matrix["MAN-NF"]),
         (0, 0): ("ND-NF", n_design["ND-NF"], n_matrix["ND-NF"])}
for (x, y), (lab, nd, nm) in cells.items():
    axA.add_patch(FancyBboxPatch((x + .06, y + .12), .80, .70, boxstyle="round,pad=0.02",
                                 fc="#eef3fa", ec="#2f4b7c", lw=1.0))
    axA.text(x + .46, y + .66, lab, ha="center", va="center", fontsize=8.5, fontweight="bold")
    axA.text(x + .46, y + .45, "n = %d" % nd, ha="center", va="center", fontsize=8)
    axA.text(x + .46, y + .26, "%d in matrix" % nm, ha="center", va="center", fontsize=6.2, color="#444")
axA.text(0.06, 1.88, "Field (magnetic field)  +", fontsize=6.6, color="#444")
axA.text(1.06, 1.88, "Field  -", fontsize=6.6, color="#444")
axA.text(0.02, 1.62, "MAN\n+", fontsize=7, color="#444", ha="center")
axA.text(0.02, 0.62, "MAN\n-", fontsize=7, color="#444", ha="center")
axA.set_title("Factorial design: HITMAN (MAN) x electric field", fontsize=8.5)
axA.text(1.0, 0.10, "2 x 2 design: %d arrays, GSE324656 (%d per arm)\n%d arrays entered the analysis matrix; excluded: %s"
         % (n_total, n_total // 4, n_matrix_total, ", ".join(sens_excl) if sens_excl else "none"),
         ha="center", va="top", fontsize=6.2, color="#333")
panel_label(axA, "A")

# ---------------- Panel B: estimands ----------------
axB = fig.add_subplot(gs[0, 1])
axB.set_axis_off()
lines = [("Estimand (frozen contrast dictionary)", True)]
for r in cont:
    prim = str(r.get("primary", "")).strip().lower() == "yes"
    txt = "%s: %s" % (r["contrast_id"], r["contrast_definition"])
    lines.append((txt + ("   [PRIMARY interaction]" if prim else ""), prim))
from fig_common import wrap_text
_wb = axB.get_position().width * fig.get_size_inches()[0] - 0.02
_hb = axB.get_position().height * fig.get_size_inches()[1]
_lhb = (7.4 * 1.30 / 72.0) / _hb
y = 0.99
for txt, prim in lines:
    fs = 7.8 if prim else 7.6
    parts = wrap_text(txt, _wb, fs)
    axB.text(0.0, y, "\n".join(parts), fontsize=fs, va="top", ha="left",
             fontweight="bold" if prim else "normal",
             color="#b2182b" if prim else "#111", linespacing=1.30)
    y -= _lhb * (len(parts) + 0.55)
axB.text(0.0, y - 0.02, "All panels/statistics in this figure set are derived from frozen\n"
                        "result tables; no re-analysis was performed.", fontsize=6.6, color="#555",
         va="top", ha="left")
panel_label(axB, "B")

# ---------------- Panel C: significant genes ----------------
axC = fig.add_subplot(gs[1, 0])
labels = ["C4\n(total HITMAN-\nassociated)", "C5\n(interaction)", "C5 stable\nsubset"]
vals = [c4_r, c5_r, stable_n]
bars = axC.bar(range(3), vals, color=["#2f4b7c", "#2f4b7c", "#7f9dc0"], width=0.6)
for i, v in enumerate(vals):
    axC.text(i, v + max(vals) * 0.025, "%s" % format(int(v), ","), ha="center", fontsize=8)
axC.set_xticks(range(3)); axC.set_xticklabels(labels, fontsize=6.8)
axC.set_ylabel("Genes at FDR < 0.05")
axC.set_ylim(0, max(vals) * 1.20)
axC.spines[["top", "right"]].set_visible(False)
axC.set_title("Differential response (R/DESeq2 engine)", fontsize=8)
axC.text(0.02, 0.98, "PyDESeq2 engine: C4 %s, C5 %s\n(same thresholds; engine counts in Supplementary Table S2)"
         % (format(int(c4_py), ","), format(int(c5_py), ",")),
         transform=axC.transAxes, fontsize=6.2, va="top", color="#555")
panel_label(axC, "C")

# ---------------- Panel D: module architecture ----------------
axD = fig.add_subplot(gs[1, 1])
axD.set_axis_off()
shared = [r["signature"] for r in cov if r["signature"].startswith("M_")]
enriched = ["C4_UP", "C4_DOWN", "C5_UP", "C5_DOWN", "C5_STABLE_UP (682)", "C5_STABLE_DOWN (920)"]
axD.text(0.0, 0.98, "Frozen module architecture", fontsize=8.5, fontweight="bold", va="top")
axD.text(0.0, 0.88, "Shared response modules (%d)" % len(shared), fontsize=7.8, fontweight="bold",
         va="top", color="#2f4b7c")
for i, s in enumerate(shared):
    axD.text(0.02, 0.80 - i * 0.058, "- " + s, fontsize=7, va="top")
axD.text(0.0, 0.80 - len(shared) * 0.058 - 0.04, "HITMAN-enriched components", fontsize=7.8,
         fontweight="bold", va="top", color="#b2182b")
for i, s in enumerate(enriched):
    axD.text(0.02, 0.80 - len(shared) * 0.058 - 0.12 - i * 0.058, "- " + s, fontsize=7, va="top")
panel_label(axD, "D")

add_caption(fig, ["Layer B per-patient state tests: NOT TESTABLE (%s)" % lb[0]["reason"]], fontsize=6.4, width=150)
paths = save_fig(fig, "Fig1")
log_values("Fig1", [
    ["A", "design arm counts (B1 master sheet)", "B1_MASTER_SAMPLE_SHEET.tsv", str(n_design)],
    ["A", "analysis-matrix arm counts", "B2_FROZEN_SAMPLE_TABLE.tsv", str(n_matrix)],
    ["C", "C4 genes FDR<0.05 (R/DESeq2)", "B35_DESEQ2_SUMMARY.tsv", c4_r],
    ["C", "C5 genes FDR<0.05 (R/DESeq2)", "B35_DESEQ2_SUMMARY.tsv", c5_r],
    ["C", "C5 stable subset genes", "C2_2_FROZEN_SIGNATURES_COVERAGE.tsv", stable_n],
    ["C", "C4/C5 genes FDR<0.05 (PyDESeq2)", "B35_DESEQ2_SUMMARY.tsv", "%s/%s" % (c4_py, c5_py)],
    ["D", "shared modules n", "C2_2_FROZEN_SIGNATURES_COVERAGE.tsv", len(shared)],
])
print("Fig1 OK", paths["tif"])
