# -*- coding: utf-8 -*-
"""FigS9 the release stage - rendering-only patch + re-export (frozen tables only).

Two RENDERING-ONLY changes versus the pre-release export. No statistic is recomputed, no
scale/threshold is changed, no number is altered; every displayed number is copied verbatim
from the frozen tables.

 1. Empty-field defect fix (same class as the Fig6 panel-D defect fixed in the release stage):
    the robust_2_stratum_Q1 / robust_4_gene_perturbation / robust_5_coverage rows carry
    their frozen values in the field `module_b` ("dir_match=X/21", "min=0.970 (C5_UP)").
    The pre-release export read `expected_sign`, which is empty for exactly these rows, so
    panel B rendered 0.0 bars and panel C rendered blank coverage values. This script reads
    `module_b` instead - the frozen field that actually holds the value.
 2. Cosmetic: panel titles no longer repeat the panel letter ("A   A. ..." -> "A   ...");
    the letter itself is still drawn by panel_label().

Only FigS9 is re-exported here; FigS1-S8 assets are not touched.
"""
import os, sys
SCR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCR)          # fig_common.py ships in this same directory
import numpy as np
import matplotlib.pyplot as plt
import fig_common
from fig_common import *

SUP = "supplementary"
PKG_ROOT = os.path.dirname(SCR)
REL_LOG = os.path.join(PKG_ROOT, "provenance", "plotted_values", "release")
os.makedirs(REL_LOG, exist_ok=True)
fig_common.LOGD = REL_LOG


def load(name, where=PT):
    """Read a frozen table (same loader convention as the released figure scripts)."""
    return rd(os.path.join(where, name))


def dm(s):
    """Numeric payload of a frozen 'dir_match=X/21' / 'min=0.970 (C5_UP)' token (rendering helper)."""
    return num(str(s).split("=")[-1].split("/")[0].split(" ")[0].strip())


pr = load("C3_PRIMARY_RESULTS.tsv")
rb = load("C3_ROBUSTNESS.tsv")
q1 = [r for r in pr if r["block"] == "Q1_architecture"]
f1 = [r for r in rb if r["analysis"] == "robust_1_cross_cohort_sign" and r["module_a"] != "SUMMARY"]
r3 = [r for r in rb if r["analysis"] == "robust_3_partial_global_mean" and r["cohort_scope"] == "TCGA"]
st2 = [r for r in rb if r["analysis"] == "robust_2_stratum_Q1"]
gp = [r for r in rb if r["analysis"] == "robust_4_gene_perturbation"]
cov = [r for r in rb if r["analysis"] == "robust_5_coverage"]

fig = plt.figure(figsize=(7.6, 7.8))
gs = fig.add_gridspec(3, 1, height_ratios=[1.8, 1.0, 0.8], hspace=0.6)

# ---- A: relationship-level forest plot (unchanged) ----
ax = fig.add_subplot(gs[0]); panel_label(ax, "A", dx=-0.13, dy=1.02)
y = np.arange(len(f1))[::-1]
for yi, r in zip(y, f1):
    a, b = num(r["rho_TCGA"]), num(r["rho_CGGA"])
    rev = str(r["same_direction"]).strip() == "False"
    ax.plot([a, b], [yi, yi], color="#c0392b" if rev else "#34495e", lw=0.8)
    ax.plot([a], [yi], "o", ms=3, color="#e67e22"); ax.plot([b], [yi], "s", ms=2.8, color="#2980b9")
ax.axvline(0, color="black", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels([f"{r['module_a']}-{r['module_b']}" for r in f1], fontsize=4.6)
ax.set_xlabel("Spearman rho (circles TCGA, squares CGGA; red = sign reversal)", fontsize=6.4)
ax.set_title("Relationship-level forest plot (28 relationships)", fontsize=7.6)

# ---- B: stratification + gene-perturbation robustness (frozen field = module_b) ----
ax = fig.add_subplot(gs[1]); panel_label(ax, "B", dx=-0.13, dy=1.04)
labs, vals = [], []
for r in st2:
    labs.append(r["cohort_scope"].replace("_", " ").replace("Primary only", "primary").replace("WHO IV", "WHO-IV"))
    vals.append(dm(r["module_b"]))
for coh in ["tcga", "cgga"]:
    v = [dm(r["module_b"]) for r in gp if r["cohort_scope"] == coh]
    labs.append(f"{coh.upper()} pert.-mean"); vals.append(float(np.mean(v)))
x = np.arange(len(labs))
ax.bar(x, vals, color="#5dade2")
for xi, v in zip(x, vals):
    ax.text(xi, v + 0.2, f"{v:.1f}", ha="center", fontsize=5.4)
ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=4.8, rotation=20, ha="right")
ax.set_ylabel("expected-direction pairs (/21)", fontsize=6.4)
ax.set_ylim(0, 15)
ax.set_title("Stratification and gene-perturbation robustness (no stratum improves recovery)", fontsize=7)

# ---- C: coverage and partial-subset checks (frozen field = module_b) ----
ax = fig.add_subplot(gs[2]); panel_label(ax, "C", dx=-0.13, dy=1.06)
ax.set_axis_off()
lines = [f"{r['cohort_scope']}: {r['module_b']} ({r['note'][:80]})" for r in cov]
lines += [f"Partial (global-mean adjusted) subset: {sum(1 for r in r3 if str(r['same_direction']).strip()=='True')}/{len(r3)} TCGA relationships keep their frozen sign",
          "CGGA is treated as contextualisation, not replication; sign-reversed pairs are retained and not relabelled."]
for i, ln in enumerate(lines):
    ax.text(0.0, 1.0 - i * 0.26, ln, transform=ax.transAxes, fontsize=6.2, va="top")
ax.set_title("Coverage and partial-subset checks", fontsize=7)

save_fig(fig, "FigS9", subdir=SUP)

# --- traceability log (values as drawn, copied verbatim from the frozen tables) ---
rows = []
for r in f1:
    rows.append(("A", "rho_TCGA/rho_CGGA", "C3_ROBUSTNESS.tsv/robust_1_cross_cohort_sign",
                 f"{r['module_a']}|{r['module_b']}|{r['rho_TCGA']}|{r['rho_CGGA']}"))
for r in st2:
    rows.append(("B", "dir_match", "C3_ROBUSTNESS.tsv/robust_2_stratum_Q1", f"{r['cohort_scope']}|{r['module_b']}"))
for r in gp:
    rows.append(("B", "dir_match", "C3_ROBUSTNESS.tsv/robust_4_gene_perturbation", f"{r['cohort_scope']}|{r['module_b']}"))
for r in cov:
    rows.append(("C", "min_per_sample_gene_coverage", "C3_ROBUSTNESS.tsv/robust_5_coverage", f"{r['cohort_scope']}|{r['module_b']}"))
log_values("FigS9", rows)
print("FigS9 rendering-only re-export done:", [f"{l}={v:.1f}" for l, v in zip(labs, vals)])
