# -*- coding: utf-8 -*-
"""Supplementary Figures S1-S9 (released figure set). Frozen tables only; no recomputation."""
import os, sys, textwrap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from fig_common import *

SUP = "supplementary"
STATES = ["MES", "AC", "OPC", "NPC", "Cycling"]

def load(name, where=PT):
    return rd(os.path.join(where, name))

def dmap(tab, key=None, pkey=None, extra=None):
    cols = list(tab[0].keys()) if tab else []
    if key is None:
        key = "delta" if "delta" in cols else "delta_mean"
    if pkey is None:
        pkey = "fdr" if "fdr" in cols else "FDR"
    d = {}
    for r in tab:
        s = r["sig"]
        if extra and extra not in s:
            continue
        d[(s, r["state"])] = (num(r[key]), num(r[pkey]))
    return d

# ============ S1 cohort QC ============
q_ss2 = load("C2_4_QC_NEFTEL_SS2.tsv")
q_10x = load("C2_4_QC_NEFTEL_10X.tsv")
q_cou = load("C2_4_COUTURIER_SAMPLE_ROLES.tsv")
fig = plt.figure(figsize=(7.4, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.1, 1.0], wspace=0.85)
adults = [r for r in q_ss2 if r["age_group"].lower().startswith("adult")]
ax = fig.add_subplot(gs[0]); panel_label(ax, "A", dx=-0.16, dy=1.05)
labs = [r["sample"] for r in q_ss2]
cells = [num(r["cells"]) for r in q_ss2]
mal = [100.0 * (num(r["malignant_cells"]) or 0) / (num(r["cells"]) or 1) for r in q_ss2]
col = ["#2874a6" if r["age_group"].lower().startswith("adult") else "#c0392b" for r in q_ss2]
ax.bar(range(len(labs)), cells, color=col)
ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=90, fontsize=4.6, color="#c0392b")
ax.set_ylabel("cells per sample", fontsize=8)
ax2 = ax.twinx(); ax2.plot(range(len(labs)), mal, "k.-", ms=3, lw=0.8)
ax2.set_ylabel("malignant fraction (%)", fontsize=8)
ax.set_title("Neftel Smart-seq2 samples (blue: adult; red: pediatric)", fontsize=7.5)
ax = fig.add_subplot(gs[1]); panel_label(ax, "B", dx=-0.24, dy=1.05)
roles = []
for r in q_cou:
    if r["role"] not in roles:
        roles.append(r["role"])
pal = plt.get_cmap("tab10")
for i, rl in enumerate(roles):
    rs = [r for r in q_cou if r["role"] == rl]
    fr = [100.0 * (num(r["malignant"]) or 0) / (num(r["cells"]) or 1) for r in rs]
    ax.scatter([num(r["cells"]) for r in rs], fr, s=16, color=pal(i % 10), label=rl)
ax.set_xscale("log"); ax.set_xlabel("cells per sample (log)", fontsize=7)
ax.set_ylabel("malignant fraction (%)", fontsize=7)
ax.legend(frameon=False, fontsize=5.4)
ax.set_title("Couturier samples by role", fontsize=7.5)
ax = fig.add_subplot(gs[2]); panel_label(ax, "C", dx=-0.30, dy=1.05)
ax.bar(range(len(q_10x)), [num(r["cells"]) for r in q_10x], color="#7dcea0")
ax.plot(range(len(q_10x)), [num(r["malignant"]) for r in q_10x], "k.-", ms=3, lw=0.8)
ax.set_xticks(range(len(q_10x))); ax.set_xticklabels([r["sample"] for r in q_10x], rotation=90, fontsize=4.6)
ax.set_ylabel("cells (bars) / malignant (line)", fontsize=6.6)
ax.set_title("Neftel 10x sensitivity", fontsize=7.5)
add_caption(fig, ["Discovery layer uses the %d adult samples (blue bars; red axis labels mark pediatric samples, which are "
                  "not part of the frozen inference set). Mean genes/cell are given in the frozen QC table." % len(adults)],
            fontsize=6.2, width=132)
save_fig(fig, "FigS1", subdir=SUP)
log_values("FigS1", [("AB", "cells/malignant", "C2_4_QC_*", f"{r['sample']}|{r['cells']}|{r['malignant_cells'] if 'malignant_cells' in r else r['malignant']}") for r in q_ss2 + q_cou])

# ============ S2 state localization detail (SREBP, TNF, p53) ============
disc = dmap(load("C2_ss2_state_tests_unified.tsv"))
la = dmap(load("C2_couturier_layerA_state_tests_all.tsv"))
SIGS = ["M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis"]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.6))
for k, (tab, ttl) in enumerate([(disc, "Discovery (20 patients)"), (la, "Layer A replication (11 patients)")]):
    ax = axes[k]; panel_label(ax, "AB"[k], dx=-0.14, dy=1.06)
    w = 0.26
    for i, s in enumerate(SIGS):
        v = [(tab.get((s, st), (None, None))[0] or 0) * 1000 for st in STATES]
        ax.bar(np.arange(len(STATES)) + (i - 1) * w, v, width=w, label=s.replace("M_", "").replace("_", " "))
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(range(len(STATES))); ax.set_xticklabels(STATES, fontsize=6.5)
    ax.set_ylabel("delta x 10^-3", fontsize=7); ax.set_title(ttl, fontsize=8)
    if k == 0:
        ax.legend(frameon=False, fontsize=5.6, ncol=1)
rows_txt = []
for s in SIGS:
    for coh, tab in [("Discovery", disc), ("LayerA", la)]:
        vals = [f"{st} d={tab.get((s, st), (None, None))[0]*1000:.2f}, FDR={tab.get((s, st), (None, None))[1]:.3g}"
                if tab.get((s, st)) else "" for st in STATES]
        rows_txt.append(f"{s} [{coh}] " + "; ".join(v for v in vals if v))
add_caption(fig, ["Per-state effect sizes and BH-FDR: " + rows_txt[0]] + rows_txt[1:], fontsize=5.4, width=170)
save_fig(fig, "FigS2", subdir=SUP)
log_values("FigS2", [("AB", "delta", "C2_ss2_state_tests_unified.tsv | C2_couturier_layerA_state_tests_all.tsv", t) for t in rows_txt])

# ============ S3 grade-B sensitivity (ER-UPR, autophagy) ============
d_cc = dmap(load("ss2_state_tests.tsv", TP))
la_cc = dmap(load("c2_11_test_A_ccdep.tsv", TP))
la_alt = dmap(load("C2_couturier_layerA_alt_state_tests.tsv"))
VAR = [("Discovery core", disc, "M_shared_ER_UPR", "M_autophagy_lysosome"),
       ("Discovery CCDEP", d_cc, "M_shared_ER_UPR_CCDEP", "M_autophagy_lysosome_CCDEP"),
       ("Layer A core", la, "M_shared_ER_UPR", "M_autophagy_lysosome"),
       ("Layer A CCDEP", la_cc, "M_shared_ER_UPR_CCDEP", "M_autophagy_lysosome_CCDEP"),
       ("Layer A log1p-mean", la_alt, "ALT_M_shared_ER_UPR", "ALT_M_autophagy_lysosome")]
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
for j, st in enumerate(["MES", "AC", "NPC"]):
    ax = axes[j]; panel_label(ax, "ABC"[j], dx=-0.16, dy=1.06)
    w = 0.38
    for i, (lab, tab, s1, s2) in enumerate(VAR):
        for k, s in enumerate([s1, s2]):
            v = (tab.get((s, st), (None, None))[0] or 0) * 1000
            ax.bar(i + (k - 0.5) * w, v, width=w, color=["#2980b9", "#e67e22"][k],
                   alpha=0.9 if i % 2 == 0 else 0.6)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(range(len(VAR))); ax.set_xticklabels(["D core", "D CCDEP", "A core", "A CCDEP", "A log1p"],
                                                      rotation=45, fontsize=5.6)
    ax.set_title(f"state = {st}", fontsize=8); ax.set_ylabel("delta x 10^-3", fontsize=7)
axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#2980b9"), plt.Rectangle((0, 0), 1, 1, color="#e67e22")],
               labels=["ER-UPR", "autophagy/lysosome"], frameon=False, fontsize=5.6)
cc = {r["signature"]: r for r in load("C2_3_CC_DEPLETION_SUMMARY.tsv")}
add_caption(fig, ["Cell-cycle-depleted gene removal (C2_3_CC_DEPLETION_SUMMARY.tsv): " +
                  "; ".join(f"{k}: {v['cellcycle_removed_n']}/{v['unique_symbols']} genes ({v['removed_fraction_pct']}%), "
                            f"remaining coverage {v['remaining_coverage_Couturier_pct']}% (Couturier)" for k, v in cc.items()),
                  "Grade-B components are reported as sensitivity-dependent only: under CCDEP the ER-UPR MES component "
                  "does not survive in Discovery; autophagy/lysosome is not significant in the bulk TCGA MES contrast "
                  "(P = 0.334). No grade upgrade is made."], fontsize=5.4, width=175)
save_fig(fig, "FigS3", subdir=SUP)
log_values("FigS3", [("ABC", "delta", "M_shared_ER_UPR/M_autophagy_lysosome variants", f"{lab}|{st}|{s}|{v}") for lab, tab, s1, s2 in VAR for s in (s1, s2) for st in STATES for v in [tab.get((s, st), (None, None))[0]]])

# ============ S4 grade-D / grade-E diagnostics ============
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.4), gridspec_kw={"width_ratios": [1, 1, 1.1]})
for j, (name, s_d, s_a) in enumerate([("RIG-I (grade D, Discovery-only)", "M_RIGI_typeI_IFN", "M_RIGI_typeI_IFN"),
                                      ("HSF1 (grade E, cohort reversal)", "M_HSF1_heat_shock", "M_HSF1_heat_shock")]):
    ax = axes[j]; panel_label(ax, "AB"[j], dx=-0.18, dy=1.06)
    vd = [(disc.get((s_d, st), (None, None))[0] or 0) * 1000 for st in STATES]
    va = [(la.get((s_a, st), (None, None))[0] or 0) * 1000 for st in STATES]
    x = np.arange(len(STATES))
    ax.bar(x - 0.19, vd, width=0.38, color="#e67e22", label="Discovery")
    ax.bar(x + 0.19, va, width=0.38, color="#2980b9", label="Layer A")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(x); ax.set_xticklabels(STATES, fontsize=6)
    ax.set_ylabel("delta x 10^-3", fontsize=7); ax.set_title(name, fontsize=7.5)
    ax.legend(frameon=False, fontsize=5.6)
ax = axes[2]; panel_label(ax, "C", dx=-0.20, dy=1.06)
lopo = load("ss2_lopo.tsv", TP)
flips = {r["sig"]: (num(r["loo_flips_top"]), r["flag_patient_sensitive"]) for r in lopo}
sel = ["AUC_M_RIGI_typeI_IFN", "AUC_M_HSF1_heat_shock", "AUC_M_RIGI_typeI_IFN_CCDEP", "AUC_M_HSF1_heat_shock_CCDEP"]
labs = [s.replace("AUC_", "") for s in sel]
vals = [flips.get(s, (0, ""))[0] or 0 for s in sel]
ax.barh(np.arange(len(sel))[::-1], vals, color=["#c0392b" if v else "#1e8449" for v in vals])
for yi, s, v in zip(np.arange(len(sel))[::-1], sel, vals):
    ax.text(v + 0.03, yi, f"{int(v)} ({flips.get(s, (0, ''))[1]})", fontsize=5.4, va="center")
ax.set_yticks(np.arange(len(sel))[::-1]); ax.set_yticklabels(labs, fontsize=5.4)
ax.set_xlabel("LOPO top-state changes (Discovery, 20 patients)", fontsize=6.4)
ax.set_xlim(0, 3)
ax.set_title("LOPO summaries", fontsize=7.5)
add_caption(fig, ["RIG-I: significant in Discovery MES but not replicated in Layer A; retained as grade D (Discovery-only). "
                  "HSF1: Discovery MES FDR 0.002 vs Layer A MES FDR 1.000 (reversal) -> grade E, cohort-dependent; no upgrade."],
            fontsize=5.8, width=165)
save_fig(fig, "FigS4", subdir=SUP)
log_values("FigS4", [("AB", "delta", "unified/LayerA state tests", ""), ("C", "loo_flips_top", "ss2_lopo.tsv", f"{s}={flips.get(s, (None,))[0]}")])

# ============ S5 alternative scoring (Layer A, log1p-mean) ============
fig = plt.figure(figsize=(7.4, 5.6))
gsS5 = fig.add_gridspec(1, 2, wspace=0.55, left=0.155, right=0.98, top=0.90, bottom=0.13)
pairs = []
for s in la_alt:
    base_sig = s[0].replace("ALT_", "")
    if (base_sig, s[1]) in la and la_alt[s][0] is not None and la[(base_sig, s[1])][0] is not None:
        pairs.append((f"{base_sig}:{s[1]}", la[(base_sig, s[1])][0] * 1000, la_alt[s][0] * 1000))
pairs.sort(key=lambda t: t[1])
HALF = int(np.ceil(len(pairs) / 2.0))
for k in range(2):
    ax = fig.add_subplot(gsS5[k])
    sub = pairs[k * HALF:(k + 1) * HALF]
    yk = np.arange(len(sub))
    ax.plot([t[1] for t in sub], yk, "o", ms=3.6, color="#e67e22", label="base (core-state AUC)")
    ax.plot([t[2] for t in sub], yk, "s", ms=3.2, color="#2980b9", label="log1p-mean (alternative)")
    for yi, t in zip(yk, sub):
        ax.plot([t[1], t[2]], [yi, yi], color="grey", lw=0.6)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(yk); ax.set_yticklabels([t[0] for t in sub], fontsize=4.6)
    ax.set_ylim(-0.8, len(sub) - 0.2)
    ax.set_xlabel("delta x 10^-3", fontsize=7)
    if k == 0:
        panel_label(ax, "A", dx=-0.24, dy=1.02)
        ax.legend(frameon=False, fontsize=6, loc="lower right")
        ax.set_title("Layer A alternative scoring sensitivity (log1p-mean) - flagged as confounded by expression\n"
                     "abundance; not used for primary inference (pairs 1-%d of %d)" % (len(sub), len(pairs)), fontsize=6.6)
    else:
        ax.set_title("same contrasts, pairs %d-%d of %d (continued)" % (k * HALF + 1, len(pairs), len(pairs)), fontsize=6.6)
add_caption(fig, ["Sign/effect agreement between base and alternative scoring is partial; the alternative metric is "
                  "reported for transparency only and no primary claim is changed."], fontsize=6.2, width=140)
save_fig(fig, "FigS5", subdir=SUP)
log_values("FigS5", [("A", "delta", "C2_couturier_layerA_alt_state_tests.tsv | C2_couturier_layerA_state_tests_all.tsv", t[0]) for t in pairs])

# ============ S6 Layer B testability ============
lb_state = load("C2_LAYERB_STATE_TESTS.tsv")
lb_qc = load("C2_LAYERB_CNV_QC.tsv")
fig = plt.figure(figsize=(7.4, 2.8))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.4, left=0.075, right=0.98, top=0.88, bottom=0.16)
ax = fig.add_subplot(gs[0]); panel_label(ax, "A", dx=-0.18, dy=1.05)
x = np.arange(len(lb_state))
ax.bar(x - 0.2, [num(r["patients"]) for r in lb_state], width=0.38, color="#2980b9", label="patients")
ax.bar(x + 0.2, [num(r["cells"]) for r in lb_state], width=0.38, color="#e67e22", label="cells")
ax.axhline(2, color="#2980b9", ls="--", lw=0.9)
ax.axhline(3, color="#e67e22", ls="--", lw=0.9)
ax.set_xticks(x); ax.set_xticklabels([r["state"] for r in lb_state], fontsize=6.2)
ax.set_ylabel("count", fontsize=7)
ax.set_title("Layer B per-state coverage (dashed = pre-specified gate)", fontsize=6.6)
ax.legend(frameon=False, fontsize=5.6)
ax = fig.add_subplot(gs[1]); panel_label(ax, "B", dx=-0.16, dy=1.05)
ax.bar(range(len(lb_qc)), [num(r["pool_malignant_n"]) for r in lb_qc], color="#7f8c8d")
ax.set_xticks(range(len(lb_qc)))
ax.set_xticklabels([f"{r['patient']}\n{r['sample']}" for r in lb_qc], rotation=90, fontsize=4.6)
ax.set_ylabel("CNV-inferred malignant cells", fontsize=7)
ax.set_title("Layer B CNV-inferred pool (infercnvpy)", fontsize=6.6)
n_pool = sum(num(r["pool_malignant_n"]) or 0 for r in lb_qc)
n_pat_pool = len(set(r["patient"] for r in lb_qc))
per_state = ", ".join("%s %s/%s" % (r["state"], r["cells"], r["patients"]) for r in lb_state)
_cap = ("Layer B CNV-inferred pool: %d patients, %s cells (infercnvpy). Pre-specified per-state gate: >= 3 cells and >= 2 patients. "
        "Recorded coverage (cells/patients): %s. The frozen decision table (C2_LAYERB_STATE_TESTS.tsv) records the per-patient state test as "
        "not applicable; that rationale matches OPC/NPC (0 cells), and Layer B state-level inference is reported as NOT TESTABLE, with no Layer B "
        "result claimed anywhere in the manuscript or figures." % (n_pat_pool, format(int(n_pool), ","), per_state))
add_caption(fig, [_cap], fontsize=5.6, width=150)
save_fig(fig, "FigS6", subdir=SUP)
log_values("FigS6", [("A", "cells/patients", "C2_LAYERB_STATE_TESTS.tsv", f"{r['state']}|{r['cells']}|{r['patients']}") for r in lb_state])

# ============ S7 random gene-set null controls (Discovery) ============
null_all = load("ss2_cellnull_ALL.tsv", TP)
sigs = []
for r in null_all:
    if "_CCDEP" in r["sig"]:
        continue
    if r["sig"] not in sigs:
        sigs.append(r["sig"])
STATES_N = ["MES", "AC", "OPC", "NPC", "Cycling"]
fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9), gridspec_kw={"width_ratios": [1.5, 1.0, 1.0]})
ax = axes[0]; panel_label(ax, "A", dx=-0.14, dy=1.06)
y = np.arange(len(sigs))[::-1]
obs_mes = [num([r for r in null_all if r["sig"] == s and r["state"] == "MES"][0]["obs"]) for s in sigs]
nul_mes = [num([r for r in null_all if r["sig"] == s and r["state"] == "MES"][0]["null_mean"]) for s in sigs]
sd_mes = [num([r for r in null_all if r["sig"] == s and r["state"] == "MES"][0]["null_sd"]) for s in sigs]
ax.errorbar(nul_mes, y, xerr=sd_mes, fmt="o", ms=3, color="#7f8c8d", label="null mean +/- SD")
ax.plot(obs_mes, y, "s", ms=3.6, color="#c0392b", label="observed")
ax.set_yticks(y); ax.set_yticklabels(sigs, fontsize=5.0)
ax.set_xlabel("mean AUC (MES state, Discovery)", fontsize=6.6)
ax.legend(frameon=False, fontsize=5.4); ax.set_title("Observed vs size/expression-matched null (MES)", fontsize=7)
ax = axes[1]; panel_label(ax, "B", dx=-0.20, dy=1.06)
M = np.array([[num([r for r in null_all if r["sig"] == s and r["state"] == st][0]["z"]) or 0 for st in STATES_N] for s in sigs])
im = ax.imshow(M, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(STATES_N))); ax.set_xticklabels(STATES_N, fontsize=5.2, rotation=45)
ax.set_yticks(range(len(sigs))); ax.set_yticklabels(sigs, fontsize=4.6)
for i in range(len(sigs)):
    for j in range(len(STATES_N)):
        p = num([r for r in null_all if r["sig"] == sigs[i] and r["state"] == STATES_N[j]][0]["emp_p"])
        if p is not None and p < 0.05:
            ax.text(j, i, "*", ha="center", va="center", fontsize=6, color="white")
ax.set_title("z vs null (* P < 0.05)", fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.7)
ax = axes[2]; panel_label(ax, "C", dx=-0.20, dy=1.06)
pb = load("ss2_pseudobulk_null.tsv", TP)
ps = []
for r in pb:
    if r["sig"] not in ps:
        ps.append(r["sig"])
MP = np.array([[num([r for r in pb if r["sig"] == s and r["state"] == st][0]["z"]) or 0 for st in STATES_N] for s in ps])
im2 = ax.imshow(MP, cmap="magma", aspect="auto")
ax.set_xticks(range(len(STATES_N))); ax.set_xticklabels(STATES_N, fontsize=5.2, rotation=45)
ax.set_yticks(range(len(ps))); ax.set_yticklabels(ps, fontsize=4.6)
ax.set_title("Pseudobulk-level null z", fontsize=7)
fig.colorbar(im2, ax=ax, shrink=0.7)
add_caption(fig, ["Cell-level panel: 20 size- and expression-matched random gene sets per signature x state (empirical P). "
                  "Pseudobulk panel: patient-pseudobulk null check. All null controls are Discovery-based and were frozen in the C2 analysis stage."],
            fontsize=5.6, width=160)
save_fig(fig, "FigS7", subdir=SUP)
log_values("FigS7", [("A", "obs/null", "ss2_cellnull_ALL.tsv", f"{s}|MES|{o}|{n}") for s, o, n in zip(sigs, obs_mes, nul_mes)])

# ============ S8 C5 diagnostic detail ============
c5 = [s for s in sigs if s.startswith("C5")]
fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.6), gridspec_kw={"width_ratios": [1.2, 1.2, 0.9]})
for k, (ttl, base, crepl, pair) in enumerate([
        ("Discovery: base vs cell-cycle-depleted", disc, d_cc, ("C5_UP", "AUC_C5_UP_CCDEP")),
        ("Layer A: base vs cell-cycle-depleted", la, la_cc, ("C5_UP", "C5_UP_CCDEP"))]):
    ax = axes[k]; panel_label(ax, "AB"[k], dx=-0.18, dy=1.06)
    for i, st in enumerate(["MES", "OPC", "NPC", "Cycling"]):
        vb = (base.get((pair[0], st), (None, None))[0] or 0) * 1000
        vc = (crepl.get((pair[1], st), (None, None))[0] or 0) * 1000
        ax.bar(i - 0.1, vb, width=0.19, color="#2980b9")
        ax.bar(i + 0.1, vc, width=0.19, color="#c0392b")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(range(4)); ax.set_xticklabels(["MES", "OPC", "NPC", "Cycling"], fontsize=6)
    ax.set_ylabel("delta x 10^-3 (C5_UP)", fontsize=6.8); ax.set_title(ttl, fontsize=7.4)
ax = axes[2]; panel_label(ax, "C", dx=-0.26, dy=1.06)
csig = ["AUC_C5_UP", "AUC_C5_UP_CCDEP", "AUC_C5_DOWN", "AUC_C5_DOWN_CCDEP",
        "AUC_C5_STABLE_UP", "AUC_C5_STABLE_UP_CCDEP", "AUC_C5_STABLE_DOWN", "AUC_C5_STABLE_DOWN_CCDEP"]
fl = [flips.get(s, (0, ""))[0] or 0 for s in csig]
ax.barh(np.arange(len(csig))[::-1], fl, color=["#c0392b" if v else "#1e8449" for v in fl])
ax.set_yticks(np.arange(len(csig))[::-1]); ax.set_yticklabels([s.replace("AUC_", "") for s in csig], fontsize=5.0)
ax.set_xlabel("LOPO top-state changes", fontsize=6.4); ax.set_xlim(0, 3)
ax.set_title("C5 LOPO stability (Discovery)", fontsize=7)
axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#2980b9"), plt.Rectangle((0, 0), 1, 1, color="#c0392b")],
               labels=["base", "CC-depleted"], frameon=False, fontsize=5.4)
add_caption(fig, ["Removal fractions: C5_UP %s%%, C5_DOWN %s%%, C5_STABLE_UP %s%%, C5_STABLE_DOWN %s%% (C2_3_CC_DEPLETION_SUMMARY.tsv). "
                  "C5 remains grade C (incompletely replicated): base MES Layer A contrast is not significant (FDR 0.0563)."
                  % (cc["C5_UP"]["removed_fraction_pct"], cc["C5_DOWN"]["removed_fraction_pct"],
                     cc["C5_STABLE_UP"]["removed_fraction_pct"], cc["C5_STABLE_DOWN"]["removed_fraction_pct"])],
            fontsize=5.6, width=160)
save_fig(fig, "FigS8", subdir=SUP)
log_values("FigS8", [("C", "loo_flips_top", "ss2_lopo.tsv", f"{s}={flips.get(s, (None,))[0]}") for s in csig])

# ============ S9 TCGA/CGGA detail ============
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
ax.set_title("A. Relationship-level forest plot (28 relationships)", fontsize=7.6)
ax = fig.add_subplot(gs[1]); panel_label(ax, "B", dx=-0.13, dy=1.04)
labs, vals = [], []
for r in st2:
    labs.append(r["cohort_scope"].replace("_", " ").replace("Primary only", "primary").replace("WHO IV", "WHO-IV"))
    vals.append(num(str(r["expected_sign"]).split("=")[-1]) or 0)
for coh in ["tcga", "cgga"]:
    v = [num(str(r["expected_sign"]).split("=")[-1]) or 0 for r in gp if r["cohort_scope"] == coh]
    labs.append(f"{coh.upper()} pert.-mean"); vals.append(float(np.mean(v)))
x = np.arange(len(labs))
ax.bar(x, vals, color="#5dade2")
for xi, v in zip(x, vals):
    ax.text(xi, v + 0.2, f"{v:.1f}", ha="center", fontsize=5.4)
ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=4.8, rotation=20, ha="right")
ax.set_ylabel("expected-direction pairs (/21)", fontsize=6.4)
ax.set_title("B. Stratification and gene-perturbation robustness (no stratum improves recovery)", fontsize=7)
ax = fig.add_subplot(gs[2]); panel_label(ax, "C", dx=-0.13, dy=1.06)
ax.set_axis_off()
lines = [f"{r['cohort_scope']}: {r['expected_sign']} ({r['note'][:80]})" for r in cov]
lines += [f"Partial (global-mean adjusted) subset: {sum(1 for r in r3 if str(r['same_direction']).strip()=='True')}/{len(r3)} TCGA relationships keep their frozen sign",
          "CGGA is treated as contextualisation, not replication; sign-reversed pairs are retained and not relabelled."]
for i, ln in enumerate(lines):
    ax.text(0.0, 1.0 - i * 0.26, ln, transform=ax.transAxes, fontsize=6.2, va="top")
ax.set_title("C. Coverage and partial-subset checks", fontsize=7)
save_fig(fig, "FigS9", subdir=SUP)
log_values("FigS9", [("A", "rho", "C3_ROBUSTNESS.tsv/C3_PRIMARY_RESULTS.tsv", f"{r['module_a']}|{r['module_b']}") for r in f1])
print("Supplementary S1-S9 done")
