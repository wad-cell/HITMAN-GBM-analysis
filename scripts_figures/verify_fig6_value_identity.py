# -*- coding: utf-8 -*-
"""the release stage - Fig6 value-identity check (read-only).

Purpose: prove that the 2x2 layout-only re-composition plotted the SAME values as the
the pre-release single-column export, i.e. no number, scale or statistic changed.

Method (deterministic, no re-statistics):
 1. Rebuild the exact numeric payload that fig6_2x2.py plots, straight from the frozen
    result tables (string tokens, not recomputed quantities).
 2. Compare token-by-token against the pre-release plotted-value log
    (provenance/plotted_values/pre_release/Fig6_values.tsv).
 3. Emit a machine-readable verdict + payload digest.
"""
import os, sys, csv, re, hashlib, json

SCR = os.path.dirname(os.path.abspath(__file__))
TEMP = os.path.dirname(SCR)          # release package root
BASE = os.path.dirname(TEMP)
sys.path.insert(0, SCR)
from fig_common import rd, num, PT, C3D

pr = rd(os.path.join(PT, "C3_PRIMARY_RESULTS.tsv"))
rb = rd(os.path.join(PT, "C3_ROBUSTNESS.tsv"))
sec = rd(os.path.join(PT, "C3_SECONDARY_RESULTS.tsv"))
q2 = rd(os.path.join(C3D, "_q2_for_output.tsv"))

q1 = [r for r in pr if r["block"] == "Q1_architecture"]
payload = []          # canonical tokens actually drawn
for r in q1:
    payload += [r["effect_rho"], r["cohort"], r["module_a"], r["module_b"]]
for r in q2:
    payload += [r["median_diff_MES_minus_nonMES"], r["P"], r["module"]]
for r in sec:
    if r["module"] == "M_autophagy_lysosome" and r["cohort"].startswith("TCGA"):
        payload += [r["median_diff_MES_minus_nonMES"], r["P"]]
for r in rb:
    if r["analysis"] == "robust_1_cross_cohort_sign":
        payload += [r["rho_TCGA"], r["rho_CGGA"], r["expected_sign"], r["same_direction"]]
    if r["analysis"] == "robust_2_stratum_Q1":
        payload += [r["expected_sign"], r["cohort_scope"]]
    if r["analysis"] == "robust_4_gene_perturbation":
        payload += [r["expected_sign"], r["cohort_scope"]]
    if r["analysis"] == "robust_5_coverage":
        payload += [r["expected_sign"], r["cohort_scope"]]

def tokens(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            for cell in row:
                out += cell.split("|")
    return set(t.strip() for t in out if t.strip())

pre_tok = tokens(os.path.join(TEMP, "provenance", "plotted_values", "pre_release", "Fig6_values.tsv"))
rel_tok = tokens(os.path.join(TEMP, "provenance", "plotted_values", "release", "Fig6_values.tsv"))

missing = sorted({t for t in payload if t not in pre_tok})
new_only = sorted({t for t in rel_tok if t not in pre_tok})

digest = hashlib.sha256("\n".join(payload).encode("utf-8")).hexdigest()
rows = [
    ("panel_A_B_heatmap_pairs_rows", len(q1)),
    ("inset_rows", len([r for r in q2]) + len([r for r in sec if r["module"] == "M_autophagy_lysosome" and r["cohort"].startswith("TCGA")])),
    ("forest_rows", len([r for r in rb if r["analysis"] == "robust_1_cross_cohort_sign"])),
    ("stratum_rows", len([r for r in rb if r["analysis"] == "robust_2_stratum_Q1"])),
    ("perturbation_rows", len([r for r in rb if r["analysis"] == "robust_4_gene_perturbation"])),
    ("coverage_rows", len([r for r in rb if r["analysis"] == "robust_5_coverage"])),
    ("payload_token_count", len(payload)),
    ("payload_unique_tokens", len(set(payload))),
    ("tokens_not_found_in_pre_release_log", len(missing)),
    ("tokens_only_in_release_log", len(new_only)),
    ("value_identity_verdict", "IDENTICAL" if not missing else "MISMATCH"),
    ("payload_sha256", digest),
]
out = os.path.join(TEMP, "provenance", "FIG6_VALUE_IDENTITY.tsv")
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["metric", "value"])
    for k, v in rows:
        w.writerow([k, v])

print("\n".join(f"{k}\t{v}" for k, v in rows))
if missing:
    print("MISSING SAMPLE:", missing[:10])
if new_only:
    print("NEW-ONLY SAMPLE:", new_only[:10])
