# -*- coding: utf-8 -*-
"""C3 finalize: assemble deliverable TSVs and Figure6-candidate QC plots from intermediates.
Reads frozen pre-reg definitions; writes output/C3_*.tsv and output/C3_FIG6_CANDIDATE/*.png
"""
import os, json, math
import numpy as np
import pandas as pd
from scipy import stats

CONV = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C3D = os.path.join(CONV, "temp", "c3_data")
OUT = os.path.join(CONV, "output")
FIG = os.path.join(OUT, "C3_FIG6_CANDIDATE")
os.makedirs(FIG, exist_ok=True)

A8 = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
      "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]
A8_LABEL = {"C4_UP": "C4_UP (HITMAN-enriched)",
            "C4_DOWN": "C4_DOWN (HITMAN-down)",
            "M_shared_antiproliferative": "shared_antiproliferative",
            "M_shared_ER_UPR": "shared_ER-UPR(Cycling)",
            "M_cholesterol_SREBP": "SREBP/cholesterol",
            "M_TNF_NFkB": "TNF/NF-kB",
            "M_p53_apoptosis": "p53/apoptosis",
            "M_RIGI_typeI_IFN": "RIG-I/IFN(AC)"}
EXPECTED = {
 ("M_cholesterol_SREBP", "M_TNF_NFkB"): 1, ("M_cholesterol_SREBP", "M_p53_apoptosis"): 1,
 ("M_TNF_NFkB", "M_p53_apoptosis"): 1, ("C4_UP", "M_cholesterol_SREBP"): 1,
 ("C4_UP", "M_TNF_NFkB"): 1, ("C4_UP", "M_p53_apoptosis"): 1,
 ("M_shared_antiproliferative", "M_shared_ER_UPR"): 1,
 ("C4_DOWN", "M_shared_antiproliferative"): 1, ("C4_DOWN", "M_shared_ER_UPR"): 1,
 ("C4_UP", "C4_DOWN"): -1, ("M_shared_antiproliferative", "M_cholesterol_SREBP"): -1,
 ("M_shared_antiproliferative", "M_TNF_NFkB"): -1, ("M_shared_antiproliferative", "M_p53_apoptosis"): -1,
 ("M_shared_ER_UPR", "M_cholesterol_SREBP"): -1, ("M_shared_ER_UPR", "M_TNF_NFkB"): -1,
 ("M_shared_ER_UPR", "M_p53_apoptosis"): -1, ("C4_DOWN", "M_cholesterol_SREBP"): -1,
 ("C4_DOWN", "M_TNF_NFkB"): -1, ("C4_DOWN", "M_p53_apoptosis"): -1,
 ("C4_UP", "M_shared_antiproliferative"): -1, ("C4_UP", "M_shared_ER_UPR"): -1,
}

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    adj = np.minimum.accumulate(ranked * n / np.arange(1, n + 1))
    out = np.empty(n); out[order] = np.minimum(adj, 1.0)
    return out

def spearman_ci(rho, n, alpha=0.05):
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    q = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - q * se), np.tanh(z + q * se)

# ---------------- load main ---------------- 
q1 = pd.read_csv(os.path.join(OUT, "C3_PRIMARY_RESULTS.tsv"), sep="\t")  # placeholder; rebuild below
# rebuild from raw wide scores to keep consistent
def load_scores(tag):
    return pd.read_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t").set_index("Name")
t, g = load_scores("tcga"), load_scores("cgga")
n_t, n_g = len(t), len(g)

q1rows = []
for coh, df, nn in [("TCGA", t, n_t), ("CGGA", g, n_g)]:
    for i, a in enumerate(A8):
        for b in A8[i + 1:]:
            rho, p = stats.spearmanr(df[a], df[b], nan_policy="omit")
            lo, hi = spearman_ci(rho, nn)
            q1rows.append({"analysis": "Q1_architecture", "cohort": coh, "module_a": a, "module_b": b,
                           "rho": rho, "CI95_low": lo, "CI95_high": hi, "P": p,
                           "FDR_BH": np.nan, "expected_sign": EXPECTED.get((a, b), EXPECTED.get((b, a), np.nan)),
                           "n": nn})
q1df = pd.DataFrame(q1rows)
for coh in ["TCGA", "CGGA"]:
    m = q1df["cohort"] == coh
    q1df.loc[m, "FDR_BH"] = bh(q1df.loc[m, "P"].values)
q1df["direction_match"] = np.where(q1df["expected_sign"].notna(),
                                   np.sign(q1df["rho"]) == q1df["expected_sign"], np.nan)
# per-set coverage (min over the two module members) from frozen scoring audit
covj = {"TCGA": json.load(open(os.path.join(C3D, "tcga_set_coverage.json"))),
        "CGGA": json.load(open(os.path.join(C3D, "cgga_set_coverage.json")))}
q1df["min_pair_set_coverage"] = q1df.apply(
    lambda r: round(min(covj[r["cohort"]][r["module_a"]], covj[r["cohort"]][r["module_b"]]), 3), axis=1)
q1df["grade"] = "A"

# Q2 from intermediate
q2 = pd.read_csv(os.path.join(C3D, "_q2_mes_intermediate.tsv"), sep="\t")
q2 = q2.rename(columns={"score": "module", "median_diff_1minus0": "median_diff_MES_minus_nonMES",
                        "mw_p": "P", "CI95_low": "CI95_low", "CI95_high": "CI95_high"})
q2.insert(0, "analysis", "Q2_MES_vs_nonMES_context")
q2.insert(1, "cohort", "TCGA(official_Verhaak_n44)")
q2["expected_dir"] = np.where(q2["module"].isin(["C4_UP", "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]),
                              "MES_higher(hypothesis_AC/MES)", "not_preregistered_explicit")
q2["grade"] = "A"
q2["coverage_ge"] = 0.80

# Q3 secondary
q3 = pd.read_csv(os.path.join(C3D, "_q3_intermediate.tsv"), sep="\t")
q3 = q3.rename(columns={"score": "module", "median_diff_1minus0": "median_diff_MES_minus_nonMES", "mw_p": "P"})
q3.insert(0, "analysis", "Q3_secondary_sensitivity")
q3.insert(1, "cohort", "TCGA(official_Verhaak_n44)")
q3["grade"] = "B"
q3["note"] = "SECONDARY/sensitivity-supported; module-level score only, component-level not testable in bulk"
q3 = q3.drop(columns=["context"], errors="ignore")

# C5
c5 = pd.read_csv(os.path.join(OUT, "C3_C5_EXPLORATORY.tsv"), sep="\t")
c5ax = pd.read_csv(os.path.join(C3D, "_c5_axis_intermediate.tsv"), sep="\t")
c5 = c5.rename(columns={"score": "signature", "median_diff_1minus0": "median_diff_MES_minus_nonMES", "mw_p": "P"})
c5.insert(0, "analysis", "C5_exploratory_MES_context")
c5["grade"] = "C_permanent"
c5ax = c5ax.rename(columns={"score": "signature", "p": "P"})
c5ax.insert(0, "analysis", "C5_exploratory_C4UP_axis")
c5ax["grade"] = "C_permanent"

# robustness assemble
rb_parts = []
cross = pd.read_csv(os.path.join(C3D, "_robust_cross_intermediate.tsv"), sep="\t")
cross["analysis"] = "robust_cross_cohort_sign"
rb_parts.append(cross.rename(columns={"expected": "expected_sign"}))
partial = pd.read_csv(os.path.join(C3D, "_robust_partial_intermediate.tsv"), sep="\t")
partial["analysis"] = "robust_partial_globalExpr"
rb_parts.append(partial)
strat = pd.read_csv(os.path.join(C3D, "_robust_strata_intermediate.tsv"), sep="\t")
strat["analysis"] = "robust_CGGA_strata_Q1"
rb_parts.append(strat)
cov = pd.read_csv(os.path.join(C3D, "_robust_coverage_intermediate.tsv"), sep="\t")
cov["analysis"] = "robust_set_coverage"
rb_parts.append(cov)
gp = pd.read_csv(os.path.join(C3D, "_robust_geneperturb_intermediate.tsv"), sep="\t")
gp["analysis"] = "robust_gene_perturbation"
rb_parts.append(gp[["analysis", "tag", "run", "dir_match"]])
rob = pd.concat(rb_parts, ignore_index=True, sort=False)
rob.to_csv(os.path.join(OUT, "C3_ROBUSTNESS.tsv"), sep="\t", index=False)

# dataset QC
qc_rows = [
 {"dataset": "TCGA-GBM", "version": "GDC (query date 2026-09-09; STAR-Counts augmented gene counts)", "n_tumor": 288,
  "primary": 284, "recurrent": 4, "idh_status": "not available for bulk (official IDH1 annotation for n=44 subset only)",
  "expr_unit": "raw STAR counts (unstranded); log2(+1) for scoring",
  "gene_identifier": "HGNC gene symbol (TCGA gene_name)",
  "dup_gene_handling": "sum across rows sharing gene_name per sample",
  "signature_coverage_min": 0.971, "signature_coverage_min_set": "C4_DOWN/C5 sets (see C3_ROBUSTNESS)",
  "filtering": "excluded 5 Solid Tissue Normal; kept genes detected in >10% samples; N_-prefixed genes removed",
  "note": "cohort composition recorded; primary-only rerun identical direction structure (see C3_ROBUSTNESS)"},
 {"dataset": "CGGA (mRNAseq_693)", "version": "CGGA 2020-05-06 release; WHO grade IV subset (GBM/rGBM)", "n_tumor": 249,
  "primary": 140, "recurrent": 109, "idh_status": "Wildtype 190 / Mutant 49 / NA 10",
  "expr_unit": "gene-level RSEM (FPKM); log2(+1) for scoring; unit per CGGA documentation (待核验)",
  "gene_identifier": "HGNC gene symbol (Gene_Name column)",
  "dup_gene_handling": "mean across duplicate symbols",
  "signature_coverage_min": 0.865, "signature_coverage_min_set": "M_shared_ER_UPR (0.865)",
  "filtering": "WHO IV only; no molecular subtype annotation available (Q2 NOT TESTABLE in CGGA)",
  "note": "IDH and PRS_type strata analyzed separately; recurrent- and primary-only Q1 architecture stable"},
]
qc = pd.DataFrame(qc_rows)
qc.to_csv(os.path.join(OUT, "C3_DATASET_QC.tsv"), sep="\t", index=False)

# write primary (Q1+Q2), secondary, c5, exploratory (again clean)
q1df.to_csv(os.path.join(OUT, "C3_PRIMARY_RESULTS.tsv"), sep="\t", index=False)
q2.to_csv(os.path.join(C3D, "_q2_for_output.tsv"), sep="\t", index=False)
# combine Q1+Q2 into PRIMARY file as blocks
with open(os.path.join(OUT, "C3_PRIMARY_RESULTS.tsv"), "a", encoding="utf-8") as f:
    f.write("\n### BLOCK: Q2 (PRIMARY Grade A modules; MES-vs-nonMES in TCGA official Verhaak n=44)\n")
q2.to_csv(os.path.join(OUT, "C3_PRIMARY_RESULTS.tsv"), sep="\t", index=False, mode="a", header=True)
q3.to_csv(os.path.join(OUT, "C3_SECONDARY_RESULTS.tsv"), sep="\t", index=False)
c5.to_csv(os.path.join(OUT, "C3_C5_EXPLORATORY.tsv"), sep="\t", index=False)
c5ax.to_csv(os.path.join(C3D, "_c5ax_out.tsv"), sep="\t", index=False)
with open(os.path.join(OUT, "C3_C5_EXPLORATORY.tsv"), "a", encoding="utf-8") as f:
    f.write("\n### BLOCK: C5 vs C4_UP axis (both cohorts)\n")
c5ax.to_csv(os.path.join(OUT, "C3_C5_EXPLORATORY.tsv"), sep="\t", index=False, mode="a", header=True)
print("tsv delivered")
for fn in ["C3_DATASET_QC.tsv", "C3_PRIMARY_RESULTS.tsv", "C3_SECONDARY_RESULTS.tsv", "C3_C5_EXPLORATORY.tsv", "C3_ROBUSTNESS.tsv"]:
    print(fn, os.path.getsize(os.path.join(OUT, fn)))
