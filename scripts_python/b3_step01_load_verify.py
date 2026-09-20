"""
B3 - step 01: load & verify (B3.1 sample table fix, B3.2 design-only validation, B3.3 round preprocessing)
GSE324656 20-sample factorial, frozen B2 framework.
"""
import sys, os, json, pickle
import numpy as np
import pandas as pd

BASE = r"${PROJECT_ROOT}"
OUT  = os.path.join(BASE, "output")
TEMP = os.path.join(BASE, "temp")
QC   = os.path.join(OUT, "B3_QC")
os.makedirs(QC, exist_ok=True)

def log(*a):
    print(*a, flush=True)

# ---------- B3.1 load sample table, fix MAN factor typo for ND-F rows ----------
tab = pd.read_csv(os.path.join(OUT, "B2_FROZEN_SAMPLE_TABLE.tsv"), sep="\t")
# authoritative condition from B1 master sheet (GSM->group_short)
b1 = pd.read_csv(os.path.join(OUT, "B1_MASTER_SAMPLE_SHEET.tsv"), sep="\t")
b1_map = dict(zip(b1["matrix_col"], b1["group_short"]))
tab["condition_authoritative"] = tab["matrix_col"].map(b1_map)
# factor encoding: ND-F condition must be MAN=ND, Field=F
expect_man = np.where(tab["condition_authoritative"].str.startswith("MAN"), "MAN", "ND")
expect_fld = np.where(tab["condition_authoritative"].str.endswith("-F"), "F", "NF")
tab["MAN_orig"] = tab["MAN"]
tab["Field_orig"] = tab["Field"]
n_bad = int((tab["MAN"] != expect_man).sum() + (tab["Field"] != expect_fld).sum())
tab["MAN"] = expect_man
tab["Field"] = expect_fld
# frozen n per group
n_expected = {"ND-NF": 6, "MAN-NF": 5, "ND-F": 4, "MAN-F": 5}
n_obs = tab.groupby("condition_authoritative").size().to_dict()
assert all(n_obs[k] == v for k, v in n_expected.items()), (n_obs, n_expected)

# ---------- load counts ----------
genes = pd.read_csv(os.path.join(TEMP, "genes.txt"), header=None, names=["gene_id"])
cols  = pd.read_csv(os.path.join(TEMP, "cols.txt"), header=None, names=["matrix_col"])
mat   = np.load(os.path.join(TEMP, "counts_mat.npy"))   # float64
log("counts shape", mat.shape, "genes", len(genes), "cols", len(cols))
assert mat.shape == (len(genes), len(cols))
assert (cols["matrix_col"].values == tab["matrix_col"].values).all()

# verify matrix col order matches sample table
# counts non-integer proportion (Salmon estimated counts)
frac_nonint = float(np.mean(~np.equal(np.mod(mat, 1), 0)))
n_genes_affected = int((np.mod(mat, 1) != 0).any(axis=1).sum())
log("non-integer entries:", round(frac_nonint * 100, 4), "% ; genes affected:", n_genes_affected)

# B3.1 report object
b31 = {
    "balanced_wording_found_in_B2_B1_files": [],   # filled by grep externally
    "sample_table_factor_fix": {
        "issue": "B2_FROZEN_SAMPLE_TABLE.tsv ND-F group (T.Cells_M/N/P/Q) MAN column was 'MAN'; conflicts with condition=ND-F (no device). Authoritative group_short from B1_MASTER_SAMPLE_SHEET confirms ND-F.",
        "action": "MAN column set to 'ND' for the 4 ND-F rows; membership unchanged; condition/group_short unchanged.",
        "rows_fixed": int((tab["MAN_orig"] != expect_man).sum()) if n_bad > 0 else 0,
        "n_bad_any": n_bad
    },
    "group_n": n_obs
}
log("B3.1 group_n", n_obs)

# ---------- B3.2 design-only validation ----------
# metadata with factors MAN (ND ref) & Field (NF ref)
meta = pd.DataFrame({
    "sample": tab["matrix_col"].values,
    "group": tab["condition_authoritative"].values,
    "MAN":   tab["MAN"].values,
    "Field": tab["Field"].values,
}, index=tab["matrix_col"].values)
meta["MAN"] = pd.Categorical(meta["MAN"], categories=["ND", "MAN"])
meta["Field"] = pd.Categorical(meta["Field"], categories=["NF", "F"])

# build design matrix with formulaic as pydeseq2 would: ~MAN*Field
try:
    from pydeseq2.dds import FormulaicContrasts
    fc = FormulaicContrasts(meta, "~MAN*Field")
    X = fc.design_matrix.copy()
    log("B3.2 design_matrix columns:", list(X.columns))
    log("X shape", X.shape)
    rank = int(np.linalg.matrix_rank(X.values))
    cond = float(np.linalg.cond(X.values))
    log("rank", rank, "cond", cond)
    # sample-to-design mapping
    X["sample"] = X.index if hasattr(X.index, "__len__") else None
except Exception as e:
    log("formulaic error", repr(e))
    X = None; rank = None; cond = None

meta.to_pickle(os.path.join(TEMP, "b3_meta.pkl"))
pd.DataFrame({"gene_id": genes["gene_id"].values}).to_csv(os.path.join(TEMP, "b3_genes.tsv"), sep="\t", index=False)
np.save(os.path.join(TEMP, "b3_counts_raw.npy"), mat)

# ---------- B3.3 round preprocessing (Strategy B) ----------
rounded = np.round(mat).astype(np.int64)
# round() halves away from zero (R style round() is banker's? R round() is IEC 60559 round-half-to-even). Check deviation only.
round_diff = rounded - mat
lib_pre  = mat.sum(axis=0)
lib_post = rounded.sum(axis=0)
absdiff = np.abs(round_diff).sum(axis=0)
maxdev  = np.abs(round_diff).max()
affected_entries = float(np.mean(round_diff != 0))
affected_genes   = int((round_diff != 0).any(axis=1).sum())
log("B3.3 affected entries", round(affected_entries*100,4), "% ; genes", affected_genes, "; maxdev", maxdev)
np.save(os.path.join(TEMP, "b3_counts_rounded.npy"), rounded)

prep = pd.DataFrame({
    "sample": cols["matrix_col"].values,
    "group": tab["condition_authoritative"].values,
    "libsize_raw": lib_pre,
    "libsize_rounded": lib_post,
    "libsize_absdiff": lib_post - lib_pre,
    "libsize_reldiff": (lib_post - lib_pre) / lib_pre,
    "entry_absdiff_sum": absdiff,
})
prep.to_csv(os.path.join(TEMP, "b3_round_diagnostics.tsv"), sep="\t", index=False)

# save objects for next steps
with open(os.path.join(TEMP, "b3_prep.pkl"), "wb") as f:
    pickle.dump({
        "sample_table": tab, "b1_map": b1_map,
        "b31": b31, "design_cols": list(X.columns) if X is not None else None,
        "design_rank": rank, "design_cond": cond,
        "frac_nonint": frac_nonint, "n_genes_affected": n_genes_affected,
        "affected_entries": affected_entries, "affected_genes": affected_genes,
        "maxdev": maxdev,
    }, f)
log("B3.2 design matrix saved; rank=%s cond=%s" % (rank, cond))
log("OK step01")
