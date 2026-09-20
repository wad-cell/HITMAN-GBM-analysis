import os, json, time
import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats

CONV = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C3D = os.path.join(CONV, "temp", "c3_data")
frozen = json.load(open(os.path.join(C3D, "frozen_sets.json"), encoding="utf-8"))
sets = frozen["sets"]

def load_matrix(tag):
    z = np.load(os.path.join(C3D, f"{tag}_matrix.npz"), allow_pickle=True)
    mat = z["mat"]; samples = list(z["samples"])
    if tag == "tcga":
        genes = json.load(open(os.path.join(C3D, "tcga_genes.json")))["gene"]
    else:
        genes = list(z["genes"])
    return mat, samples, genes

def run(tag):
    mat, samples, genes = load_matrix(tag)
    # log transform (counts: log2(count+1); fpkm: log2(fpkm+1)); keep values finite
    X = np.log2(mat.astype(np.float64) + 1.0)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    df = pd.DataFrame(X.T, index=genes, columns=samples)
    # filter zero-expression rows entirely? keep as is for ranking stability but drop all-zero rows (uninformative ties)
    df = df[(df.sum(axis=1) > 0)]
    # restrict gene_sets to those with enough members present
    gene_universe = set(df.index)
    gsets = {}
    for k, gs in sets.items():
        present = [g for g in gs if g in gene_universe]
        gsets[k] = present
        print(tag, k, "n_set", len(gs), "present_in_matrix", len(present),
              "coverage=%.3f" % (len(present) / len(gs)))
    with open(os.path.join(C3D, f"{tag}_set_coverage.json"), "w") as f:
        json.dump({k: len(v) / len(sets[k]) for k, v in gsets.items()}, f)
    res = gp.ssgsea(data=df, gene_sets=gsets, outdir=None, sample_norm_method="rank",
                    min_size=2, max_size=20000, permutation_num=0, verbose=False)
    long = res.res2d[["Name", "Term", "ES"]].copy()
    long.to_csv(os.path.join(C3D, f"{tag}_ssgsea_long.tsv"), sep="\t", index=False)
    wide = long.pivot_table(index="Name", columns="Term", values="ES")
    wide.to_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t")
    print(tag, "samples scored:", wide.shape)
    return wide

t0 = time.time()
tcga_score = run("tcga")
print("TCGA ssgsea done", time.time() - t0)
t0 = time.time()
cgga_score = run("cgga")
print("CGGA ssgsea done", time.time() - t0)
