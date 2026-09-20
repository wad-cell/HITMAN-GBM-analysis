import os, json, time
import numpy as np
import pandas as pd

ws = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # workspace root via conv dir? this file is in temp/c3_data -> root = workspace? see below
# file lives at <ws>\temp\c3_data\prep_matrices.py -> conv root = dirname(dirname(dirname(abspath)))
CONV = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEMP = os.path.join(CONV, "temp")
C3D = os.path.join(TEMP, "c3_data")
OUT = os.path.join(CONV, "output")
os.makedirs(C3D, exist_ok=True)

t0 = time.time()
# ---------------- TCGA ----------------
star_dir = os.path.join(C3D, "tcga_star")
files = sorted([f for f in os.listdir(star_dir) if f.endswith(".tsv")])
print("TCGA files:", len(files))
first = pd.read_csv(os.path.join(star_dir, files[0]), sep="\t", comment="#", dtype=str)
first = first[["gene_name", "gene_type", "unstranded"]].copy()
first.columns = ["gene_name", "gene_type", "v0"]
first["gene_name"] = first["gene_name"].fillna("")
keep_mask = first["gene_name"].notna() & (first["gene_name"] != "") & (~first["gene_name"].str.startswith("N_"))
first = first[keep_mask].reset_index(drop=True)
# aggregate duplicate gene_name in the first file to establish unique gene rows
agg0 = first.groupby("gene_name", sort=True).agg(gene_type=("gene_type", "first")).reset_index()
gene_names = agg0["gene_name"].tolist()
gene_types = agg0["gene_type"].tolist()
G = len(gene_names)
gidx = {g: i for i, g in enumerate(gene_names)}
print("TCGA unique gene rows:", G)
M = np.zeros((len(files), G), dtype=np.float32)
samples = []
for j, fn in enumerate(files):
    sample = fn[:-4]
    samples.append(sample)
    df = pd.read_csv(os.path.join(star_dir, fn), sep="\t", comment="#",
                     usecols=["gene_name", "unstranded"], dtype={"gene_name": str, "unstranded": np.float64})
    df["gene_name"] = df["gene_name"].fillna("")
    df = df[df["gene_name"].notna() & (df["gene_name"] != "") & (~df["gene_name"].str.startswith("N_"))]
    s = df.groupby("gene_name", sort=False)["unstranded"].sum()
    # map into row vector
    for name, v in s.items():
        i = gidx.get(name)
        if i is not None:
            M[j, i] = v
    if (j + 1) % 50 == 0:
        print(f"  parsed {j+1}/{len(files)} {time.time()-t0:.0f}s", flush=True)
# keep genes with nonzero count in >10% samples
det = (M > 0).sum(axis=0)
m_keep = det > max(5, len(files) * 0.1)
M = M[:, m_keep]
gene_names2 = [g for g, k in zip(gene_names, m_keep) if k]
gene_types2 = [g for g, k in zip(gene_types, m_keep) if k]
print("TCGA final genes:", M.shape, time.time()-t0)
np.savez_compressed(os.path.join(C3D, "tcga_matrix.npz"), mat=M, samples=np.array(samples, dtype=object))
with open(os.path.join(C3D, "tcga_genes.json"), "w") as f:
    json.dump({"gene": gene_names2, "type": gene_types2}, f)

# ---------------- CGGA ----------------
expr = pd.read_csv(os.path.join(C3D, "CGGA.mRNAseq_693.RSEM-genes.txt"), sep="\t")
expr = expr.rename(columns={"Gene_Name": "gene"})
expr = expr[expr["gene"].notna()].groupby("gene", sort=True).mean().reset_index()
clin = pd.read_csv(os.path.join(C3D, "CGGA.mRNAseq_693_clinical.txt"), sep="\t")
clin.columns = [c.strip() for c in clin.columns]
clin = clin.rename(columns={"CGGA_ID": "id", "Grade": "grade", "Histology": "hist", "PRS_type": "prs_type",
                            "IDH_mutation_status": "idh", "1p19q_codeletion_status": "code19",
                            "OS": "os", "Censor (alive=0; dead=1)": "censor"})
# WHO IV subset
iv = clin[clin["grade"].astype(str).str.replace("WHO ", "").str.strip().isin(["IV", "IV "] ) | clin["grade"].astype(str).str.contains("IV")].copy()
iv = clin[clin["grade"].astype(str).str.contains("IV", na=False)].copy()
print("CGGA total samples:", clin.shape[0], "WHO IV:", iv.shape[0])
col_iv = [c for c in iv["id"].tolist() if c in expr.columns]
print("CGGA IV samples in expr:", len(col_iv))
E = expr[["gene"] + col_iv].copy()
Emat = E.drop(columns=["gene"]).to_numpy(dtype=np.float32).T
print("CGGA expr shape:", Emat.shape, "genes:", E.shape[0])
np.savez_compressed(os.path.join(C3D, "cgga_matrix.npz"), mat=Emat, samples=np.array(col_iv, dtype=object),
                    genes=E["gene"].tolist())
iv[["id", "prs_type", "hist", "grade", "age" if "age" in iv.columns else "id", "idh", "code19", "os", "censor"]].to_csv(
    os.path.join(C3D, "cgga_clinical_iv.tsv"), sep="\t", index=False)
print("done", time.time() - t0)
