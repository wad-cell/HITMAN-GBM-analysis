import os, json
import pandas as pd
from scipy import stats

CONV = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEMP = os.path.join(CONV, "temp")
C3D = os.path.join(TEMP, "c3_data")
OUT = os.path.join(CONV, "output")

def load_gmt(path):
    sets = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            sets[p[0]] = [g for g in p[2:] if g]
    return sets
hm = load_gmt(os.path.join(TEMP, "b3_gslib_HALLMARK.gmt"))
rt = load_gmt(os.path.join(TEMP, "b3_gslib_REACTOME.gmt"))

def uni(*names):
    s = set()
    for src, nm in names:
        s |= set(src[nm])
    return sorted(s)

modules = {
    "M_shared_antiproliferative": uni((hm, "E2F Targets"), (hm, "Myc Targets V1"), (hm, "G2-M Checkpoint")),
    "M_shared_ER_UPR": uni((hm, "Unfolded Protein Response"), (rt, "Unfolded Protein Response (UPR) R-HSA-381119")),
    "M_cholesterol_SREBP": uni((hm, "Cholesterol Homeostasis"), (rt, "Cholesterol Biosynthesis R-HSA-191273"), (rt, "Regulation Of Cholesterol Biosynthesis By SREBP (SREBF) R-HSA-1655829")),
    "M_TNF_NFkB": uni((hm, "TNF-alpha Signaling via NF-kB"), (rt, "TNFR1-induced NFkappaB Signaling Pathway R-HSA-5357956")),
    "M_p53_apoptosis": uni((hm, "p53 Pathway"), (hm, "Apoptosis")),
    "M_RIGI_typeI_IFN": uni((rt, "DDX58/IFIH1-mediated Induction Of Interferon-Alpha/Beta R-HSA-168928"), (rt, "Interferon Alpha/Beta Signaling R-HSA-909733")),
    "M_HSF1_heat_shock": uni((rt, "HSF1 Activation R-HSA-3371511"), (rt, "HSF1-dependent Transactivation R-HSA-3371571"), (rt, "Regulation Of HSF1-mediated Heat Shock Response R-HSA-3371453")),
    "M_autophagy_lysosome": uni((rt, "Macroautophagy R-HSA-1632852"), (rt, "Lysosome Vesicle Biogenesis R-HSA-432720")),
}
expected = {"M_shared_antiproliferative": 480, "M_shared_ER_UPR": 151, "M_cholesterol_SREBP": 121, "M_TNF_NFkB": 220,
            "M_p53_apoptosis": 340, "M_RIGI_typeI_IFN": 136, "M_HSF1_heat_shock": 99, "M_autophagy_lysosome": 146}
assert all(len(modules[k]) == v for k, v in expected.items()), "module rebuild mismatch"

# C4/C5/stable from frozen signature gene table (C2_2)
sig = pd.read_csv(os.path.join(OUT, "C2_2_FROZEN_SIGNATURE_GENES.tsv"), sep="\t")
sig["symbol"] = sig["symbol"].fillna("")
sig = sig[sig["symbol"] != ""]
sig = sig.drop_duplicates(subset=["signature", "symbol"])
sets = {}
for name in ["C4_UP", "C4_DOWN", "C5_UP", "C5_DOWN", "C5_STABLE_UP", "C5_STABLE_DOWN"]:
    sets[name] = sorted(sig.loc[sig["signature"] == name, "symbol"].unique().tolist())
    print(name, "n unique symbols:", len(sets[name]))
# merge
all_sets = {**modules, **sets}
direction = {}
for k in modules:
    direction[k] = 0  # directionless module; signed analysis uses context-agnostic score
direction["C4_UP"] = 1; direction["C4_DOWN"] = -1
direction["C5_UP"] = 1; direction["C5_DOWN"] = -1
direction["C5_STABLE_UP"] = 1; direction["C5_STABLE_DOWN"] = -1
with open(os.path.join(C3D, "frozen_sets.json"), "w") as f:
    json.dump({"sets": all_sets, "direction": direction}, f)
# write a tsv registry for the report
rows = []
for k in all_sets:
    rows.append({"set": k, "n_unique_symbols": len(all_sets[k]), "direction": direction[k]})
pd.DataFrame(rows).to_csv(os.path.join(C3D, "frozen_set_registry.tsv"), sep="\t", index=False)
print("frozen set registry saved; total sets:", len(all_sets))
