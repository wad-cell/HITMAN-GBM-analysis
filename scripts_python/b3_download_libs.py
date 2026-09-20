"""download gene set libraries (Enrichr) and cache to temp; quick check"""
import os, traceback
import gseapy as gp
BASE = r"${PROJECT_ROOT}"
TEMP = os.path.join(BASE, "temp")
LOG = os.path.join(TEMP, "b3_gselib_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)
libs = {
    "HALLMARK": "MSigDB_Hallmark_2020",
    "REACTOME": "Reactome_2022",
    "GOBP": "GO_Biological_Process_2021",
}
try:
    for key, lib in libs.items():
        gs = gp.get_library(name=lib)
        log(key, "sets", len(gs))
        path = os.path.join(TEMP, f"b3_gslib_{key}.gmt")
        with open(path, "w", encoding="utf-8") as f:
            for name, genes in gs.items():
                f.write(name + "\t\t" + "\t".join(genes) + "\n")
        log(key, "saved", path)
    log("LIB DONE")
except Exception:
    log("LIB EXC", traceback.format_exc())
