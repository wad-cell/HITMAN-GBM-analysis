"""map ENSG -> gene symbol via mygene (batched POST), resume-capable, retries"""
import os, json, time, urllib.request, traceback
import pandas as pd
BASE = r"${PROJECT_ROOT}"
TEMP = os.path.join(BASE, "temp")
LOG = os.path.join(TEMP, "b3_map_log.txt")
OUTP = os.path.join(TEMP, "b3_ensg_symbol_map.tsv")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)
try:
    genes = pd.read_csv(os.path.join(TEMP, "b3_genes.tsv"), sep="\t")["gene_id"].tolist()
    core = [g.split(".")[0] for g in genes]
    uniq = sorted(set(core))
    done = set()
    if os.path.exists(OUTP):
        df = pd.read_csv(OUTP, sep="\t", dtype=str)
        done = set(df["gene_id"].tolist())
        log("resume with existing", len(done))
    BATCH = 500
    pending = [g for g in uniq if g not in done]
    log("pending", len(pending), "of", len(uniq))
    fh = open(OUTP, "a", encoding="utf-8") if done else open(OUTP, "w", encoding="utf-8")
    if not done:
        fh.write("gene_id\tsymbol\n")
    ok = 0; fail = 0
    for i in range(0, len(pending), BATCH):
        batch = pending[i:i+BATCH]
        q = json.dumps({"q": batch, "scopes": "ensembl.gene", "fields": "symbol", "species": "human"})
        got = False
        for attempt in range(4):
            try:
                req = urllib.request.Request("https://mygene.info/v3/query", data=q.encode(),
                                             headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    res = json.loads(r.read().decode())
                for hit in res:
                    qid = hit.get("query"); sym = hit.get("symbol")
                    if qid and sym:
                        fh.write(f"{qid}\t{sym}\n")
                        ok += 1
                got = True
                break
            except Exception as e:
                log("retry", i, "attempt", attempt, str(e)[:100])
                time.sleep(2 + 4*attempt)
        if not got:
            fail += len(batch)
            log("FAILED batch", i)
        if (i // BATCH) % 4 == 0:
            fh.flush()
            log("progress", i + len(batch), "/", len(pending), "ok", ok, "fail", fail)
        time.sleep(0.1)
    fh.close()
    log("MAP DONE total ok", ok, "fail", fail)
except Exception:
    log("MAP EXC", traceback.format_exc())
