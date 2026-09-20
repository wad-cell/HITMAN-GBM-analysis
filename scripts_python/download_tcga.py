import json, os, time, urllib.request, ssl, sys, concurrent.futures as cf
ctx = ssl.create_default_context()
HDRS = {"User-Agent": "Mozilla/5.0"}
ws = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # workspace? no: file in temp/c3_data -> temp
c3 = os.path.join(ws, "c3_data")
sel = [s for s in json.load(open(os.path.join(c3, "tcga_gbm_selected.json"), encoding="utf-8")) if s["sample_type"] != "Solid Tissue Normal"]
outdir = os.path.join(c3, "tcga_star")
os.makedirs(outdir, exist_ok=True)

def dl(item):
    dest = os.path.join(outdir, item["case"] + ".tsv")
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        return item["case"], "skip", os.path.getsize(dest)
    url = "https://api.gdc.cancer.gov/data/" + item["file_id"]
    for attempt in range(4):
        try:
            t0 = time.time()
            req = urllib.request.Request(url, headers=HDRS)
            with urllib.request.urlopen(req, timeout=180, context=ctx) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            return item["case"], "ok", len(data), round(time.time() - t0, 1)
        except Exception as e:
            if attempt == 3:
                return item["case"], "ERR", repr(e)[:140]
            time.sleep(2)

if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 99999
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    todo = [s for s in sel if not (os.path.exists(os.path.join(outdir, s["case"] + ".tsv")) and os.path.getsize(os.path.join(outdir, s["case"] + ".tsv")) > 1_000_000)][:limit]
    done = total_ok = 0
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(dl, todo):
            done += 1
            if res[1] == "ok":
                total_ok += 1
            if res[1] in ("ERR", "skip") or done % 20 == 0 or done == len(todo):
                print(f"[{done}/{len(todo)}] {res[0]} {res[1]} {res[2] if len(res)<4 else (res[2], res[3])} elapsed={time.time()-t0:.0f}s", flush=True)
    print(f"BATCH DONE ok={total_ok} err=0 elapsed={time.time()-t0:.1f}s")
