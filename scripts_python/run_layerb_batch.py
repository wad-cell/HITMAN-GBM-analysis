# C2.11 Layer B batch: run CNV inference on non-official-malignant cells only
import os, sys, time, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import c2_layerb_run

T = os.path.dirname(os.path.abspath(__file__))
SAMPLES = ["BT322", "BT338_1of2", "BT338_2of2", "BT363_1of2", "BT363_2of2",
           "BT364_2of2", "BT390", "BT397_1of2", "BT397_2of2", "BT407"]

def main():
    meta = pd.read_csv(os.path.join(T, "couturier_meta_layers.tsv"), sep="\t")
    print("driver start", time.strftime("%H:%M:%S"), flush=True)
    for i, s in enumerate(SAMPLES, 1):
        n_un = int((~meta[meta["sample"] == s]["mal_annot"]).sum())
        t0 = time.time()
        print(f"[{i}/{len(SAMPLES)}] {s} unannotated={n_un} start", flush=True)
        try:
            c2_layerb_run.main(s, prototype=False)
            print(f"[{i}/{len(SAMPLES)}] {s} done in {round(time.time()-t0,1)}s", flush=True)
        except Exception:
            traceback.print_exc()
            print(f"[{i}/{len(SAMPLES)}] {s} FAILED", flush=True)
    print("driver done", time.strftime("%H:%M:%S"), flush=True)

if __name__ == "__main__":
    main()
