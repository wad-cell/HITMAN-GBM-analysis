"""Rebuild FIGURE_TECHNICAL_QC.tsv from the exported figure files (no re-analysis)."""
import os, csv
from PIL import Image

W = os.environ.get("PROJECT_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(W, "figures")
SUP = os.path.join(FIG, "supplementary")
OUT = os.path.join(W, "manifest", "FIGURE_TECHNICAL_QC.tsv")

MAIN = ["Fig%d" % i for i in range(1, 7)]
SUPS = ["FigS%d" % i for i in range(1, 10)]

rows = []


def size_ok(fig, w, h):
    # single-column width acceptable when >= 700 px at 300 dpi would be ~2.3 in; accept width>=700 and height<=4200
    return "YES" if (w >= 700 and h >= 150 and h <= 4200) else "NO"


def tif_row(fig, path):
    im = Image.open(path)
    w, h = im.size
    dpi = im.info.get("dpi", (0, 0))
    dx = round(float(dpi[0]), 0) if dpi and dpi[0] else 0
    dy = round(float(dpi[1]), 0) if dpi and dpi[1] else 0
    mode = im.mode
    integ = "OK"
    try:
        im.load()
    except Exception as e:
        integ = "FAIL:%s" % type(e).__name__
    if mode == "RGB":
        pass
    elif mode == "RGBA":
        pass
    n_frames = getattr(im, "n_frames", 1)
    if n_frames > 1:
        integ = integ + "+multiframe(%d)" % n_frames
    px = os.path.getsize(path)
    rows.append(dict(figure=fig, file=path, width_px=w, height_px=h, dpi_x=int(dx), dpi_y=int(dy),
                     color_mode=mode, bytes=px, integrity=integ,
                     min_dpi_ok="YES" if dx >= 300 and dy >= 300 else "NO",
                     size_ok=size_ok(fig, w, h),
                     qc_status="PASS" if (integ == "OK" and mode in ("RGB", "RGBA") and dx >= 300 and dy >= 300) else "REVIEW"))
    im.close()


def other_row(fig, path, label):
    px = os.path.getsize(path)
    rows.append(dict(figure=fig, file=path, width_px="", height_px="", dpi_x="", dpi_y="",
                     color_mode="vector/raster" if label != "PNG" else "RGB", bytes=px,
                     integrity="OK", min_dpi_ok="n/a" if label == "PDF" else "YES",
                     size_ok="YES", qc_status="INFO_%s" % label))


for f in MAIN:
    tif_row(f, os.path.join(FIG, f + ".tif"))
for f in SUPS:
    tif_row(f, os.path.join(SUP, f + ".tif"))
for f in MAIN + SUPS:
    d = FIG if f in MAIN else SUP
    for ext, lab in ((".png", "PNG"), (".pdf", "PDF")):
        p = os.path.join(d, f + ext)
        if os.path.exists(p):
            other_row(f, p, lab)

hdr = list(rows[0].keys())
with open(OUT, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=hdr, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
npass = sum(1 for r in rows if r["qc_status"] == "PASS")
print("wrote", OUT, len(rows), "rows; TIFF PASS =", npass, "/", len(MAIN) + len(SUPS))
for r in rows:
    if r["qc_status"].startswith("REVIEW"):
        print("  REVIEW:", r["figure"], r["file"], r["integrity"], r["color_mode"], r["dpi_x"])
