# -*- coding: utf-8 -*-
"""release layout QC (read-only, plotting-only).

Re-uses the QC logic (scripts_figures/_layout_qc_figures.py) but audits only the
the release stage figure script.  Reports per figure: clipped text, overlapping annotations,
crowding of long tick-label sets, panel letters present, canvas size.
No data are read or recomputed; matplotlib objects are only inspected.
"""
import runpy, sys, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SC)
NUM = re.compile(r"^[<>=]?\s*[-+]?[\d.,]+\s*\*?$")
LETTER = re.compile(r"^[A-Z]$")
# after the layout pass, panel letters are folded into the axes header as "A   <title>"
LETTERPRE = re.compile(r"^([A-Z])\s{2,}\S")
REPORT = []
FIGIX = {"n": 0}

FIGLIST = ["fig6_2x2.py"]
QC_DIR = os.path.join(os.path.dirname(SC), "provenance", "qc_layout")
os.makedirs(QC_DIR, exist_ok=True)
OUT = os.path.join(QC_DIR, "layout_qc_fig6_2x2.txt")


def _renderer(fig):
    try:
        return fig.canvas.get_renderer()
    except Exception:
        return fig._get_renderer()


def _owner(fig, t):
    for i, ax in enumerate(fig.axes):
        if t is getattr(ax, "title", None):
            return "ax%d.title[%s]" % (i, (ax.get_title() or "")[:28])
        if t in list(ax.get_xticklabels()):
            return "ax%d.xtick[%s]" % (i, (ax.get_title() or "")[:28])
        if t in list(ax.get_yticklabels()):
            return "ax%d.ytick[%s]" % (i, (ax.get_title() or "")[:28])
        if t in list(ax.texts):
            return "ax%d.text[%s]" % (i, (ax.get_title() or "")[:28])
    return "fig.text"


def audit(fig, tag):
    if id(fig) in FIGIX:
        return
    FIGIX[id(fig)] = True
    fig.canvas.draw()
    r = _renderer(fig)
    W, H = fig.canvas.get_width_height()
    FIGIX["n"] += 1
    fid = "F%02d" % FIGIX["n"]
    titles = [a.get_title() for a in fig.axes if a.get_title()]
    outside, items = [], []
    for t in fig.findobj(matplotlib.text.Text):
        s = (t.get_text() or "").strip()
        if not s or not t.get_visible():
            continue
        items.append((s, t, t.get_window_extent(renderer=r)))
    for s, t, bb in items:
        if bb.x1 > W + 1 or bb.x0 < -1 or bb.y1 > H + 1 or bb.y0 < -1:
            outside.append("CLIPPED  %-46s | %-42s | x=[%.0f,%.0f] y=[%.0f,%.0f] W=%d H=%d"
                           % (s[:46].replace("\n", " "), _owner(fig, t), bb.x0, bb.x1, bb.y0, bb.y1, W, H))
    over = []
    big = [(s, bb) for s, t, bb in items if (len(s) >= 8 and not NUM.match(s)) or LETTER.match(s)]
    for i in range(len(big)):
        for j in range(i + 1, len(big)):
            sa, a = big[i]
            sb, b = big[j]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > 2 and oy > 2:
                area = ox * oy
                small = min((a.x1 - a.x0) * (a.y1 - a.y0), (b.x1 - b.x0) * (b.y1 - b.y0))
                if small > 0 and area / small > 0.20:
                    over.append("OVERLAP  '%s' <-> '%s' (%.0f%%)" % (sa[:44], sb[:44], 100 * area / small))
    crowd = []
    for i, ax in enumerate(fig.axes):
        yl = [t for t in ax.get_yticklabels() if (t.get_text() or "").strip()]
        if len(yl) >= 12:
            bbs = sorted([(t.get_window_extent(renderer=r).y0, t.get_window_extent(renderer=r).y1)
                          for t in yl])
            gaps = [bbs[k + 1][0] - bbs[k][1] for k in range(len(bbs) - 1)]
            if gaps and min(gaps) < 0.5:
                crowd.append("CROWD    ax%d %d yticklabels min_gap=%.2fpx" % (i, len(yl), min(gaps)))
    # legibility metrics: smallest font actually drawn, and the tightest vertical
    # gap between consecutive y-tick labels of any panel (print size = canvas size)
    fsz = sorted({round(float(t.get_fontsize()), 2) for s, t, bb in items if s})
    gaps_px = []
    for ax in fig.findobj(matplotlib.axes.Axes):
        yl = [t for t in ax.get_yticklabels() if (t.get_text() or "").strip()]
        if len(yl) >= 8:
            bbs = sorted([(t.get_window_extent(renderer=r).y0, t.get_window_extent(renderer=r).y1)
                          for t in yl])
            g = [bbs[k + 1][0] - bbs[k][1] for k in range(len(bbs) - 1)]
            if g:
                gaps_px.append((min(g), len(yl)))
    tight = min(gaps_px) if gaps_px else (float("nan"), 0)
    leg = ("min_font_pt=%.1f min_ytick_gap_px=%.2f (n=%d)" % (fsz[0], tight[0], tight[1]))
    letters = sorted({s for s, t, bb in items if LETTER.match(s)} |
                     {LETTERPRE.match(s).group(1) for s, t, bb in items if LETTERPRE.match(s)})
    keep = [x for x in letters if x in "ABCDEFG"]
    REPORT.append("=== %s | %s | canvas=%dx%d clipped=%d overlap=%d crowd=%d letters=%s | %s | %s"
                  % (fid, tag, W, H, len(outside), len(over), len(crowd), "".join(keep) or "NONE", leg,
                     " / ".join(t[:40] for t in titles[:2])))
    for x in outside + over + crowd:
        REPORT.append("    " + x)
    plt.close(fig)


def main():
    for fn in FIGLIST:
        _o = plt.Figure.savefig
        plt.Figure.savefig = lambda self, *a, **k: audit(self, fn)
        try:
            runpy.run_path(os.path.join(SC, fn), run_name="__main__")
        except Exception as e:
            REPORT.append("!! %s failed: %s: %s" % (fn, type(e).__name__, e))
        plt.Figure.savefig = _o
        plt.close("all")
    txt = "\n".join(REPORT)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
