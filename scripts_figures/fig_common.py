"""Figure reproduction helpers for the released figure set (plotting only).

Reads FROZEN result tables only. No statistics are recomputed here.
Frozen sources (resolved relative to PROJECT_ROOT, see the path block below):
  processed_tables/*.tsv           frozen result tables
  frozen_intermediates/*.tsv       frozen analysis intermediates
"""
import os, csv, json, datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Path resolution - no machine-specific absolute paths anywhere in this package.
# PROJECT_ROOT may be exported explicitly; otherwise it is derived from the location
# of this file (the package layout is fixed):
#   <PROJECT_ROOT>/scripts_figures/<this file>
#   <PROJECT_ROOT>/processed_tables/                frozen result tables
#   <PROJECT_ROOT>/frozen_intermediates/            frozen analysis intermediates
#   <PROJECT_ROOT>/figures/[supplementary/]         figure masters
#   <PROJECT_ROOT>/provenance/plotted_values/       plotted-value logs
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
BASE = PROJECT_ROOT
PT   = os.path.join(BASE, "processed_tables")          # frozen result tables
TP   = os.path.join(BASE, "frozen_intermediates")      # frozen analysis intermediates
C3D  = TP
FIGD = os.path.join(BASE, "figures")                   # figure masters
SUPD = os.path.join(FIGD, "supplementary")
SCRD = os.path.join(BASE, "scripts_figures")
LOGD = os.path.join(BASE, "provenance", "plotted_values", "release")
for d in (FIGD, SUPD, LOGD):
    os.makedirs(d, exist_ok=True)

DPI_TIF = 400
DPI_PNG = 300

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def rd(path):
    """Read a frozen TSV into a list of dicts."""
    with open(path, encoding="utf-8", errors="replace") as f:
        return [r for r in csv.DictReader(f, delimiter="\t")]


def num(x):
    if x is None:
        return None
    s = str(x).strip().replace("*", "").replace(",", "")
    if s in ("", "NA", "NaN", "nan", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def panel_label(ax, letter, dx=-0.085, dy=1.035):
    """Draw a panel letter and register it so that the layout pass can fold it
    into the axes title (plotting-only; no content change)."""
    t = ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom", ha="left")
    ax._panel_letter = letter
    ax._panel_letter_text = t
    return t


def text_panel(ax, lines, fontsize=7.5, x=0.0, y=1.0, dy=0.075, wrap=None):
    ax.set_axis_off()
    for i, ln in enumerate(lines):
        ax.text(x, y - i * dy, ln, transform=ax.transAxes, fontsize=fontsize,
                va="top", ha="left", family="DejaVu Sans")


# ---------------------------------------------------------------------------
# Layout-only post-processing.  These helpers never touch data: they only wrap,
# shift or slightly resize already-plotted text so that nothing is clipped at
# the canvas edge and no two annotations collide.
# ---------------------------------------------------------------------------
import textwrap as _tw
from matplotlib.text import Text as _Text


def _wrap(s, width_in, fs):
    n = max(16, int(width_in / (0.55 * fs / 72.0)))
    return _tw.wrap(str(s), n, break_long_words=False, break_on_hyphens=False) or [str(s)]


def wrap_text(s, width_in, fs):
    """Public wrapper: break a string so it fits `width_in` inches at `fs` pt."""
    return _wrap(s, width_in, fs)


def _renderer(fig):
    fig.canvas.draw()
    try:
        return fig.canvas.get_renderer()
    except AttributeError:  # pragma: no cover
        from matplotlib.backend_bases import FigureCanvasBase
        return FigureCanvasBase(fig).get_renderer()


def reflow_rows(ax, rows, fontsize=6.4, top=1.0, x=0.0, lead=1.30, row_gap=0.5,
                width_in=None):
    """Draw a list of single-line rows as wrapped rows inside `ax`.

    Rendering-only: the strings are unchanged, they are only broken across
    lines so that nothing leaves the figure canvas.  Returns the line count.
    """
    fig = ax.figure
    if width_in is None:
        width_in = ax.get_position().width * fig.get_size_inches()[0] - 0.02
    h_in = ax.get_position().height * fig.get_size_inches()[1]
    lh_ax = (fontsize * lead / 72.0) / h_in
    y = top
    n_lines = 0
    for s in rows:
        parts = _wrap(s, width_in, fontsize)
        ax.text(x, y, "\n".join(parts), transform=ax.transAxes, fontsize=fontsize,
                va="top", ha="left", linespacing=lead)
        y -= lh_ax * (len(parts) + row_gap)
        n_lines += len(parts)
    return n_lines


def fit_rows_in_ax(ax, rows, fs_max=6.4, fs_min=5.2, top=0.97, lead=1.28,
                   row_gap=0.35, width_in=None, pad=0.02):
    """Draw wrapped rows inside `ax`, auto-shrinking the font so the whole block
    stays inside the axes box.  Rendering-only: strings are never changed.
    Returns (fontsize_used, row_count)."""
    fig = ax.figure
    if width_in is None:
        width_in = ax.get_position().width * fig.get_size_inches()[0] - pad
    h_in = ax.get_position().height * fig.get_size_inches()[1]
    fs = fs_max
    while True:
        n = sum(len(_wrap(s, width_in, fs)) for s in rows)
        need = (n + row_gap * len(rows)) * fs * lead / 72.0
        if need <= h_in * top or fs <= fs_min + 0.01:
            break
        fs = round(fs - 0.2, 2)
    n = reflow_rows(ax, rows, fontsize=fs, top=top, lead=lead, row_gap=row_gap,
                    width_in=width_in)
    return fs, n


def _title_wrap_width(fig, ax, p, min_in=2.6):
    """Usable width (inches) for a left-aligned panel header.

    The header starts at the axes' left edge.  It may extend past the axes
    when the axes is narrow (e.g. a square heatmap inside a tall row), but it
    must stop before the next substantial axes to its right and before the
    canvas edge.  No content is altered - only where the wrap breaks.
    """
    W_in = fig.get_size_inches()[0]
    left_in = p.x0 * W_in
    right_lim = W_in - left_in - 0.04
    for o in fig.axes:
        if o is ax:
            continue
        q = o.get_position()
        if q.width * W_in < 0.5:                # colour-bar / hairline axes
            continue
        if q.x0 > p.x0 + 1e-6 and q.y0 < p.y1 - 1e-6 and q.y1 > p.y0 + 1e-6:
            right_lim = min(right_lim, q.x0 * W_in - left_in - 0.04)
    return max(p.width * W_in - 0.03, min(min_in, right_lim))


def fit_titles(fig, pad=2.5):
    """Fold each panel letter into its axes header (left aligned, wrapped to the
    usable width) and draw it as a free text above the axes, so that letters
    and titles can neither collide with each other nor descend into the panel.
    The original centre title is cleared to avoid a duplicated title."""
    W_in = fig.get_size_inches()[0]
    H_in = fig.get_size_inches()[1]
    n = 0
    for ax in list(fig.axes):
        letter = getattr(ax, "_panel_letter", None)
        if letter is None:
            continue
        lab = getattr(ax, "_panel_letter_text", None)
        ttl_c = ax.get_title()                  # centre title
        ttl_l = ax.get_title(loc="left")        # pre-existing left title
        ttl = ttl_c or ttl_l
        if not ttl:
            continue
        src = ax.title if ttl_c else ax._left_title
        fs = src.get_fontsize()
        p = ax.get_position()
        lines = _wrap(ttl, _title_wrap_width(fig, ax, p), fs)
        ax.set_title("", loc="left")            # never render a real title again
        if ttl_c:
            ax.title.set_text("")
        if lab is not None:
            try:
                lab.remove()
            except Exception:
                pass
        dy = (pad / 72.0) / max(p.height * H_in, 1e-6)
        ax.text(0.0, 1.0 + dy, letter + "   " + "\n".join(lines),
                transform=ax.transAxes, ha="left", va="bottom", fontsize=fs,
                clip_on=False, zorder=6)
        ax._panel_letter_text = None
        n += 1
    return n


def shrink_overflow_texts(fig, floor=5.2):
    """Shrink single-line annotations whose bbox would leave the canvas."""
    W, H = fig.canvas.get_width_height()
    r = _renderer(fig)
    groups = {}
    for t in fig.findobj(_Text):
        s = (t.get_text() or "")
        if not s.strip() or "\n" in s:
            continue
        fs = float(t.get_fontsize())
        if fs <= floor + 0.01 or not t.get_visible():
            continue
        bb = t.get_window_extent(renderer=r)
        if bb.x1 <= W - 2 and bb.x0 >= 2 and bb.y0 >= 2 and bb.y1 <= H - 2:
            continue
        need = 1.0
        if bb.x1 > W - 2:
            need = bb.width / max(1.0, (W - 2) - bb.x0)
        if bb.x0 < 2 or bb.y1 > H - 2 or bb.y0 < 2:
            need = max(need, 1.05)
        key = (round(t.get_position()[0], 2), round(fs, 1), str(t.get_transform()))
        e = groups.setdefault(key, [fs, []])
        e[0] = min(e[0], fs / max(need, 1e-6))
        e[1].append(t)
    k = 0
    for (x, fs0, tr), (fs_new, texts) in groups.items():
        fs_new = max(floor, min(fs0, fs_new))
        if fs_new < fs0 - 0.05:
            for t in texts:
                t.set_fontsize(fs_new)
            k += 1
    return k


def autowrap_overflow_texts(fig, min_len=30, pad=3.0):
    """Wrap any single-line annotation that is wider than the canvas.

    Rendering-only: the string content is unchanged, it is only broken across
    lines so that the whole annotation fits inside the exported canvas.
    """
    W, H = fig.canvas.get_width_height()
    r = _renderer(fig)
    n = 0
    for t in list(fig.findobj(_Text)):
        s = (t.get_text() or "")
        if len(s) < min_len or "\n" in s or not t.get_visible():
            continue
        if abs(float(t.get_rotation()) % 360.0) > 0.5:
            continue
        bb = t.get_window_extent(renderer=r)
        if bb.x1 <= W - pad:
            continue
        avail = max(80.0, (W - pad) - max(bb.x0, pad))
        cpl = max(16, int(len(s) * avail / max(bb.width, 1.0)) - 2)
        t.set_text("\n".join(_tw.wrap(s, cpl, break_long_words=False, break_on_hyphens=False)))
        n += 1
    return n


def clamp_texts(fig, pad=2.0):
    """Last resort: nudge any remaining out-of-canvas text back inside."""
    W, H = fig.canvas.get_width_height()
    r = _renderer(fig)
    moved = 0
    for t in fig.findobj(_Text):
        if not (t.get_text() or "").strip():
            continue
        bb = t.get_window_extent(renderer=r)
        dx = dy = 0.0
        if bb.x0 < pad:
            dx = pad - bb.x0
        if bb.x1 > W - pad:
            dx = (W - pad) - bb.x1
        if bb.y0 < pad:
            dy = pad - bb.y0
        if bb.y1 > H - pad:
            dy = (H - pad) - bb.y1
        if abs(dx) < 0.5 and abs(dy) < 0.5:
            continue
        tr = t.get_transform()
        x, y = t.get_position()
        px0, py0 = tr.transform((x, y))
        px1, _ = tr.transform((x + 1.0, y))
        _, py2 = tr.transform((x, y + 1.0))
        sx, sy = px1 - px0, py2 - py0
        nx = x + dx / sx if abs(sx) > 1e-9 else x
        ny = y + dy / sy if abs(sy) > 1e-9 else y
        t.set_position((nx, ny))
        moved += 1
    return moved


def prune_ticks(fig):
    """Drop major-tick positions that fall outside an axis' view limits.

    Matplotlib may keep locator ticks that sit marginally outside the limits
    (e.g. a '3' on an axis ending at 2.40); their labels are then painted off
    the canvas. Only auto-formatted (non-categorical) axes are touched, tick
    positions inside the limits are kept verbatim. Rendering-only fix.
    """
    from matplotlib.ticker import FixedLocator, FixedFormatter
    n = 0
    for ax in list(fig.axes):
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            try:
                if isinstance(axis.get_major_formatter(), FixedFormatter):
                    continue
                locs = list(axis.get_majorticklocs())
                keep = [v for v in locs if lim[0] - 1e-9 <= v <= lim[1] + 1e-9]
                if len(keep) == len(locs) or not keep:
                    continue
                axis.set_major_locator(FixedLocator(keep))
                n += 1
            except Exception:
                continue
    return n


def save_fig(fig, name, subdir=""):
    """Write TIFF (>=300 dpi, RGB, LZW) + PDF + PNG masters. Returns dict of paths."""
    prune_ticks(fig)
    fit_titles(fig)
    autowrap_overflow_texts(fig)
    shrink_overflow_texts(fig)
    clamp_texts(fig)
    out = os.path.join(FIGD, subdir) if subdir else FIGD
    os.makedirs(out, exist_ok=True)
    tif = os.path.join(out, name + ".tif")
    pdf = os.path.join(out, name + ".pdf")
    png = os.path.join(out, name + ".png")
    fig.savefig(tif, dpi=DPI_TIF, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    fig.savefig(png, dpi=DPI_PNG, facecolor="white")
    plt.close(fig)
    # enforce 8-bit RGB (drop alpha) and keep 400 dpi tag
    try:
        from PIL import Image
        im = Image.open(tif)
        if im.mode != "RGB":
            bg = Image.new("RGB", im.size, (255, 255, 255))
            im = im.convert("RGBA") if im.mode == "P" else im
            bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
            im = bg
        im.save(tif, format="TIFF", compression="tiff_lzw", dpi=(DPI_TIF, DPI_TIF))
    except Exception as e:  # pragma: no cover
        print("WARN post-process tif:", e)
    try:
        from PIL import Image
        im = Image.open(png)
        if im.mode != "RGB":
            bg = Image.new("RGB", im.size, (255, 255, 255))
            bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
            im = bg
        im.save(png, format="PNG", dpi=(DPI_PNG, DPI_PNG))
    except Exception as e:  # pragma: no cover
        print("WARN post-process png:", e)
    return {"tif": tif, "pdf": pdf, "png": png}


def log_values(fig_name, rows):
    """Append plotted values to the traceability log (metric, source, value)."""
    p = os.path.join(LOGD, fig_name + "_values.tsv")
    new = not os.path.exists(p)
    with open(p, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        if new:
            w.writerow(["figure", "panel", "plotted_metric", "source_file", "plotted_value"])
        for r in rows:
            w.writerow([fig_name] + list(r))
    return p


def star(fdr, p=None):
    v = fdr if fdr is not None else p
    if v is None:
        return ""
    return "***" if v < 0.001 else "**" if v < 0.01 else "*" if v < 0.05 else ""


STAMP = "released figure set - frozen tables only"


def add_caption(fig, lines, fontsize=5.6, width=150, x=0.02):
    """Add a wrapped caption band INSIDE the figure canvas (never clipped).

    Plotting-only helper: the canvas is extended downward by the caption height
    and every existing axes is translated up so that its size in inches is
    unchanged. No data are read, recomputed or restyled.
    """
    import textwrap
    wrapped = []
    for ln in lines:
        wrapped += textwrap.wrap(str(ln), width) or [str(ln)]
    h = fig.get_size_inches()[1]
    w = fig.get_size_inches()[0]
    lh = fontsize * 1.45 / 72.0
    band = lh * (len(wrapped) + 1.8)
    fig.set_size_inches(w, h + band)
    for ax in fig.axes:
        p = ax.get_position()
        ax.set_position([p.x0, (p.y0 * h + band) / (h + band), p.width,
                         p.height * h / (h + band)])
    fig.text(x, (band - lh * 0.9) / (h + band), "\n".join(wrapped),
             fontsize=fontsize, va="top", ha="left")
    return band
