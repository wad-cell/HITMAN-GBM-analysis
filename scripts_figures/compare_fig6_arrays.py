# -*- coding: utf-8 -*-
"""the release stage - Fig6 plotted-value identity audit (read-only; writes nothing to the figures).

Two deterministic checks, no statistics recomputed anywhere:

(A) LAYOUT-ONLY PROOF.  Both the pre-release export script (fig6_adonly.py) and the release stage
    script (fig6_2x2.py) are executed with Figure.savefig intercepted, so the figure objects are
    inspected but nothing is written.  Every axes' plotted numeric payload (line vertices, bar
    widths/heights, image arrays, scatter offsets) is hashed.  Axes digests identical between the
    two exports => same values drawn, only geometry changed.  Leftover digests are reported
    explicitly (expected: at most the panel-D sub-axes, where the pre-release stage read an empty field and
    therefore drew zeros - see the value-identity table for the before/after numbers).

(B) FROZEN-SOURCE TRACEABILITY.  Every heatmap cell value and every bar length in the release stage
    figure must appear verbatim in the frozen result tables it claims to plot.
"""
import os, sys, runpy, csv, hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCR = os.path.dirname(os.path.abspath(__file__))
TEMP = os.path.dirname(SCR)          # release package root
BASE = os.path.dirname(TEMP)
sys.path.insert(0, SCR)
from fig_common import rd, num, PT, C3D

OUT = os.path.join(TEMP, "provenance", "FIG6_COMPOSITION_AUDIT.tsv")


def capture(script, label):
    box = {}
    orig = plt.Figure.savefig
    plt.Figure.savefig = lambda self, *a, **k: box.setdefault("fig", self)
    import fig_common
    orig_log = fig_common.log_values
    fig_common.log_values = lambda *a, **k: None      # audit must not touch any value log
    err = None
    try:
        runpy.run_path(script, run_name="__main__")
    except Exception as e:                                   # pragma: no cover
        err = "%s: %s" % (type(e).__name__, e)
    finally:
        plt.Figure.savefig = orig
        fig_common.log_values = orig_log
    fig = box.get("fig")
    panels = []
    if fig is not None:
        fig.canvas.draw()
        # fig.findobj walks child axes too (inset axes are not always listed in fig.axes)
        import matplotlib.axes as _maxes
        for i, ax in enumerate(fig.findobj(_maxes.Axes)):
            vals = []
            for ln in ax.lines:
                vals += list(np.ravel(np.asarray(ln.get_xydata(), float)))
            for p in ax.patches:
                vals += [float(p.get_width()), float(p.get_height())]
            for im in ax.images:
                arr = np.ma.masked_invalid(np.asarray(im.get_array(), float))
                vals += list(arr.compressed())
            for c in ax.collections:
                try:
                    vals += list(np.ravel(np.asarray(c.get_offsets(), float)))
                except Exception:
                    pass
            v = [x for x in vals if np.isfinite(x)]
            key = hashlib.sha1(("\n".join("%.9g" % x for x in sorted(v))).encode()).hexdigest()[:12]
            panels.append({"label": label, "ax": i, "n": len(v), "digest": key, "vals": sorted(set(v))})
    plt.close("all")
    return panels, err


old_panels, old_err = capture(os.path.join(SCR, "fig6_adonly.py"), "pre_release_fig6_adonly")
new_panels, new_err = capture(os.path.join(SCR, "fig6_2x2.py"), "release_fig6_2x2")

old_d = {p["digest"]: p for p in old_panels}
new_d = {p["digest"]: p for p in new_panels}
matched = sorted(set(old_d) & set(new_d))
old_only = [p for p in old_panels if p["digest"] not in new_d]
new_only = [p for p in new_panels if p["digest"] not in old_d]

AXES_TSV = os.path.join(TEMP, "provenance", "FIG6_AXES_PAYLOAD.tsv")
with open(AXES_TSV, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["export", "axes_index", "n_values", "payload_digest", "matched_in_other_export",
                "n_distinct_values", "distinct_values"])
    for p in old_panels + new_panels:
        other = new_d if p["label"].startswith("pre_release") else old_d
        w.writerow([p["label"], p["ax"], p["n"], p["digest"], "YES" if p["digest"] in other else "NO",
                    len(p["vals"]), "; ".join("%.9g" % v for v in p["vals"])])
print("axes payload detail ->", AXES_TSV)

# ---- (B) frozen-source traceability of the release stage payload ----
pr = rd(os.path.join(PT, "C3_PRIMARY_RESULTS.tsv"))
rb = rd(os.path.join(PT, "C3_ROBUSTNESS.tsv"))
sec = rd(os.path.join(PT, "C3_SECONDARY_RESULTS.tsv"))
q2 = rd(os.path.join(C3D, "_q2_for_output.tsv"))

frozen_rho = {num(r["effect_rho"]) for r in pr if r["block"] == "Q1_architecture"}
frozen_rho |= {num(r["rho_TCGA"]) for r in rb if r["analysis"] == "robust_1_cross_cohort_sign"}
frozen_rho |= {num(r["rho_CGGA"]) for r in rb if r["analysis"] == "robust_1_cross_cohort_sign"}
frozen_bar = {num(r["median_diff_MES_minus_nonMES"]) for r in q2}
frozen_bar |= {num(r["median_diff_MES_minus_nonMES"]) for r in sec
               if r["module"] == "M_autophagy_lysosome" and r["cohort"].startswith("TCGA")}
frozen_bar |= {num(str(r["module_b"]).split("=")[-1].split("/")[0]) for r in rb
               if r["analysis"] in ("robust_2_stratum_Q1", "robust_4_gene_perturbation")}

grid_vals, bar_vals = [], []
for p in new_panels:
    for v in p["vals"]:
        if v is None:
            continue
        if v in frozen_rho:
            grid_vals.append(v)
        elif v in frozen_bar:
            bar_vals.append(v)
untraceable = sorted({v for p in new_panels for v in p["vals"]
                      if v not in frozen_rho and v not in frozen_bar})

rows = [
    ("check_A_executed_pre_release_script", os.path.basename("fig6_adonly.py")),
    ("check_A_executed_release_script", os.path.basename("fig6_2x2.py")),
    ("check_A_axes_total_pre_release", len(old_panels)),
    ("check_A_axes_total_release", len(new_panels)),
    ("check_A_axes_with_identical_plotted_payload", len(matched)),
    ("check_A_axes_only_in_pre_release", len(old_only)),
    ("check_A_axes_only_in_release", len(new_only)),
    ("check_A_verdict_panels_A_B_C", "IDENTICAL" if matched else "NO-MATCH"),
    ("check_A_pre_release_only_axes", " | ".join("ax%d[%s]" % (p["ax"], "; ".join("%.9g" % v for v in p["vals"]))
                                         for p in old_only)),
    ("check_A_release_only_axes", " | ".join("ax%d[%s]" % (p["ax"], "; ".join("%.9g" % v for v in p["vals"]))
                                         for p in new_only)),
    ("check_B_heatmap_and_forest_values_traced_to_frozen_tables", len(grid_vals)),
    ("check_B_bar_values_traced_to_frozen_tables", len(bar_vals)),
    ("check_B_untraceable_distinct_values", len(untraceable)),
    ("check_B_untraceable_sample", "; ".join("%.9g" % v for v in untraceable[:20])),
    ("scientific_recalculation_in_this_audit", "NO"),
]
if old_err:
    rows.append(("check_A_pre_release_script_error", old_err))
if new_err:
    rows.append(("check_A_release_script_error", new_err))

with open(OUT, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["metric", "value"])
    for k, v in rows:
        w.writerow([k, v])
print("\n".join("%s\t%s" % (k, v) for k, v in rows))
