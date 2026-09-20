# Figure reproduction instructions

The released figures are renders of frozen result tables. Re-running the figure scripts performs **no analysis**:
each script only reads a frozen table and draws the panels.

## 1. Environment

- Python 3.11 (see `environment/python_version.txt`).
- Matplotlib / NumPy / pandas / SciPy / Pillow / pypdf / infercnvpy versions: `environment/python_environment_freeze.txt`.
- R is not required for the figures (R code is retained for the analysis stages only).

## 2. Run order (from the package root)

```
python scripts_figures/fig1.py
python scripts_figures/fig2_ac_only.py
python scripts_figures/fig3.py
python scripts_figures/fig4.py
python scripts_figures/fig5.py
python scripts_figures/fig6_adonly.py          # Fig6 panels A-D
python scripts_figures/fig6_2x2.py             # Fig6 2x2 re-composition (master)
python scripts_figures/figS.py                 # Supplementary FigS1-S9 panel build
python scripts_figures/figS9_render.py         # released FigS9 render
```

Outputs are written in place to `figures/` (main figures) and `figures/supplementary/` (supplementary figures),
as TIFF (400 dpi, RGB), PNG and PDF. Per-panel plotted values are logged by the scripts
(see `provenance/plotted_values/`).

To keep the released files untouched, copy the package to a scratch directory and run there.

## 3. Verification

- `provenance/plotted_values/*_values.tsv` — per-panel values as emitted by the released scripts.
- `provenance/qc_assets/FIGURE_VALUE_QC.tsv` — comparison of the regenerated logs against the stage-recorded logs.
- `provenance/qc_assets/FIGURE_TECHNICAL_QC.tsv` — dpi / size / colour-mode / integrity checks of every asset.
- `provenance/qc_assets/FIGURE_FIG6_2X2_VALUE_IDENTITY.tsv` — panel-by-panel value identity for the Fig6
  2x2 re-composition against the panel render.
- `manifest/SHA256SUMS.tsv` — checksums of every file in the package.
