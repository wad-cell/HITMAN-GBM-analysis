# Environment note

- `python_version.txt` / `python_environment_freeze.txt`: frozen Python environment (3.11.x) used for the final
  figure assets; the freeze lists the exact matplotlib / numpy / pandas / scipy / Pillow / pypdf / infercnvpy versions.
- `software_versions_bulk_rnaseq.tsv`, `software_versions_microarray_r.tsv`, `software_versions_R_sessionInfo.txt`:
  versions recorded for the analysis stages (Python and R).

## Figure rendering note

The figure masters were produced by the frozen figure scripts. During release preparation the figure scripts were
re-executed as a reproducibility check and the per-panel plotted values were compared with the values recorded when
the masters were generated:

- Fig2 and FigS9 reproduce byte-identically / are identical in plotted values;
- for the remaining figures the plotted-value sets are identical, while minor raster layout differences (caption band
  and axis-tick spacing) may appear if a different matplotlib patch version is used.

`scientific_recalculation` is NO for all panels: no analysis, statistic or scale is recomputed when regenerating the
figures. Any layout difference is cosmetic only and does not affect the plotted values, the panel letters or the
scientific content.
