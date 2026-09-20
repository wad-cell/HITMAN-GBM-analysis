# Code and Results Release Package — HITMAN factorial-specificity / comparative-electrotherapy study

Release package assembled for submission. The reported scientific content is fixed at the versions described in the manuscript; this package contains **no new analyses**.
All items below are either (a) retained analysis code, (b) fixed processed result tables, or (c) environment manifests.
No analysis was re-run, re-fitted or re-derived during release preparation. Two layout-only operations were performed on the figure masters (Fig6 2x2 re-composition; FigS9 render-only re-export) with the plotted values verified unchanged, and the figure scripts were re-executed as a reproducibility check into a separate scratch directory without modifying the released figures.

Source workspace root is referred to as `${PROJECT_ROOT}` throughout. All hard-coded local absolute paths in the
retained scripts were replaced with the placeholders `${PROJECT_ROOT}` / `${USER_HOME}`; local usernames were removed.
No credential, token, password, or API key was present in the retained scripts (verified by pattern scan; 0 hits).

---

## 1. Package structure

```
CODE_RELEASE/
├── README.md                          # this file
├── scripts_python/                    # retained Python analysis code (pipeline order by stage prefix)
├── scripts_R/                         # retained R analysis code
├── environment/                       # software/version manifests
├── genesets/                          # fixed gene-set definitions actually used
├── processed_tables/                  # fixed processed, non-identifiable result tables (70 files)
├── qc_figures/                        # QC diagnostic plots retained for provenance (not manuscript figures)
└── docs/
    └── vendor/_infercnv_src.py        # reference source excerpt kept for CNV-QC provenance
```
The itemised inventory (path, category, status, SHA-256, source path, notes) is in
`CODE_RELEASE_MANIFEST.tsv` (435 data rows; 231 shipped files INCLUDED, 25 EXCLUDED).

---

## 2. Data sources (public only)

| Layer | Data | Accession |
|---|---|---|
| Bulk RNA-seq, primary | GSE324656 | GEO |
| Microarray, cell-line electrotherapy | E-MTAB-9043 | ArrayExpress |
| Single-cell (10x + Smart-seq2) | GSE131928 (Neftel), Couturier et al. layer A | GEO |
| Human-cohort context | TCGA-GBM (GDC), CGGA | GDC / CGGA |

No patient-identifiable data are included. All processed tables derive from public datasets.

---

## 3. Environment

| File | Content |
|---|---|
| `environment/software_versions_bulk_rnaseq.tsv` | package versions for the bulk RNA-seq layer |
| `environment/software_versions_microarray_r.tsv` | package versions for the microarray (R/limma) layer |
| `environment/software_versions_R_sessionInfo.txt` | R `sessionInfo()` for the DESeq2 engine-verification layer |
| `environment/python_environment_freeze.txt` | `pip freeze` snapshot of the packaging-time Python environment (105 entries) |
| `environment/python_version.txt` | Python interpreter version |

Analyses were executed on Windows with Python 3.11.x and R; the release package is not containerised.
No external R environment is bundled — R-side scripts require a working R installation matching
`software_versions_R_sessionInfo.txt`.

---

## 4. Reproduction order

Scripts are retained under their original names; the pipeline order follows the stage prefixes.
`${PROJECT_ROOT}` must be set to a local working directory before execution.

**Environment and software versions**
`scripts_R/versions.R`, `scripts_python/b3_download_libs.py`

**B1/B2 — sample-sheet / metadata audit**
Fixed outputs (e.g. `B1_MATRIX_AUDIT.tsv`, `B2_ANALYSIS_MANIFEST.tsv`, `B2_CONTRAST_DICTIONARY.tsv`) are included,
but their **producer scripts are not present in the workspace** (see §6). Treated as fixed inputs downstream.

**B3 — bulk RNA-seq primary analysis**
`b3_step01_load_verify.py` → `b3_step02_run.py` → `b3_step03_sensitivity.py` → `b3_step04_diag.py` →
(**step05 MISSING**) → `b3_step06_arch.py` → `b3_step07_gsea.py` → `b3_step08_gsea_summary.py` →
`b3_step09_evidence.py` → `b3_step10_gsea_points.py`; helper `b3_gene_map2.py`.

**B35 — differential-engine verification**
`scripts_R/B35_run_deseq2.R`, `scripts_R/B35_engine_verification.R`, `B35_compare_engines.py`,
`B35_final_annot.py`, `B35_deepdive.py`.

**C1 — cell-line microarray layer**
`scripts_R/C1_preprocess.R` → `C1_preprocess2.R` → `C1_preprocess3.R` → `C1_preprocess4.R` →
`C1_limma.R` → `C1_gsea.R` → `C1_dir.R` + `C1_dir_rrho.R`; `scripts_R/C1_versions.R`.

**C2 — single-cell layer**
`c2_layerb_prep.py`, `c2_layerb_run.py`, `scRNA_lib.py`, `null_batch.py`, `run_layerb_batch.py`,
then `c2_11_stepA.py` → `c2_11_stepB.py` → `c2_11_stepB2.py` → `c2_11_stepC.py`.

**C3 — human-cohort contextualisation (TCGA/CGGA)**
`prep_matrices.py`, `prep_sets.py`, `score_ssgsea.py`, `download_tcga.py`, then
`c3_main_analysis.py`, `c3_robust_part1.py`, `c3_robust_genegene.py`, `c3_rebuild_primary.py`,
`c3_rebuild_robust.py`, `c3_finalize_tsv.py`, `c3_figures.py`.

**Intermediates not shipped.** Large intermediate objects are deliberately excluded (e.g.
`C1_raw_combined.rds` 18.0 MB, `C1_eset_rma.rds` 8.1 MB, `C1_eset_rma_annotated.rds` 8.5 MB,
`B35_R_C1..C5.rds` ≈1.3–1.4 MB each, `c2_11_pool_csc.npz` 93.9 MB, single-cell matrices up to ~1 GB).
They are regenerable from the public sources via the C1/C2 scripts, or available from the authors on request.

---

## 5. Gene sets

| File | Role |
|---|---|
| `genesets/C1_directional_sets.gmt` | fixed directional gene sets used in the C1 microarray pipeline |
| `genesets/frozen_signature_genes.tsv` | fixed gene sets used for single-cell module scoring |

External MSigDB/GO libraries cached during analysis (HALLMARK, REACTOME, GOBP) are **not redistributed**
(licensing); they must be obtained from the original provider and re-cached under the same filenames.
Recorded SHA-256 of the cached copies used:
`b3_gslib_HALLMARK.gmt` 37dec27de2228508…, `b3_gslib_REACTOME.gmt` c3ea9df73ad5a442…,
`b3_gslib_GOBP.gmt` 5318d0920aec9a34… (full hashes in `CODE_RELEASE_MANIFEST.tsv`).

---

## 6. Items declared MISSING (not re-created, not invented)

1. `temp/b3_step05_*.py` — pipeline numbering gap (step04 → step06 in the retained sequence).
2. Producer scripts for 41 fixed result tables (e.g. `B1_MATRIX_AUDIT.tsv`, `B2_ANALYSIS_MANIFEST.tsv`,
   `B35_C4/C5_ENGINE_COMPARISON.tsv`, `B35_R_MODEL_MATRIX.csv`, `C1_*_vs_CTRL.tsv`, `C1_PATHWAY_CONCORDANCE.tsv`,
   `C2_12_MASTER_EVIDENCE.tsv`, `C2_11_*`, `C2_*layerA*`, `C4_MASTER_EVIDENCE.tsv`, `C4_NEGATIVE_RESULTS.tsv`,
   `GSE324656_*.tsv`, and others — full list in the manifest with `status = MISSING`).
   The tables themselves **are** included; the generating code was not retained in the workspace.
   No substitute script was written, because any re-written script would not be guaranteed to reproduce the
   frozen numbers exactly.
3. Figure-assembly code for Fig1–Fig5 — no assembly script is retained in the workspace.
4. Repository DOI (GitHub/Zenodo) — none assigned or reserved; must not be fabricated.

---

## 7. Excluded files (25) and why

Excluded categories recorded in the manifest: environment probes (`b3_probe_*.py`, `b3_smoke*.py`, `b3_netprobe.py`,
`test_read.R`), debug/diagnostic probes (`b3_debug_*.py`, `b3_diag_design.py`, `C1_diag*.R`),
superseded copies (`b3_gene_map.py`, `B35_engine_verification_20260907_224843_132.R`),
internal working/report builders (`b3_rebuild_summary.py`, `B35_build_deliverables.py`, `_show_anchors.py`),
and internal manuscript-formatting tooling — none contain analysis logic required
to reproduce released results.

---

## 8. Sanitisation record

| Action | Count |
|---|---|
| scripts with `${PROJECT_ROOT}` replacement | 24 |
| local-username occurrences removed | 0 residual (verified) |
| credential/token/password/API-key patterns | 0 found, 0 residual |
| residual identifier scan over the released package (local username, local user-profile paths, workspace user ID, vendor name, e-mail, bearer tokens) | 0 hits |
| internal pipeline-stage vocabulary removed at assembly (`Stage <label>` in script headers, README and one released table annotation; drafting notes) | 19 occurrences in 12 files, 0 residual (verified) |

---

## 9. Licence and reuse

Released under the **MIT License** — see `LICENSE` in the package root. The license covers the code and the
processed, non-identifiable result tables. Copyright line: `Copyright (c) 2026 The HITMAN study authors`
(confirm the final copyright holder before publication).
The external gene-set libraries in §5 remain under their own terms.


---

## 10. Package-hygiene changes applied at assembly

Internal project-workflow wording was removed from this release: the pipeline-stage label prefix ("Stage ...") was
deleted from script headers/comments, from this README and from the annotation column of one released table, and
"frozen ..." wording in descriptive text was changed to "fixed ...". Script logic, parameters, inputs and result
values were not modified; all Python scripts were re-compiled after editing, and the SHA-256 values in
`CODE_RELEASE_MANIFEST.tsv` were refreshed for every touched file. A full before/after record is provided in
`CODE_RELEASE_PATCH_LOG.tsv` (delivered with this package, not inside it). Filenames that carry the internal layer
numbering (`b3_*`, `B35_*`, `C1_*`, `c2_11_*`, `c3_*`, `frozen_signature_genes.tsv`) were deliberately left
unchanged so that the package stays traceable to the analysis code and to the manifest.

---

## 11. Figure provenance

`FIGURE_PROVENANCE.tsv` (package root) lists, for every main and supplementary figure panel:
`figure`, `panel`, `source_frozen_table`, `plotting_script`, `output_file`, `scientific_recalculation`, `QC_status`.
`scientific_recalculation` is **NO** for all 44 panels: every released figure panel is a rendering of the frozen
result tables, produced without re-analysis or re-statistics.

## 12. Reproduction

See `RUN_REPRODUCTION.md`. All figure scripts read only the frozen tables in `processed_tables/` and
`frozen_intermediates/` and write `figures/` and `figures/supplementary/` in place.

## 13. Environment notes

See `environment/ENVIRONMENT_NOTE.md`; the environment freeze is `environment/python_environment_freeze.txt`.
The released masters are the assets produced by the frozen figure scripts; the regenerated plotted-value logs in
`provenance/plotted_values/` reproduce all plotted values.

## 14. Repository / archive preparation

- Layout prepared for a GitHub repository mirrored to Zenodo: `README.md`, `LICENSE`, `CITATION.cff`,
  `.zenodo.json`, `manifest/` (checksums), `docs/`.
- Repository: https://github.com/wad-cell/HITMAN-GBM-analysis (private until publication).
- No DOI is asserted in this package: repository DOI is assigned at release time.
- Author list / affiliation in `CITATION.cff` and `.zenodo.json` are final (Jing Deng, Ying Wang, Jiateng Zeng,
  Lianghong Yu, Hongliang Ge — Department of Neurosurgery, Neurosurgery Research Institute, The First Affiliated
  Hospital, Fujian Medical University, Fuzhou 350005, Fujian, China). ORCID, GitHub username(s), release date and
  the final copyright holder remain to be supplied by the authors at publication.
- Manuscript statement (Data Availability), compliant wording template — no figure-history language is required:

  > "All processed, non-identifiable result tables, analysis code and figure-generation scripts are provided in the
  > release package / repository accompanying this manuscript (`[repository URL]`, `[DOI]`)."

## 15. Release statements

- No figure panel was re-analysed, re-statisticised or re-interpreted for this release.
- Panel letters and panel content were verified against the plotting scripts, the rendered assets and the manuscript
  legends (main figures) or the supplementary legend enumeration.
- The manuscript scientific content was not modified.
