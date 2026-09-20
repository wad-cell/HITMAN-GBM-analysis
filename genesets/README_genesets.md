# Gene-set definitions — provenance

| File | Content | Role |
|---|---|---|
| `C1_directional_sets.gmt` | fixed directional (up/down) gene sets | C1 microarray layer directional testing |
| `frozen_signature_genes.tsv` | frozen signature gene lists | Layer B single-cell module scoring |

## External libraries (NOT redistributed)

The following external libraries were downloaded during analysis and cached locally. They are **not** included
in this package for licensing reasons. Obtain them from the original provider and cache them under the same
filenames in the analysis working directory:

| Cached filename | Size (bytes) | SHA-256 of the copy used |
|---|---|---|
| `b3_gslib_HALLMARK.gmt` | 44286 | `37dec27de2228508… (full hash in CODE_RELEASE_MANIFEST.tsv)` |
| `b3_gslib_REACTOME.gmt` | 752852 | `c3ea9df73ad5a442… (full hash in CODE_RELEASE_MANIFEST.tsv)` |
| `b3_gslib_GOBP.gmt` | 1509840 | `5318d0920aec9a34… (full hash in CODE_RELEASE_MANIFEST.tsv)` |

## Notes

- Gene-set membership was fixed before the release; no gene set was modified during packaging.
- The `.gmt` files included here are the definitions actually consumed by the pipelines listed in the
  package README; no reconstruction from memory was performed.
