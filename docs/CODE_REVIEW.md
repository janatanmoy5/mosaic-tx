# Script-by-script verification and scientific code review

## Review standard

Every repository script was checked for its declared input, output, biological replicate, statistical
model, failure behavior and current execution status. Python files pass byte-code compilation; R files
pass parser validation. Existing unit tests and `make check` pass on the local Mac. “Syntax verified”
does not mean a data workflow was executed.

## Executed scripts

| Script | Review outcome | Execution evidence |
|---|---|---|
| `build_public_data.py` | Checksum validation, sample-count validation and explicit dataset restriction are appropriate. | Previously executed for GSE50084; compiled files and manifests present. |
| `build_bulk_cohorts.py` | Dataset choices are explicit; delegates validation to tested library functions. Source is overly compressed and should be reformatted before a methods release. | Previously executed for GSE36059 and GSE275824. |
| `02b_annotate_bioc.R` | Uses platform-specific Bioconductor mappings and averages duplicate probes. “First” mapping for multi-mapped probes is a sensitivity limitation. | Executed for all three bulk cohorts; annotation manifests present. |
| `03_bulk_de_gsea.R` | Correct moderated two-group models and BH correction. GSE50084 lacks covariate adjustment and depends on historical labels. | Executed; DE, GSEA, volcanoes and manifest present. |
| `04_external_validation.R` | Correctly tests ABMR versus non-rejection and separately versus TCMR. It does not establish ABMR specificity unless effect direction is evaluated. | Executed; the pilot atlas now enforces direction. |
| `05_treatment_response.R` | Patient blocking and placebo-adjusted difference-in-differences are appropriate. With 22 patients, estimates and intervals matter more than gene-level significance. | Executed at W24 and W52. |
| `06_published_score_trajectories.R` | Patient random intercept and treatment-by-time interactions are appropriate; BH is applied across declared score/time tests. | Executed; six published scores evaluated. |
| `07_score_mechanism_modules.py` | Correct within-cohort z-score construction and coverage audit. It requires the declared Python scientific environment. | Smoke-tested with synthetic expression; superseded for the pilot by script 08. |
| `08_build_robustness_atlas.R` | Fully annotated. Separates platforms, freezes modules, uses patient-blocked longitudinal models and direction-aware evidence grades. | Executed end to end on three cohorts; manifest and tables present. |
| `09_make_pilot_figures.R` | Fully annotated. Reads saved results without recomputing statistics and exports PNG/PDF. | Executed; combined and individual panels present. |

## Syntax-verified but not yet dataset-executed

| Script | Blocking condition and review decision |
|---|---|
| `10_build_scanpy_object.py` | Accepts 10x H5/MTX or H5AD, validates human NovaSeq metadata and preserves counts. E-MTAB-12051 is distributed as a large count CSV, so a chunked CSV-to-sparse compiler is still required. Do not call this workflow complete. |
| `11_sc_qc_cluster.py` | Performs basic Scanpy normalization, HVG, PCA, neighbors, UMAP and Leiden. Dataset-specific doublet, ambient RNA, batch and QC thresholds are not implemented; therefore it remains scaffolding. |
| `12_patient_pseudobulk.py` | Correctly aggregates raw counts and enforces minimum cells. Downstream edgeR/DESeq2 modeling and patient/cell-type design are still required. |
| `13_grouped_nested_ml.py` | Keeps patients together and nests tuning, but lacks external validation, bootstrap intervals and small-fold safeguards. It must not be run for the 22-patient treatment cohort. |

## Library modules and tests

| Component | Review outcome |
|---|---|
| `mosaic_tx/geo.py` | GEO parsing, verified download and GSE50084 compilation are covered by unit tests. |
| `mosaic_tx/geo_multi.py` | GSE36059 group preservation and GSE275824 patient/treatment/time requirements are tested. |
| `mosaic_tx/io.py` | Shared checksum, atomic-download and manifest utilities; no destructive operations. |
| `tests/test_geo.py` | Covers expected compilation invariants. |
| `tests/test_geo_multi.py` | Protects nephrectomy separation and longitudinal identifiers. |

## Defects found by real execution

1. The first pilot run referenced a factor-dependent coefficient name; execution stopped. The code
   now uses the invariant second coefficient in a two-group design.
2. The first figure run assumed a non-existent `SE` column and a Cairo/X11 PDF device. Confidence
   intervals now use the moderated effect/t statistic and the native PDF device.
3. The first evidence grade treated any significant ABMR-versus-TCMR difference as specificity,
   including effects higher in TCMR. Tier A now requires significant positive ABMR enrichment, and
   negative-control modules receive `CONTROL`, never a biological evidence grade.

## Remaining mandatory work

- Reformat and annotate the older compressed acquisition/R scripts.
- Add unit tests for module scoring, evidence-tier direction, and pseudobulk aggregation.
- Implement and test chunked conversion for E-MTAB-12051.
- Freeze dataset-specific single-cell QC and batch-handling specifications before execution.
- Add continuous integration tests using small fixture matrices for scripts 08 and 09.
- Do not enable ML, drug reversal or QSAR targets until their upstream evidence gates pass.
