# Three-GEO bulk-expression comparison

## Current project type

The executed project is a **bulk transcriptomic study using three human kidney-transplant microarray
series**. It is not a single-cell or spatial analysis. Scripts `10`–`13` are future gated utilities;
they are syntax-checked by the main pipeline but are not used to generate the current results because
no eligible single-cell count matrix has been compiled.

## Why these cohorts are compared but not merged

| GEO series | Public design | Platform | Role | Within-study comparison |
|---|---|---|---|---|
| GSE50084 | 115 profiles: 61 biopsies and 54 blood samples | GPL6244 | DSA-state discovery; primary analysis uses 61 biopsies | DSA+/rejection vs DSA+/non-rejection; DSA+/non-rejection vs DSA-/non-rejection |
| GSE36059 | 411 kidney-biopsy profiles | GPL570 | independent replication and disease specificity | ABMR vs non-rejecting; ABMR vs TCMR |
| GSE275824 | 63 serial biopsies from 22 randomized patients | GPL13667 | randomized molecular treatment response | felzartamab-minus-placebo change from baseline at weeks 24 and 52 |

The raw expression matrices are never concatenated: the arrays, sampling designs and biological
questions differ. Each study is annotated and modeled independently. Cross-study comparison uses:

1. Same-direction standardized module effects.
2. Within-study effect estimates and confidence intervals.
3. Hallmark GSEA normalized enrichment scores (NES) and BH FDR.
4. Explicit ABMR-versus-TCMR specificity in GSE36059.
5. Patient-blocked difference-in-differences in GSE275824.

Integrated outputs are `results/atlas/module_effects.csv`,
`results/figures/selected_hallmark_gsea_atlas.csv`, and
`results/tables/Table_2_primary_estimands.csv`.

## Scripts responsible for the datasets

- `scripts/build_public_data.py`: download, checksum verification and compilation of GSE50084.
- `scripts/build_bulk_cohorts.py`: download, verification and compilation of GSE36059 and GSE275824.
- `geo_utils.py`: explicit sample-label parsers and data-contract checks.
- `scripts/02b_annotate_bioc.R`: platform-specific probe-to-gene annotation.
- `scripts/03_bulk_de_gsea.R`: GSE50084 differential expression and GSEA.
- `scripts/04_external_validation.R`: GSE36059 replication and specificity GSEA.
- `scripts/05_treatment_response.R`: GSE275824 patient-blocked longitudinal DE and GSEA.
- `scripts/06_published_score_trajectories.R`: mixed models of published response scores.
- `scripts/08_build_robustness_atlas.R`: standardized integration without pooling expression.
- `scripts/15_make_dataset_result_figures.R` and `16_make_estimand_figure.R`: executed-data plots.

## Adding another GEO bulk dataset

Adding an accession is not only adding a URL. The required reviewed changes are:

1. Add accession, platform, expected sample count, URL and pinned SHA-256 to `configs/datasets.json`.
2. Add a cohort-specific phenotype mapping in `geo_utils.py`; unknown labels must stop execution.
3. Add the acquisition command to the `data` target in `Makefile`.
4. Add platform annotation to the `annotate` target.
5. Add a dataset-specific R analysis with declared case, control, covariates and biological unit.
6. Add effect/GSEA outputs to `scripts/08_build_robustness_atlas.R` only when the estimand is comparable.
7. Add sample-count and artifact assertions to `scripts/14_validate_pipeline_outputs.py`.
8. Add tests for label preservation and exclusion rules.

A study should not be pooled merely to increase sample size. It must contribute a declared role:
discovery, independent replication, rejection specificity, treatment response, or cell localization.
