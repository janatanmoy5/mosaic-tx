# MOSAIC-Tx

**M**ulticohort **O**mics **S**tudy of **A**ntibody-associated **I**njury and its **C**ourse after
**T**ransplantation.

MOSAIC-Tx is a public-data reanalysis of kidney-allograft antibody-mediated rejection (ABMR). It
asks one biological question:

> Which ABMR expression programs are reproducible, rejection-specific, cell-localizable, and sufficiently robust across reference definitions and treatment states to justify downstream drug-reversal analysis?

This repository does not claim patient-specific treatment prediction, clinical validation, wet-lab
validation, causal target discovery, or validated drug repurposing. It compiles published,
de-identified datasets and tests prespecified contrasts with reproducible statistics.

## Pipeline outputs

All generated plots, statistical tables, sample labels, manifests, logs, and checksums are written under `results/`. No separate manuscript folder is created.

## Why this question is supported by available data

GSE50084 contains DSA-negative, DSA-positive/no-rejection, and DSA-positive/rejection samples.
GSE36059 supplies a larger independent ABMR/TCMR/non-rejection biopsy cohort. GSE275824 contains
longitudinal randomized placebo/felzartamab biopsies with published ABMR, IFN-gamma, NK-cell and
endothelial molecular scores. E-MTAB-12051 and GSE145927 can localize replicated programs to cell
types, but their modest patient counts preclude cell-level machine-learning claims.

## Verified scope

- [x] Checksum-verified acquisition of three public bulk cohorts
- [x] 589 expression profiles compiled in total
- [x] Platform-specific Bioconductor probe-to-gene annotation
- [x] Moderated differential expression and Hallmark GSEA
- [x] Independent ABMR-versus-non-rejection replication
- [x] ABMR-versus-TCMR specificity analysis
- [x] Patient-blocked randomized treatment-response analysis
- [x] Published-score trajectory analysis
- [ ] Single-cell patient-level pseudobulk localization
- [x] Cross-cohort effect-size synthesis with explicit TCMR specificity testing
- [ ] Externally validated patient-level biomarker model
- [ ] Evidence-qualified drug-signature reversal
- [ ] Target-specific QSAR (only if the activity-data gate passes)

The new single-cell branch is restricted to human NovaSeq 6000 data and builds a canonical Scanpy
`AnnData` object with optional MuData export. Stable allograft, non-rejection injury and native healthy
kidney are separate reference domains. Spatial transcriptomics remains exploratory. Drug reversal,
QSAR, docking and molecular dynamics are downstream gates—not assumed deliverables.

## Start here

1. [`docs/ABSTRACT.md`](docs/ABSTRACT.md) — structured project abstract
2. [`docs/MANUSCRIPT_DRAFT.md`](docs/MANUSCRIPT_DRAFT.md) — Introduction, Methods, Results, Discussion, and figure legends
3. [`docs/STUDY_BLUEPRINT.md`](docs/STUDY_BLUEPRINT.md) — biological question, estimands, models, figures, and decision gates
4. [`docs/DATASET_COMPARISON.md`](docs/DATASET_COMPARISON.md) — cohort scope and adding another GEO series
5. [`docs/ANALYSIS_SPECIFICATION.md`](docs/ANALYSIS_SPECIFICATION.md) — frozen contrasts and statistics
6. [`docs/RESULTS.md`](docs/RESULTS.md) — executed findings and claim boundaries
7. [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) — claims the data cannot support
8. [`docs/PIPELINE_EXECUTION.md`](docs/PIPELINE_EXECUTION.md) — one-command execution and outputs

## Run the canonical manuscript pipeline

From the repository root, run exactly one command:

```bash
bash scripts/run_three_dataset_pipeline.sh
```

The script executes acquisition and checksum verification, platform annotation, all declared statistical models, GSEA, module scoring, inferential figures, sample-level verification, tests, artifact validation, and manuscript packaging in order. It stops immediately if any stage fails.

Primary plots are written to `results/figures/`; audit tables are written to `results/tables/`. `Table_S1_sample_labels.csv` contains all 589 public profiles with accession, original title, analysis label, inclusion status, treatment/timepoint where applicable, and the reason for inclusion or exclusion.

Equivalent stage-by-stage execution:

```bash
make data
make annotate
make bulk
make scores
make atlas
make check
```

The tracked analysis results were generated from the three public cohorts on the local Mac. Single-cell,
ML, drug reversal, QSAR, docking and MD remain gated and are not represented as completed.

## Cohort and biospecimen contract

| Cohort | Profiles used | Primary specimen | Role | Prespecified comparisons |
|---|---:|---|---|---|
| GSE50084 | 61 biopsies; 54 blood profiles retained separately | Kidney-allograft biopsy | DSA-conditioned discovery | DSA+/rejection vs DSA+/no rejection; DSA+/no rejection vs DSA-/no rejection |
| GSE36059 | 381 primary-analysis profiles | Kidney-allograft biopsy | independent replication and specificity | ABMR vs non-rejecting; ABMR vs TCMR |
| GSE275824 | 63 biopsies from 22 randomized patients | Serial kidney-allograft biopsy | randomized molecular response | felzartamab-minus-placebo change from baseline at weeks 24 and 52 |

Blood and biopsy expression are never pooled. GSE50084 blood is retained in the audit table but is
not used in the primary tissue analysis. GSE36059 nephrectomy and mixed-rejection profiles remain in
metadata but are excluded from primary contrasts. Only GSE275824 supports connected trajectories.

## Pipeline execution order

| Order | Program | Purpose | Principal output |
|---:|---|---|---|
| 1 | `make data` | download, checksum and compile the GEO series | processed sample tables and probe matrices |
| 2 | `scripts/02b_annotate_bioc.R` | platform-specific probe-to-gene mapping | gene matrices and annotation manifests |
| 3 | `scripts/03_bulk_de_gsea.R` | GSE50084 moderated DE and Hallmark GSEA | discovery DE/GSEA tables |
| 4 | `scripts/04_external_validation.R` | replication and ABMR-versus-TCMR specificity | validation DE/GSEA tables |
| 5 | `scripts/05_treatment_response.R` | patient-blocked treatment-by-time contrasts | week-24 and week-52 effects |
| 6 | `scripts/06_published_score_trajectories.R` | mixed-effects analysis of deposited MMDx programs | interaction estimates and manifest |
| 7 | `scripts/08_build_robustness_atlas.R` | harmonize module effects across cohorts | `results/atlas/` |
| 8 | `scripts/15_make_dataset_result_figures.R` | cohort, distribution, forest and MMDx displays | Figure 1 and audit CSVs |
| 9 | `scripts/16_make_estimand_figure.R` | inferential synthesis, GSEA and platform-specific PCA | Figure 2 and PCA coordinates |
| 10 | `scripts/17_build_audit_tables.py` | build sample, cohort, estimand and module tables | `results/tables/` |
| 11 | `make check` plus output validator | tests, syntax checks and artifact contracts | fail-fast validation |
| 12 | `shasum` | hash principal outputs | `results/output_checksums.sha256` |

## Statistical interpretation

- Differential expression uses limma moderated models with Benjamini-Hochberg FDR within each contrast.
- Hallmark GSEA uses preranked fgsea; direction follows the printed numerator-minus-comparator label.
- Mechanism modules are prespecified within-cohort mean gene-z scores and cannot establish cell origin.
- GSE275824 inference uses a patient-blocked treatment-by-time interaction.
- Figure 1, dataset_results_main.pdf, combines four cohort panels (A-D) and four aligned GSE275824 patient-trajectory panels (E-H).
- Longitudinal panels use the original deposited ABMRact, ABMRIFNG, ABMRNK and ABMREndo values and report patient-blocked treatment-by-time estimates with BH-adjusted FDR.
- Figure 2, estimand_results_main.pdf, contains the inferential synthesis, pathway atlas and platform-specific PCA.
- These are the only two graphical deliverables; plotted values and coordinates remain as audit CSVs.
- PCA is platform-specific, uses 1,000 most-variable genes, and is descriptive rather than predictive.
- KNN is not used as clustering; supervised prediction requires a separate target and external validation.

## Results directory

- `results/bulk/`: cohort-specific DE, GSEA and treatment-model results.
- `results/atlas/`: cross-cohort mechanism effects and run manifest.
- `results/figures/`: PNG/PDF figures and exact plotted values or coordinates.
- `results/tables/`: cohort composition, estimands, sample labels and module coverage.
- `results/logs/three_dataset_pipeline.log`: complete execution record.
- `results/output_checksums.sha256`: integrity hashes for principal deliverables.

## Reproduce verified analyses

```bash
python3 scripts/build_public_data.py --dataset GSE50084
python3 scripts/build_bulk_cohorts.py --dataset GSE36059
python3 scripts/build_bulk_cohorts.py --dataset GSE275824

Rscript scripts/02b_annotate_bioc.R data/processed/GSE50084/expression_probe_matrix.csv.gz data/processed/GSE50084/samples.csv GPL6244 data/processed/GSE50084
Rscript scripts/02b_annotate_bioc.R data/processed/GSE36059/expression_probe_matrix.csv.gz data/processed/GSE36059/samples.csv GPL570 data/processed/GSE36059
Rscript scripts/02b_annotate_bioc.R data/processed/GSE275824/expression_probe_matrix.csv.gz data/processed/GSE275824/samples.csv GPL13667 data/processed/GSE275824

Rscript scripts/03_bulk_de_gsea.R data/processed/GSE50084/gene_expression.csv.gz data/processed/GSE50084/samples.csv results/bulk/GSE50084
Rscript scripts/04_external_validation.R data/processed/GSE36059/gene_expression.csv.gz data/processed/GSE36059/samples.csv results/bulk/GSE36059
Rscript scripts/05_treatment_response.R data/processed/GSE275824/gene_expression.csv.gz data/processed/GSE275824/samples.csv results/bulk/GSE275824
Rscript scripts/06_published_score_trajectories.R data/processed/GSE275824/samples.csv results/bulk/GSE275824

python3 scripts/07_score_mechanism_modules.py --expression data/processed/GSE36059/gene_expression.csv.gz --output results/bulk/GSE36059/mechanism_module_scores.csv

python3 -m unittest discover -s tests -v
```

Generic single-cell entry points (after downloading an eligible study and preparing barcode-level
metadata) are `scripts/10_build_scanpy_object.py`, `scripts/11_sc_qc_cluster.py`, and
`scripts/12_patient_pseudobulk.py`. `scripts/13_grouped_nested_ml.py` accepts only frozen,
patient-level features and keeps patients together during nested validation.

Required R packages are listed in [`environment.yml`](environment.yml). Raw and derived assay files
are gitignored and rebuilt from their public sources. The code is research software, not a medical
device or clinical decision system.
