# Feasibility and scalability audit

## What has been executed on this Mac

| Stage | Input size | Peak-scale task | Status |
|---|---:|---|---|
| GSE50084 | 13 MB compressed | 115 samples × 29,096 probes | Passed |
| GSE36059 | 115 MB compressed | 411 samples × 54,675 probes | Passed |
| GSE275824 | 4.9 MB compressed | 63 samples × 49,293 probes | Passed |
| Gene annotation | three array platforms | 19,434–21,358 genes | Passed |
| Limma/GSEA | three cohorts | cohort-level contrasts | Passed |

The bulk workflow uses roughly hundreds of megabytes on disk and completes on a modern Mac. It is
scalable for these datasets but is not a general clinical platform.

## Single-cell feasibility

- E-MTAB-12051: 16 kidney-allograft biopsies; official count CSV is 3,092,166,324 bytes plus metadata.
  A dense CSV import is unsuitable for a low-memory laptop. Planned solution: chunked conversion to a
  sparse on-disk format, then patient-level pseudobulk. Recommended: 32–64 GB RAM and 15–25 GB free disk.
- GSE145927: 81,139 cells from five transplant biopsies; processed archive approximately 247 MB plus
  2.1 MB metadata. Technically feasible on a 16–32 GB Mac, but biologically underpowered for flexible
  patient-level models.

No deep-learning model is justified by these donor counts. Cell-level train/test splits are prohibited.

## Spatial feasibility

Verified public transplant spatial data are too small or disease-mismatched for a central ABMR claim.
GSE252137 has 12 TCMR/non-TCMR ROIs. It can support localization or specificity checks only. Spatial
analysis is therefore optional and exploratory.

## ML feasibility

The bulk cohorts can support penalized pathway-level models with nested cohort-aware evaluation, but
not a high-capacity deep network. Because labels, platforms and cohort definitions differ, external
validation is more important than optimizing an internal cross-validation score. ML is deferred until
feature definitions and cross-cohort mappings are frozen.

## Docking/MD feasibility

No target currently passes a biological target-selection gate. Docking and MD are therefore not part
of the executable primary pipeline. Adding them now would be technically possible but scientifically
unjustified.
