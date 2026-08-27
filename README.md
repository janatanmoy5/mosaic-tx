# MOSAIC-Tx

MOSAIC-Tx is a reproducible public-data analysis of antibody-mediated kidney-allograft rejection. It separates DSA exposure, damaging rejection, independent replication, TCMR specificity, and randomized molecular response.

## Start here

1. [Abstract](docs/ABSTRACT.md)
2. [Project summary](docs/PROJECT_SUMMARY.md)
3. [Methods](docs/METHODS.md)
4. [Results](docs/RESULTS.md)

## Main outputs

- [Figure 1](results/figures/dataset_results_main.pdf)
- [Figure 2](results/figures/estimand_results_main.pdf)
- [Cohort composition](results/tables/Table_1_cohort_composition.csv)
- [Primary estimands](results/tables/Table_2_primary_estimands.csv)
- [Auditable sample labels](results/tables/Table_S1_sample_labels.csv)
- [Complete derived results](results/)

## Run the verified pipeline

```bash
bash scripts/run_three_dataset_pipeline.sh
```

The pipeline compiles three public cohorts, maps probes to genes, performs moderated differential expression and Hallmark GSEA, scores prespecified modules, estimates randomized longitudinal effects, generates figures and tables, runs tests, validates outputs, and writes checksums.

## Principal conclusion

Immune and endothelial programs distinguish DSA-positive rejection from DSA-positive non-rejection and replicate against non-rejection. Major immune programs are shared with or stronger in TCMR and are not ABMR-specific biomarkers. Randomized data support transient week-24 suppression of immune and deposited MMDx programs, with unresolved endothelial change and no supported week-52 effect.

## Scope

This repository contains completed public bulk-transcriptomic analyses. It does not contain completed single-cell, spatial, clinical machine-learning, drug-repositioning, QSAR, docking, molecular-dynamics, wet-laboratory, or clinical-validation results. It does not support treatment decisions.

Citation metadata are provided in `CITATION.cff`. See `LICENSE` for reuse terms.
