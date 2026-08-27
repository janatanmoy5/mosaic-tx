# Three-dataset pipeline execution record

## One-command entry point

```bash
cd /Users/janat/Desktop/project_plan/files/mosaic-tx
scripts/run_three_dataset_pipeline.sh
```

The runner uses `set -euo pipefail`, checks command and R-package dependencies before starting,
captures terminal output in `results/logs/three_dataset_pipeline.log`, and returns success only after
strict output validation.

## Stages executed

1. Checksum-verify and compile GSE50084, GSE36059, and GSE275824.
2. Map platform probes using `hugene10sttranscriptcluster.db`, `hgu133plus2.db`, and `hgu219.db`.
3. Run prespecified moderated differential-expression models and Hallmark GSEA.
4. Run patient-blocked randomized difference-in-differences at weeks 24 and 52.
5. Fit published molecular-score mixed models with patient random intercepts.
6. Construct nine frozen module scores independently within every platform.
7. Build the direction-aware robustness and ABMR-specificity atlas.
8. Render the two-panel atlas and four-panel main manuscript figure in PNG and PDF.
9. Run Python unit tests, Python/R syntax checks, GEO compilation validation, and artifact validation.
10. Write SHA-256 checksums for principal tables and figures.

## Successful local execution

The complete command was executed successfully on the Desktop repository. Verified inputs were:

| Dataset | Profiles | Platform | Role |
|---|---:|---|---|
| GSE50084 | 115 total; 61 biopsy profiles in primary contrasts | GPL6244 | DSA-state discovery |
| GSE36059 | 411 | GPL570 | ABMR replication and TCMR specificity |
| GSE275824 | 63 biopsies from 22 randomized patients | GPL13667 | longitudinal treatment response |

Annotation produced 19,998, 21,358, and 19,434 unique genes, respectively. The validator requires
large DE tables, 12 published-score interaction rows, 63 module-effect rows, nine robustness rows,
the exact three-dataset manifest, and both raster and vector manuscript figures.

## Main outputs

- `results/atlas/module_effects.csv`
- `results/atlas/robustness_atlas.csv`
- `results/atlas/module_coverage.csv`
- `results/atlas/run_manifest.json`
- `results/figures/estimand_results_main.pdf`
- `results/figures/dataset_results_main.pdf`
- `results/tables/Table_S1_sample_labels.csv`
- `results/output_checksums.sha256`

## Nonfatal warnings observed

- Installed `data.table`, `fgsea`, and `msigdbr` packages were built under a slightly newer R 4.4
  patch release than the local R 4.4.2 runtime.
- Bioconductor reported expected one-to-many probe-to-symbol mappings; the declared annotation script
  retains its prespecified first mapping and averages probes per symbol.
- `fgseaMultilevel` reported extremely small or potentially overestimated pathway P values and a very
  small number of tied statistics. Exact gene/module conclusions do not depend on those pathway tails.

These warnings are recorded rather than hidden. They did not cause a failed validation stage.

## Scope boundary

This command is the complete **three-public-bulk-dataset analysis**. Single-cell compilation, patient
pseudobulk, ML, drug reversal, QSAR, docking, and MD require additional inputs and separate evidence
gates. They are not silently simulated or described as executed by this runner.
