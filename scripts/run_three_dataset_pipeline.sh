#!/usr/bin/env bash

# Execute the complete, public three-dataset MOSAIC-Tx bulk analysis with fail-fast behavior.
set -euo pipefail

# Resolve the repository from this script so the command works from any current directory.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# Check executable dependencies before starting a potentially long acquisition or analysis stage.
for command_name in python3 Rscript make shasum; do
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "Missing required command: $command_name" >&2
    exit 2
  }
done

# Check the R packages required by acquisition-independent inference and figure generation.
Rscript -e 'required <- c("data.table","AnnotationDbi","jsonlite","limma","fgsea","msigdbr","ggplot2","nlme","patchwork"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))'

# Record all terminal output while preserving nonzero exit codes from every pipeline stage.
mkdir -p results/logs results/figures
log_path="results/logs/three_dataset_pipeline.log"
exec > >(tee "$log_path") 2>&1
echo "MOSAIC-Tx three-dataset pipeline started"

# Run acquisition/checksum validation, annotation, inference, score models, and atlas generation.
make data
make annotate
make bulk
make scores
make atlas

# Render inferential summary panels and the biopsy/patient-level primary results figure.
Rscript scripts/15_make_dataset_result_figures.R . results/figures
Rscript scripts/16_make_estimand_figure.R . results/atlas results/figures
python3 scripts/17_build_audit_tables.py --root . --output results/tables

# Run unit/syntax checks and strict artifact validation.
make check
python3 scripts/14_validate_pipeline_outputs.py --root .

# Record checksums for the principal machine-readable and graphical outputs.
shasum -a 256 \
  results/atlas/module_effects.csv \
  results/atlas/robustness_atlas.csv \
  results/atlas/run_manifest.json \
  results/figures/dataset_results_main.png \
  results/figures/dataset_results_main.pdf \
  results/figures/estimand_results_main.png \
  results/figures/estimand_results_main.pdf \
  results/figures/cross_cohort_replication_specificity.csv \
  results/figures/selected_hallmark_gsea_atlas.csv \
  results/figures/individual_sample_module_scores.csv \
  results/figures/GSE275824_patient_trajectories.csv \
  results/figures/GSE275824_change_from_baseline.csv \
  results/figures/GSE36059_validation_effects.csv \
  results/figures/GSE50084_descriptive_pca.csv \
  results/figures/GSE36059_descriptive_pca.csv \
  results/tables/Table_1_cohort_composition.csv \
  results/tables/Table_2_primary_estimands.csv \
  results/tables/Table_S1_sample_labels.csv \
  results/tables/Table_S2_module_definitions_and_coverage.csv \
  results/tables/audit_manifest.json \
  > results/output_checksums.sha256

echo "MOSAIC-Tx three-dataset pipeline completed successfully"
