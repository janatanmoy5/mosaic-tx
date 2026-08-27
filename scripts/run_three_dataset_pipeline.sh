#!/usr/bin/env bash

# Execute the complete, public three-dataset MOSAIC-Tx bulk pilot with fail-fast behavior.
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
mkdir -p results/logs results/pilot/figures
log_path="results/logs/three_dataset_pipeline.log"
exec > >(tee "$log_path") 2>&1
echo "MOSAIC-Tx three-dataset pipeline started"

# Run acquisition/checksum validation, annotation, inference, score models, and atlas generation.
make data
make annotate
make bulk
make scores
make pilot

# Render the four-panel main figure after statistical results have been frozen to disk.
Rscript scripts/09b_make_manuscript_figure.R . results/pilot/atlas results/pilot/figures

# Run unit/syntax checks and strict artifact validation.
make check
python3 scripts/14_validate_pipeline_outputs.py --root .

# Record checksums for the principal machine-readable and graphical outputs.
shasum -a 256 \
  results/pilot/atlas/module_effects.csv \
  results/pilot/atlas/robustness_atlas.csv \
  results/pilot/atlas/run_manifest.json \
  results/pilot/figures/main_figure_multicohort.png \
  results/pilot/figures/main_figure_multicohort.pdf \
  > results/pilot/output_checksums.sha256

echo "MOSAIC-Tx three-dataset pipeline completed successfully"
