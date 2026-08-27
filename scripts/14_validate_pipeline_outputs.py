#!/usr/bin/env python3
"""Fail the three-dataset pipeline if required artifacts are absent, empty, or malformed."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def require_nonempty(path: Path) -> None:
    """Require a regular output file with at least one byte."""
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"required output is missing or empty: {path}")


def require_csv_rows(path: Path, minimum_rows: int) -> None:
    """Require a parseable CSV containing a header and the declared minimum data rows."""
    require_nonempty(path)
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        rows = sum(1 for _ in reader)
    if not header or rows < minimum_rows:
        raise RuntimeError(f"CSV failed row validation: {path} ({rows} rows)")


def main() -> None:
    """Validate core cohort, inference, atlas, and manuscript-figure artifacts."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()

    # Validate processed matrices and cohort manifests for all three declared datasets.
    for dataset in ("GSE50084", "GSE36059", "GSE275824"):
        require_nonempty(root / "data" / "processed" / dataset / "gene_expression.csv.gz")
        require_csv_rows(root / "data" / "processed" / dataset / "samples.csv", 20)
        require_nonempty(root / "data" / "processed" / dataset / "annotation_manifest.json")

    # Validate representative differential-expression and treatment outputs.
    require_csv_rows(root / "results/bulk/GSE50084/de_rejection_vs_dsa_negative_control.csv", 1000)
    require_csv_rows(root / "results/bulk/GSE36059/de_ABMR_vs_non_rejecting.csv", 1000)
    require_csv_rows(root / "results/bulk/GSE275824/de_difference_in_differences_W24.csv", 1000)
    require_csv_rows(root / "results/bulk/GSE275824/published_score_interactions.csv", 12)

    # Validate the complete module atlas dimensions and deterministic manifest.
    require_csv_rows(root / "results/pilot/atlas/module_effects.csv", 63)
    require_csv_rows(root / "results/pilot/atlas/robustness_atlas.csv", 9)
    manifest_path = root / "results/pilot/atlas/run_manifest.json"
    require_nonempty(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if set(manifest["datasets"]) != {"GSE50084", "GSE36059", "GSE275824"}:
        raise RuntimeError("pilot manifest does not contain exactly the three declared cohorts")

    # Validate both raster and vector forms of the main paper figure.
    require_nonempty(root / "results/pilot/figures/main_figure_multicohort.png")
    require_nonempty(root / "results/pilot/figures/main_figure_multicohort.pdf")
    print("Pipeline output validation passed for all three datasets.")


if __name__ == "__main__":
    main()
