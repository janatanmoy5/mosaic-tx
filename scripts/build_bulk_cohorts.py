#!/usr/bin/env python3
"""Acquire and compile the two non-discovery bulk GEO cohorts."""

import argparse
import json
import sys
from pathlib import Path

# Import the standalone compiler from the repository root without installation.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geo_utils import compile_geo, download_verified


def main() -> None:
    # Only cohorts with explicit phenotype parsers are accepted.
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, choices=["GSE36059", "GSE275824"])
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    # The registry pins URL, checksum, platform and expected sample count.
    specification = json.loads((ROOT / "configs/datasets.json").read_text())["datasets"][args.dataset]
    raw_path = ROOT / "data/raw" / args.dataset / f"{args.dataset}_series_matrix.txt.gz"
    output_path = ROOT / "data/processed" / args.dataset

    # A corrupt cache or changed source stops execution before analysis.
    download_verified(specification["url"], raw_path, specification["sha256"], args.force_download)

    # Cohort-specific mapping retains all original diagnostic or trial strata.
    manifest = compile_geo(raw_path, output_path, args.dataset, specification["expected_samples"])
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
