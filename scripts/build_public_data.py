#!/usr/bin/env python3
"""Download, checksum, compile and validate the GSE50084 discovery cohort."""

import argparse
import json
import sys
from pathlib import Path

# Use the standalone GEO compiler without installing a Python package.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geo_utils import compile_gse50084, download_verified, validate_compiled


def main() -> int:
    # This entry point accepts only the cohort whose phenotype labels it understands.
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="GSE50084", choices=["GSE50084"])
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    # Resolve checksum-pinned source and deterministic output locations.
    specification = json.loads((ROOT / "configs/datasets.json").read_text())["datasets"][args.dataset]
    raw_path = ROOT / "data/raw" / args.dataset / f"{args.dataset}_series_matrix.txt.gz"
    output_path = ROOT / "data/processed" / args.dataset

    print(f"Downloading/verifying {args.dataset} from NCBI GEO...")
    download_verified(specification["url"], raw_path, specification["sha256"], args.force_download)

    print("Compiling metadata and submitted normalized expression matrix...")
    compile_gse50084(raw_path, output_path, specification["expected_samples"])
    manifest = validate_compiled(output_path, specification["expected_samples"])
    print(json.dumps(manifest, indent=2))
    print(f"Validated outputs: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
