#!/usr/bin/env python3
"""Download, checksum, compile, and validate the GSE50084 discovery cohort."""

import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from mosaic_tx.geo import compile_gse50084, download_verified, validate_compiled

def main():
 p=argparse.ArgumentParser(); p.add_argument("--dataset",default="GSE50084",choices=["GSE50084"]); p.add_argument("--force-download",action="store_true"); a=p.parse_args()
 spec=json.loads((ROOT/"configs/datasets.json").read_text())["datasets"][a.dataset]
 raw=ROOT/"data/raw"/a.dataset/f"{a.dataset}_series_matrix.txt.gz"; out=ROOT/"data/processed"/a.dataset
 print(f"Downloading/verifying {a.dataset} from NCBI GEO...")
 download_verified(spec["url"],raw,spec["sha256"],a.force_download)
 print("Compiling metadata and submitted normalized expression matrix...")
 compile_gse50084(raw,out,spec["expected_samples"]); manifest=validate_compiled(out,spec["expected_samples"])
 print(json.dumps(manifest,indent=2)); print(f"Validated outputs: {out}"); return 0
if __name__=="__main__": raise SystemExit(main())
