#!/usr/bin/env python3
"""Download, verify, and compile the declared multicohort GEO expression datasets."""

import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geo_utils import compile_geo, download_verified
p=argparse.ArgumentParser();p.add_argument('--dataset',required=True,choices=['GSE36059','GSE275824']);p.add_argument('--force-download',action='store_true');a=p.parse_args()
spec=json.loads((ROOT/'configs/datasets.json').read_text())['datasets'][a.dataset]
raw=ROOT/'data/raw'/a.dataset/f'{a.dataset}_series_matrix.txt.gz';out=ROOT/'data/processed'/a.dataset
download_verified(spec['url'],raw,spec['sha256'],a.force_download)
print(json.dumps(compile_geo(raw,out,a.dataset,spec['expected_samples']),indent=2))
