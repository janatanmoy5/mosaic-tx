#!/usr/bin/env python3
"""Aggregate raw single-cell counts to patient/sample/cell-type pseudobulk profiles."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--cell-type-column", default="cell_type")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min-cells", type=int, default=20)
    args = ap.parse_args()

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    obj = sc.read_h5ad(args.input)
    keys = ["dataset_id", "patient_id", "sample_id", "phenotype", args.cell_type_column]
    absent = [x for x in keys if x not in obj.obs]
    if absent:
        raise ValueError(f"missing grouping columns: {absent}")
    counts = obj.layers.get("counts")
    if counts is None:
        raise ValueError("counts layer is required for pseudobulk")

    matrices, rows = [], []
    frame = obj.obs.reset_index(drop=False).rename(columns={obj.obs.index.name or "index": "barcode"})
    for values, idx in frame.groupby(keys, observed=True, sort=True).indices.items():
        if len(idx) < args.min_cells:
            continue
        summed = counts[idx].sum(axis=0)
        matrices.append(sparse.csr_matrix(summed))
        rows.append(dict(zip(keys, values)) | {"n_cells": len(idx)})
    if not matrices:
        raise ValueError("no pseudobulk group passed min-cells")
    X = sparse.vstack(matrices).tocsr()
    pb = ad.AnnData(X=X, obs=pd.DataFrame(rows), var=obj.var.copy())
    pb.layers["counts"] = X.copy()
    pb.uns["replicate_unit"] = "patient/sample within cell type"
    pb.uns["warning"] = "Do not treat multiple cell types from one patient as independent folds."
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pb.write_h5ad(args.output, compression="gzip")


if __name__ == "__main__":
    main()
