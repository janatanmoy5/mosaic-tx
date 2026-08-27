#!/usr/bin/env python3
"""Build a provenance-rich AnnData object from public human NovaSeq 6000 scRNA-seq data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REQUIRED = {
    "dataset_id", "sample_id", "patient_id", "phenotype", "tissue",
    "species", "platform", "reference_domain",
}


def validate_metadata(obs) -> None:
    missing = sorted(REQUIRED - set(obs.columns))
    if missing:
        raise ValueError(f"metadata missing required columns: {', '.join(missing)}")
    if obs[list(REQUIRED)].isna().any().any():
        raise ValueError("required metadata contains missing values")
    if set(obs["species"].astype(str)) != {"Homo sapiens"}:
        raise ValueError("primary object must contain Homo sapiens only")
    if set(obs["platform"].astype(str)) != {"Illumina NovaSeq 6000"}:
        raise ValueError("primary object must contain Illumina NovaSeq 6000 only")
    if (obs["patient_id"].astype(str).str.strip() == "").any():
        raise ValueError("patient_id cannot be empty: cells are not independent replicates")


def read_input(path: Path, fmt: str):
    import scanpy as sc

    if fmt == "10x-h5":
        return sc.read_10x_h5(path)
    if fmt == "10x-mtx":
        return sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    if fmt == "h5ad":
        return sc.read_h5ad(path)
    raise ValueError(f"unsupported format: {fmt}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--format", choices=["10x-h5", "10x-mtx", "h5ad"], required=True)
    ap.add_argument("--metadata", type=Path, required=True,
                    help="CSV indexed by cell barcode with required provenance columns")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--min-cells", type=int, default=3)
    ap.add_argument("--max-mito-pct", type=float, default=25.0)
    args = ap.parse_args()

    import anndata as ad
    import pandas as pd
    import scanpy as sc

    obj = read_input(args.input, args.format)
    obj.var_names_make_unique()
    meta = pd.read_csv(args.metadata, index_col=0, dtype=str)
    missing_barcodes = obj.obs_names.difference(meta.index)
    if len(missing_barcodes):
        raise ValueError(f"metadata absent for {len(missing_barcodes)} cell barcodes")
    obj.obs = obj.obs.join(meta.loc[obj.obs_names])
    validate_metadata(obj.obs)

    sc.pp.filter_cells(obj, min_genes=args.min_genes)
    sc.pp.filter_genes(obj, min_cells=args.min_cells)
    obj.var["mt"] = obj.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(obj, qc_vars=["mt"], inplace=True)
    obj = obj[obj.obs["pct_counts_mt"] <= args.max_mito_pct].copy()
    obj.layers["counts"] = obj.X.copy()

    obj.uns["mosaic_tx_build"] = {
        "input": str(args.input), "format": args.format,
        "min_genes": args.min_genes, "min_cells": args.min_cells,
        "max_mito_pct": args.max_mito_pct,
        "note": "Doublet and ambient-RNA handling must be documented per dataset.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    obj.write_h5ad(args.output, compression="gzip")

    manifest = args.output.with_suffix(".manifest.json")
    manifest.write_text(json.dumps({
        "cells": obj.n_obs, "genes": obj.n_vars,
        "patients": int(obj.obs["patient_id"].nunique()),
        "samples": int(obj.obs["sample_id"].nunique()),
        "platforms": sorted(obj.obs["platform"].unique().tolist()),
        "output": str(args.output),
        "anndata_version": ad.__version__, "scanpy_version": sc.__version__,
    }, indent=2) + "\n")

    try:
        import muon as mu
        mdata = mu.MuData({"rna": obj})
        mdata.write_h5mu(args.output.with_suffix(".h5mu"))
    except ImportError:
        pass


if __name__ == "__main__":
    main()
