#!/usr/bin/env python3
"""Normalize, embed and cluster a MOSAIC-Tx AnnData object without assigning automatic truth labels."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--resolution", type=float, default=0.6)
    ap.add_argument("--n-hvg", type=int, default=3000)
    args = ap.parse_args()

    import scanpy as sc

    obj = sc.read_h5ad(args.input)
    if "counts" not in obj.layers:
        raise ValueError("counts layer is required")
    sc.pp.normalize_total(obj, target_sum=1e4)
    sc.pp.log1p(obj)
    obj.raw = obj
    sc.pp.highly_variable_genes(obj, n_top_genes=min(args.n_hvg, obj.n_vars),
                                flavor="seurat", subset=False)
    sc.tl.pca(obj, use_highly_variable=True)
    sc.pp.neighbors(obj, n_neighbors=15, n_pcs=min(40, obj.obsm["X_pca"].shape[1]))
    sc.tl.umap(obj)
    sc.tl.leiden(obj, resolution=args.resolution, key_added="leiden")
    sc.tl.rank_genes_groups(obj, "leiden", method="wilcoxon", key_added="cluster_markers")
    obj.uns["cell_type_status"] = (
        "unreviewed: assign labels only after marker, negative-marker, reference-map and cluster review"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    obj.write_h5ad(args.output, compression="gzip")


if __name__ == "__main__":
    main()
