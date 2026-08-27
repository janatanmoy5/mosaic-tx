#!/usr/bin/env python3
"""Score prespecified MOSAIC-Tx mechanism modules in a gene-by-sample expression matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expression", type=Path, required=True,
                    help="CSV/CSV.GZ with gene symbols in first column and samples in remaining columns")
    ap.add_argument("--modules", type=Path, default=Path("configs/mechanism_modules.json"))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    expr = pd.read_csv(args.expression, index_col=0)
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.groupby(level=0).mean()
    cfg = json.loads(args.modules.read_text())
    threshold = float(cfg["rules"]["minimum_gene_coverage_fraction"])
    rows, coverage = {}, []
    for name, spec in cfg["modules"].items():
        declared = [g.upper() for g in spec["genes"]]
        present = [g for g in declared if g in expr.index]
        fraction = len(present) / len(declared)
        coverage.append({"module": name, "n_declared": len(declared),
                         "n_present": len(present), "coverage_fraction": fraction,
                         "present_genes": ";".join(present)})
        if fraction < threshold:
            continue
        block = expr.loc[present].astype(float)
        sd = block.std(axis=1, ddof=1).replace(0, np.nan)
        z = block.sub(block.mean(axis=1), axis=0).div(sd, axis=0)
        rows[name] = z.mean(axis=0, skipna=True)
    if not rows:
        raise ValueError("no mechanism module met the configured gene-coverage threshold")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).rename_axis("sample_id").to_csv(args.output)
    pd.DataFrame(coverage).to_csv(args.output.with_suffix(".coverage.csv"), index=False)


if __name__ == "__main__":
    main()
