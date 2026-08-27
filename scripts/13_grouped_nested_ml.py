#!/usr/bin/env python3
"""Leakage-resistant patient-grouped nested CV for frozen biomarker features."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", type=Path, required=True, help="CSV: rows are biopsies/patients")
    ap.add_argument("--target", required=True)
    ap.add_argument("--group", default="patient_id")
    ap.add_argument("--feature-list", type=Path, required=True,
                    help="one prespecified feature name per line")
    ap.add_argument("--positive-label", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--outer-folds", type=int, default=5)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from sklearn.calibration import calibration_curve
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
    from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    df = pd.read_csv(args.features)
    names = [x.strip() for x in args.feature_list.read_text().splitlines() if x.strip()]
    absent = [x for x in names + [args.target, args.group] if x not in df]
    if absent:
        raise ValueError(f"missing columns: {absent}")
    y = (df[args.target].astype(str) == args.positive_label).astype(int).to_numpy()
    groups = df[args.group].astype(str).to_numpy()
    if len(np.unique(groups)) != len(groups):
        pass  # expected for serial biopsies; grouped folds prevent leakage
    if min(np.bincount(y)) < args.outer_folds:
        raise ValueError("too few samples in a class for requested outer folds")

    X = df[names]
    outer = StratifiedGroupKFold(n_splits=args.outer_folds, shuffle=True, random_state=2026)
    pred = np.full(len(df), np.nan)
    folds, coefficients = [], []
    for fold, (train, test) in enumerate(outer.split(X, y, groups), start=1):
        inner = StratifiedGroupKFold(n_splits=max(2, args.outer_folds - 1), shuffle=True,
                                     random_state=3000 + fold)
        pipe = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", LogisticRegression(penalty="elasticnet", solver="saga",
                                          max_iter=10000, random_state=fold)),
        ])
        grid = GridSearchCV(pipe, {"model__C": [0.01, 0.1, 1, 10],
                                   "model__l1_ratio": [0.0, 0.5, 1.0]},
                            scoring="average_precision", cv=inner, n_jobs=-1)
        grid.fit(X.iloc[train], y[train], groups=groups[train])
        pred[test] = grid.predict_proba(X.iloc[test])[:, 1]
        coefs = grid.best_estimator_.named_steps["model"].coef_[0]
        coefficients.append(dict(zip(names, map(float, coefs))))
        folds.append({"fold": fold, "n_train": len(train), "n_test": len(test),
                      "best_parameters": grid.best_params_})

    result = {
        "analysis_unit": "row with patient-grouped splitting",
        "model": "elastic-net logistic regression",
        "features_frozen_from": str(args.feature_list),
        "n_rows": len(df), "n_patients": int(len(np.unique(groups))),
        "auroc": float(roc_auc_score(y, pred)),
        "auprc": float(average_precision_score(y, pred)),
        "brier": float(brier_score_loss(y, pred)),
        "folds": folds,
        "coefficient_nonzero_frequency": {
            name: float(np.mean([abs(c[name]) > 1e-12 for c in coefficients])) for name in names
        },
        "limitations": [
            "Internal nested cross-validation is not external clinical validation.",
            "Calibration and confidence intervals require adequate independent patient numbers.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    pd.DataFrame({"observed": y, "predicted": pred, args.group: groups}).to_csv(
        args.output.with_suffix(".predictions.csv"), index=False)


if __name__ == "__main__":
    main()
