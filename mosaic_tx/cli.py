"""Command-line interface."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from . import __version__
from .io import read_rows, sha256, write_csv
from .scoring import Weights, rank_interventions


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mosaic-tx")
    subparsers = parser.add_subparsers(dest="command", required=True)
    rank = subparsers.add_parser("rank", help="rank interventions per patient")
    rank.add_argument("--disease", required=True)
    rank.add_argument("--perturbations", required=True)
    rank.add_argument("--safety", required=True)
    rank.add_argument("--output", required=True)
    rank.add_argument("--audit", required=True)
    rank.add_argument("--infection-weight", type=float, default=0.20)
    rank.add_argument("--nephrotoxicity-weight", type=float, default=0.25)
    rank.add_argument("--interaction-weight", type=float, default=0.15)
    rank.add_argument("--minimum-coverage", type=float, default=0.25)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    disease_columns = ("patient_id", "cell_state", "feature", "effect", "confidence")
    perturbation_columns = ("intervention", "cell_state", "feature", "effect", "evidence")
    safety_columns = ("intervention", "infection_risk", "nephrotoxicity", "interaction_risk")
    paths = [args.disease, args.perturbations, args.safety]
    weights = Weights(
        infection=args.infection_weight,
        nephrotoxicity=args.nephrotoxicity_weight,
        interaction=args.interaction_weight,
        minimum_coverage=args.minimum_coverage,
    )
    rankings = rank_interventions(
        read_rows(args.disease, disease_columns),
        read_rows(args.perturbations, perturbation_columns),
        read_rows(args.safety, safety_columns),
        weights,
    )
    write_csv(args.output, rankings)
    audit = {
        "schema_version": "1.0",
        "software_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {str(Path(path)): sha256(path) for path in paths},
        "weights": weights.__dict__,
        "rows_written": len(rankings),
        "research_use_only": True,
    }
    audit_path = Path(args.audit)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

