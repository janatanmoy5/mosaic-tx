"""Transparent reference implementation of cell-state-aware reversal scoring."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import math

from .io import DataContractError, finite_number


@dataclass(frozen=True)
class Weights:
    infection: float = 0.20
    nephrotoxicity: float = 0.25
    interaction: float = 0.15
    minimum_coverage: float = 0.25

    def validate(self) -> None:
        values = (self.infection, self.nephrotoxicity, self.interaction)
        if any(value < 0 for value in values):
            raise ValueError("safety weights must be non-negative")
        if not 0 <= self.minimum_coverage <= 1:
            raise ValueError("minimum_coverage must be between 0 and 1")


def _bounded(row: dict[str, str], field: str, context: str) -> float:
    value = finite_number(row, field, context)
    if not 0 <= value <= 1:
        raise DataContractError(f"{context}: {field} must be between 0 and 1")
    return value


def rank_interventions(
    disease_rows: list[dict[str, str]],
    perturbation_rows: list[dict[str, str]],
    safety_rows: list[dict[str, str]],
    weights: Weights = Weights(),
) -> list[dict[str, object]]:
    """Rank each intervention per patient using weighted negative cosine similarity.

    Disease effects are signed deviations from a stable reference. Perturbation effects are signed
    changes caused by an intervention. Perfect opposition yields reversal=1; perfect alignment=-1.
    Confidence/evidence values weight each overlapping feature. Missing evidence is not imputed.
    """
    weights.validate()
    disease: dict[str, dict[tuple[str, str], tuple[float, float]]] = defaultdict(dict)
    for index, row in enumerate(disease_rows, start=2):
        context = f"disease row {index}"
        key = (row["cell_state"].strip(), row["feature"].strip())
        patient = row["patient_id"].strip()
        if not patient or not all(key):
            raise DataContractError(f"{context}: identifiers cannot be blank")
        if key in disease[patient]:
            raise DataContractError(f"{context}: duplicate patient/cell_state/feature")
        disease[patient][key] = (
            finite_number(row, "effect", context),
            _bounded(row, "confidence", context),
        )

    perturbations: dict[str, dict[tuple[str, str], tuple[float, float]]] = defaultdict(dict)
    for index, row in enumerate(perturbation_rows, start=2):
        context = f"perturbation row {index}"
        drug = row["intervention"].strip()
        key = (row["cell_state"].strip(), row["feature"].strip())
        if not drug or not all(key):
            raise DataContractError(f"{context}: identifiers cannot be blank")
        if key in perturbations[drug]:
            raise DataContractError(f"{context}: duplicate intervention/cell_state/feature")
        perturbations[drug][key] = (
            finite_number(row, "effect", context),
            _bounded(row, "evidence", context),
        )

    safety: dict[str, tuple[float, float, float]] = {}
    for index, row in enumerate(safety_rows, start=2):
        context = f"safety row {index}"
        drug = row["intervention"].strip()
        if drug in safety:
            raise DataContractError(f"{context}: duplicate intervention")
        safety[drug] = (
            _bounded(row, "infection_risk", context),
            _bounded(row, "nephrotoxicity", context),
            _bounded(row, "interaction_risk", context),
        )

    results: list[dict[str, object]] = []
    for patient, patient_state in sorted(disease.items()):
        total_mass = sum(abs(effect) * confidence for effect, confidence in patient_state.values())
        if total_mass == 0:
            raise DataContractError(f"patient {patient}: disease signature has zero weighted mass")
        for drug, drug_state in sorted(perturbations.items()):
            if drug not in safety:
                raise DataContractError(f"intervention {drug}: missing safety row")
            overlap = patient_state.keys() & drug_state.keys()
            overlap_mass = sum(
                abs(patient_state[key][0]) * patient_state[key][1] for key in overlap
            )
            coverage = overlap_mass / total_mass
            numerator = disease_norm = drug_norm = 0.0
            for key in overlap:
                disease_effect, confidence = patient_state[key]
                drug_effect, evidence = drug_state[key]
                reliability = confidence * evidence
                numerator += reliability * disease_effect * drug_effect
                disease_norm += reliability * disease_effect**2
                drug_norm += reliability * drug_effect**2
            denominator = math.sqrt(disease_norm * drug_norm)
            reversal = -numerator / denominator if denominator else 0.0
            infection, kidney, interaction = safety[drug]
            safety_penalty = (
                weights.infection * infection
                + weights.nephrotoxicity * kidney
                + weights.interaction * interaction
            )
            eligible = coverage >= weights.minimum_coverage
            priority = reversal * coverage - safety_penalty if eligible else float("nan")
            results.append(
                {
                    "patient_id": patient,
                    "intervention": drug,
                    "priority_score": round(priority, 6) if eligible else "",
                    "reversal_score": round(reversal, 6),
                    "coverage": round(coverage, 6),
                    "safety_penalty": round(safety_penalty, 6),
                    "overlap_features": len(overlap),
                    "eligible": str(eligible).lower(),
                }
            )
    results.sort(
        key=lambda row: (
            row["patient_id"],
            -(float(row["priority_score"]) if row["priority_score"] != "" else -math.inf),
            row["intervention"],
        )
    )
    return results

