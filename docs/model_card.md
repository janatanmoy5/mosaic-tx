# Model card: MOSAIC-Tx reference scorer 0.1

## Intended use

Method development and retrospective research prioritization of candidate interventions for
laboratory testing. Users are computational biologists working with transplant clinicians,
pharmacists, pathologists, statisticians, and data-governance teams.

## Prohibited use

Diagnosis, treatment selection, prescribing, dosing, replacing Banff/clinical adjudication,
estimating individual benefit, or use outside an approved research protocol.

## Current model

The reference implementation computes confidence/evidence-weighted negative cosine similarity
between disease and perturbation effects within matched cell-state features. It multiplies the
reversal score by measured disease-state coverage and subtracts separate, configurable safety
penalties. Low-coverage pairs are ineligible. It does not impute missing perturbation effects.

## Principal limitations

- Anticorrelation is not causality and may reflect toxic suppression.
- In vitro perturbations may not transport to a graft under immunosuppression.
- Cell labels, reference choices, tissue dissociation, treatment timing, and batch effects can
  dominate the score.
- A normalized safety penalty is not a clinical risk estimate and cannot encode patient-specific
  contraindications.
- The simple score has no dose, pharmacokinetic, temporal, combination, or spatial mechanism.
- Underrepresented ancestry, age, sex, organ, center, and rare endotypes can produce hidden failure.

## Required reporting

Always show reversal, coverage, safety components, evidence provenance, uncertainty, subgroup
performance, missingness, and abstentions. Never report a candidate rank without the model version
and input hashes.

