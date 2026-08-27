# Project summary

## Working title

**Reference-robust and cell-state-aware qualification of molecular signatures for computational drug repositioning in antibody-mediated kidney-allograft rejection**

The binding question-to-estimand and figure-to-claim specification is in `STUDY_BLUEPRINT.md`.

## Central question

Which ABMR expression programs are reproducible, rejection-specific, cell-localizable, and sufficiently robust across reference definitions and treatment states to justify downstream drug-reversal analysis?

## Executed analysis status

The current bulk analysis found reproducible rejection-associated NK, FCGR3A/Fc, cytotoxic, IFN/chemokine and endothelial modules across two cohorts, but no module was positively enriched in ABMR versus TCMR. Short-term felzartamab-associated reductions were detected for NK, cytotoxic and IFN/chemokine modules at week 24; no module retained an FDR-significant difference at week 52. These findings qualify rather than validate inputs for drug reversal. See `RESULTS.md`.

## Biological model to test—not assume

The study tests the following explicit chain:

```text
DSA → endothelial antibody engagement → FCGR3A/CD16-positive NK and monocyte activation
→ cytotoxicity plus IFN-gamma/chemokines → microvascular endothelial injury
→ persistent versus treatment-responsive graft injury
```

Every arrow has an observable proxy, limitation and falsification test in `MECHANISM_MAP.md`. The final
trajectory is summarized by a two-component model:

1. **Immune-effector activity:** IFN-gamma-responsive, NK-cell/Fc-receptor and antibody-associated
   programs that may decrease during CD38-directed therapy.
2. **Tissue-injury persistence:** endothelial, microvascular and parenchymal injury programs that may
   remain when immune-effector activity falls.

The model is rejected if the components are not reproducible across independent cohorts, cannot be
localized at patient level in single-cell data, or do not show distinguishable longitudinal behavior.
No transcriptomic association is called causal.

## Aim 1 — Separate DSA exposure from damaging rejection

Primary GSE50084 contrasts:

- DSA-positive/rejection versus DSA-positive/no-rejection
- DSA-positive/no-rejection versus DSA-negative/no-rejection
- DSA-positive/rejection versus DSA-negative/no-rejection

The difficult DSA-positive/no-rejection group is central, not a disposable control. Success requires
interpretable effect estimates and pathway results; a null contrast is retained and may falsify a
simple binary model.

## Aim 2 — Replicate ABMR biology and test rejection specificity

In GSE36059, compare ABMR with non-rejecting allografts and separately with TCMR. Nephrectomies and
mixed rejection remain distinct. A program is considered replicated only when direction agrees,
FDR is controlled in the validation cohort, and the effect is not merely a pan-rejection signal.

## Aim 3 — Test immune-versus-tissue response during randomized therapy

In GSE275824, analyze baseline, week 24 and week 52 using patient blocking and a placebo-adjusted
change-from-baseline contrast. Analyze both genome-wide expression and the study-provided ABMR,
IFN-gamma, NK and endothelial scores. The key test is whether immune and endothelial/tissue scores
show different trajectories; the analysis does not estimate individual treatment benefit.

## Integrated conclusion

The executed analysis supports rejection-associated immune and endothelial programs, but it does not support ABMR specificity for the principal immune modules. Randomized data support a short-term week-24 reduction in NK, cytotoxic, interferon-gamma, and deposited MMDx ABMR programs, without a statistically resolved endothelial effect or supported week-52 persistence.

The most defensible biological model is therefore a dynamically suppressible immune component accompanied by a less clearly treatment-responsive endothelial component. This is a hypothesis for prospective validation, not proof of cell origin, causality, clinical utility, or drug efficacy.

## Scope and next evidence gate

The current completed work is the three-cohort public bulk-transcriptomic analysis. Single-cell localization, spatial analysis, clinical machine learning, drug reversal, QSAR, docking, molecular dynamics, wet-laboratory validation, and clinical validation remain uncompleted and must not be presented as results.

The next scientific priority is independent review of sample labels and statistical models, followed by a suitable external cohort and patient-level cell-resolved validation.
