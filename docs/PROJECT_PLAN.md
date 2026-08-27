# MOSAIC-Tx project plan

## Working title

**Reference-robust and cell-state-aware qualification of molecular signatures for computational drug repositioning in antibody-mediated kidney-allograft rejection**

The binding question-to-estimand and figure-to-claim specification is in `STUDY_BLUEPRINT.md`.

## Central question

Which ABMR expression programs are reproducible, rejection-specific, cell-localizable, and sufficiently robust across reference definitions and treatment states to justify downstream drug-reversal analysis?

## Executed pilot status

The current bulk pilot found reproducible rejection-associated NK, FCGR3A/Fc, cytotoxic, IFN/chemokine and endothelial modules across two cohorts, but no module was positively enriched in ABMR versus TCMR. Short-term felzartamab-associated reductions were detected for NK, cytotoxic and IFN/chemokine modules at week 24; no module retained an FDR-significant difference at week 52. These findings qualify rather than validate inputs for drug reversal. See `PILOT_RESULTS.md`.

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

## Aim 4 — Localize replicated programs to cells

Use human NovaSeq 6000 E-MTAB-12051 as the primary biopsy cell-localization resource and GSE216203
as a circulating-NK mechanistic sensitivity dataset. GSE145927 was generated on HiSeq and is retained
as a technology sensitivity analysis, not merged into the NovaSeq discovery object.
Inference is performed by patient-level pseudobulk within prespecified cell classes. Cells are not
independent replicates. Because E-MTAB-12051 is a 3.09 GB count CSV and GSE145927 includes only five
biopsy donors, this aim requires a separate 32–64 GB RAM environment and is not yet completed.

The canonical object is Scanpy `AnnData`, with raw counts retained and optional MuData export. Stable
allograft, non-rejection injury and native healthy kidney remain distinct reference domains. Cell
labels require markers, negative markers, reference mapping and manual cluster review; KIR expression
alone cannot define NK cells.

## Aim 5 — Build and externally challenge a parsimonious biomarker

Freeze features only after bulk replication and patient-level cell localization. Train penalized
logistic regression as the primary model, with all observations from a patient kept in the same fold.
Feature selection and tuning occur inside nested group-aware cross-validation. The decisive test is an
untouched study or leave-one-study-out validation, including calibration and specificity versus TCMR.
If independent-patient numbers are inadequate, report association and do not train a classifier.

## Aim 6 — Generate intervention hypotheses under explicit gates

Query replicated up/down signatures against perturbational resources and retain compounds only when
reversal is robust across contexts and target pharmacology agrees with the supported NK–Fc–endothelial
mechanism. QSAR proceeds only for a single tractable target with an adequate, comparable public human
activity series and scaffold-held-out evaluation. Docking/MD remain optional structural follow-up,
not biological validation. See `TRANSLATIONAL_EXTENSION.md` and `DRUG_QSAR_PROTOCOL.md`.

## Spatial evidence

Spatial data are not a primary validation aim. GSE252137 contains 12 TCMR/non-TCMR ROIs and can test
anatomical plausibility or specificity, not establish ABMR replication. The previously cited
GSE280559 is unrelated to transplantation and has been removed. ABMR spatial case studies with one
patient are literature support only.

## Primary endpoints

- Aim 1: moderated log-fold changes and Hallmark enrichment for the three declared contrasts.
- Aim 2: independent direction concordance and ABMR-versus-TCMR specificity.
- Aim 3: placebo-adjusted longitudinal interaction for published molecular scores and expression.
- Aim 4: patient-level pseudobulk direction concordance in prespecified cell classes.
- Aim 5: grouped nested validation, calibration and untouched-study performance.
- Aim 6: robust signature reversal; target-specific QSAR only after the chemical-data gate.

## Success and stop criteria

Proceed to mechanistic synthesis only if a program is supported in at least two independent cohorts
and is not explained solely by TCMR. Stop target prioritization if cell localization is inconsistent,
treatment direction is unstable, or results depend on one reference definition. Do not initiate
docking or MD without a replicated target, known directional biology, a tractable structure and a
reference-ligand benchmark.

## Expected paper figures

1. Dataset architecture, contrasts and evidence gates.
2. DSA-state continuum in GSE50084, including the null difficult-control result.
3. External ABMR replication and ABMR-versus-TCMR specificity.
4. Randomized longitudinal immune, NK and endothelial score trajectories.
5. Patient-level single-cell localization, only after successful computation.
6. Cross-cohort evidence synthesis and falsification tests.
7. Patient-level biomarker validation, only if sample size permits.
8. Drug-reversal evidence ladder and QSAR feasibility decision—not an efficacy claim.

## Journal positioning

A public-data-only study should be judged after completion. Reasonable aspirational venues include
American Journal of Transplantation, Kidney International Reports, Communications Biology, or a
strong computational biology journal. Nature Communications or Cell Reports Medicine would require
substantial multimodal replication and a clearly new biological conclusion. NEJM and primary
Nature/Cell are not realistic without major prospective or experimental evidence.
