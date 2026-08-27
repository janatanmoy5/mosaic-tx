# MOSAIC-Tx biological and statistical blueprint

## Biological question

Does DSA-associated antibody-mediated kidney-allograft rejection (ABMR) contain a reproducible
Fc/cytotoxic immune-effector program that is distinguishable from both uncomplicated DSA exposure
and T-cell-mediated rejection (TCMR), and is that program preferentially reduced by randomized
CD38-directed treatment relative to endothelial/tissue injury?

This is a qualification study, not a biomarker-discovery contest. The unit of inference is a biopsy
for cross-sectional cohorts and a patient for longitudinal treatment analyses.

## Why the question matters

DSA is neither necessary nor sufficient evidence of active graft injury. Cytotoxic, interferon and
endothelial transcripts can also occur in TCMR and other inflammatory injuries. A drug-reversal
signature is interpretable only if it separates: (i) antibody exposure from tissue injury, (ii) ABMR
from shared rejection biology, and (iii) treatment-responsive immune activity from persistent injury.

## Primary estimands and decision rules

| Question | Cohort and comparison | Estimand | Statistical model | Required evidence | Failure interpretation |
|---|---|---|---|---|---|
| Is DSA alone sufficient? | GSE50084: DSA+/no rejection vs DSA-/no rejection | Standardized module mean difference | limma moderated two-group model | FDR < 0.05 and coherent direction | DSA exposure alone does not activate the program |
| What distinguishes rejection among DSA+ grafts? | GSE50084: DSA+/rejection vs DSA+/no rejection | Standardized module mean difference | limma moderated two-group model | FDR < 0.05 | Module does not distinguish damaging rejection from tolerated DSA |
| Does the program replicate? | GSE36059: ABMR vs non-rejecting | Standardized module mean difference | limma moderated two-group model | Direction agreement and FDR < 0.05 | Discovery association is not externally reproducible |
| Is it ABMR-specific? | GSE36059: ABMR vs TCMR | Standardized module mean difference | limma moderated two-group model | Positive effect and FDR < 0.05 | Program is shared rejection biology or TCMR-enriched |
| Does treatment change it? | GSE275824: felzartamab vs placebo, baseline to week 24/week 52 | Difference-in-differences | patient-blocked limma; mixed model for published scores | interaction FDR < 0.05, direction consistent with suppression | No supported treatment-associated molecular reduction |
| Is immune response preferential to injury response? | GSE275824: NK/IFN effects compared with endothelial effects | Contrast of standardized treatment interactions | score-specific interaction estimates and confidence intervals | immune effects are more negative with non-overlapping uncertainty supported by a direct contrast | Apparent response is global, uncertain, or not compartment-selective |

BH false-discovery-rate correction is applied within each declared contrast family. Effect sizes and
95% confidence intervals are primary; thresholded significance is secondary. No cross-platform raw
expression values are pooled.

## Prespecified biological programs

- **Fc-effector capacity:** FCGR3A, FCER1G, TYROBP, SYK, LCP2 and FCRL3.
- **Cytotoxic execution:** GNLY, PRF1, GZMB, GZMH, CTSW and CCL5.
- **IFN-gamma/chemokine response:** IFNG, STAT1, IRF1, CXCL9, CXCL10, CXCL11, GBP1 and IDO1.
- **Endothelial activation/injury:** SELE, VCAM1, ICAM1, PLVAP, ACKR1, ESM1, HLA-DRA and GBP1.
- **Lineage/admixture controls:** NK identity, T-cell and monocyte modules.
- **Tissue injury/repair control:** HAVCR1, LCN2, SPP1, SOX9, KRT8, KRT18, VIM and TIMP1.

These are outcome-independent modules. Bulk expression cannot determine whether FCGR3A originates
from NK cells or monocytes; cell-specific claims require donor-level single-cell pseudobulk.

## Figure-to-claim contract

1. **Primary Figure A:** three GSE50084 effect-size contrasts. Answers whether DSA exposure and
   rejection are separable.
2. **Primary Figure B:** discovery-versus-validation effect concordance, with ABMR-versus-TCMR
   specificity encoded. Answers whether replication is specific rather than pan-rejection.
3. **Primary Figure C:** Hallmark normalized enrichment scores across all declared contrasts.
   Answers which broader biological processes support or contradict the module result.
4. **Primary Figure D:** randomized published-score interaction estimates at weeks 24 and 52.
   Answers magnitude, uncertainty and durability of treatment-associated change.
5. **Supplementary QC:** PCA and individual distributions. These diagnose cohort structure; they do
   not establish the biological conclusion.

## Current conclusion permitted by the executed pilot

The data support reproducible immune-effector and endothelial activation in rejection relative to
non-rejection. They do not currently support ABMR specificity: the ABMR-versus-TCMR comparison is
null or TCMR-directed for the principal immune programs. Randomized felzartamab data support a
short-term reduction in published ABMR, IFN and NK scores, without convincing week-52 persistence.
Thus, drug-repositioning input should be described as a **treatment-responsive shared rejection
signature**, pending single-cell localization and a direct immune-versus-injury response contrast.

## Advancement gates

- Do not call an NK-cell mechanism until the program replicates in patient-level NK pseudobulk and
  is distinguishable from FCGR3A-positive monocytes and cytotoxic T cells.
- Do not train a clinical classifier without an untouched patient-level external cohort and
  calibration assessment.
- Do not perform QSAR until one human molecular target, action direction and sufficiently large
  comparable activity series have passed the target gate.
- Docking and molecular dynamics can explain a chemically qualified lead; they cannot rescue weak
  transplant biology or establish efficacy.
