# Executed analysis: reference-aware qualification of ABMR-associated transcript programs

## Status

The bulk analysis was executed on the local Mac using the processed public matrices for GSE50084,
GSE36059 and GSE275824. The run comprised 589 expression profiles in total, including 22 patients
with serial randomized-treatment biopsies. The complete effect table, module coverage, evidence grades, run manifest, two principal figures, and source plotting tables are stored under `results/`. Methods and project interpretation are maintained in `docs/METHODS.md` and `docs/PROJECT_SUMMARY.md`.

This analysis is a stress test of molecular inputs for drug repositioning. It is not an independent
drug-discovery claim and does not duplicate the experimental aims of R01AI180118.

## Principal result

No prespecified biological module achieved Tier A, which required:

1. FDR-significant, same-direction association in GSE50084 and GSE36059; and
2. significant positive enrichment in ABMR compared with TCMR.

The NK, FCGR3A/Fc, cytotoxic, IFN/chemokine and endothelial modules replicated strongly in ABMR
versus non-rejection, but they were not positively enriched in ABMR versus TCMR. Several immune
modules were significantly higher in TCMR. Therefore these bulk signatures are reproducible
rejection-associated programs, not ABMR-specific mechanisms.

## Selected module effects

Effects are within-cohort standardized module-score differences. Positive values favor the first
group in each contrast. FDR is Benjamini-Hochberg adjusted across the nine declared modules within a
contrast.

| Module | DSA+ ABMR vs DSA− non-rejection | Independent ABMR vs non-rejection | ABMR vs TCMR | Felzartamab ΔΔ W24 | Felzartamab ΔΔ W52 | Tier |
|---|---:|---:|---:|---:|---:|---|
| NK identity | 0.459; FDR 3.17e−7 | 0.539; FDR 3.92e−8 | −0.343; FDR 0.060 | −0.914; FDR 0.025 | −0.057; FDR 0.973 | B |
| FCGR3A/Fc signaling | 0.526; FDR 5.19e−8 | 0.423; FDR 2.98e−5 | −0.896; FDR 2.33e−6 | −0.295; FDR 0.494 | −0.227; FDR 0.745 | B |
| Cytotoxicity | 0.416; FDR 9.37e−6 | 0.508; FDR 2.37e−6 | −0.642; FDR 0.00134 | −0.915; FDR 0.025 | −0.003; FDR 0.992 | B |
| IFNγ/chemokines | 1.001; FDR 4.93e−9 | 0.704; FDR 9.86e−11 | −0.767; FDR 2.34e−5 | −1.122; FDR 0.0090 | −0.538; FDR 0.284 | B |
| Endothelial activation/injury | 0.858; FDR 2.45e−8 | 0.496; FDR 3.18e−10 | −0.159; FDR 0.165 | −0.539; FDR 0.103 | −0.501; FDR 0.284 | B |
| Tissue injury/repair | 0.027; FDR 0.827 | 0.209; FDR 0.040 | −0.022; FDR 0.882 | −0.046; FDR 0.892 | −0.521; FDR 0.284 | C |

## Interpretation

- DSA-positive rejection shows a large IFN/chemokine and endothelial response relative to
  DSA-negative non-rejection.
- DSA-positive/no-rejection versus DSA-negative/no-rejection yielded no FDR-significant module after
  multiplicity correction. This does not establish accommodation, but it rejects a simple assumption
  that DSA positivity alone produces the full bulk rejection program.
- All major immune modules replicate in ABMR versus non-rejection.
- FCGR3A/Fc, cytotoxicity and IFN/chemokine modules are significantly lower in ABMR than TCMR in
  GSE36059. Bulk immune activation therefore cannot be presented as ABMR-specific.
- At week 24, felzartamab is associated with placebo-adjusted reductions in NK identity,
  cytotoxicity and IFN/chemokine modules. The FCGR3A/Fc module is not significant.
- At week 52, none of the nine modules has an FDR-significant treatment difference. This is compatible
  with loss of the short-term immune effect, but it is not proof of recurrence in individual patients.
- Endothelial activation is directionally reduced at weeks 24 and 52 but is not statistically resolved.
  The analysis therefore cannot claim persistent endothelial injury or equivalence between arms.

## What can be presented

The appropriate analysis message is:

> Cross-cohort transcriptomic modules robustly capture rejection-associated immune and endothelial
> activity, but fail ABMR-versus-TCMR specificity. Short-term CD38-directed treatment preferentially
> suppresses NK/cytotoxic/IFN programs. These findings motivate cell-state localization and
> reference-aware filtering before expression-reversal drug prioritization.

## What cannot be presented

- The modules are not validated clinical biomarkers.
- Bulk FCGR3A expression does not establish NK-cell or monocyte causality.
- The analysis does not measure antibody binding to endothelium.
- The randomized cohort is too small for a treatment-response ML model.
- A nonsignificant endothelial treatment contrast does not demonstrate biological persistence.
- No drug candidate, combination, QSAR model, docking result or efficacy claim was generated.

## Next evidence gate

Cell-resolved pseudobulk analysis must determine whether the bulk FCGR3A, cytotoxic and IFN signals
arise from NK cells, monocytes, CD8 T cells or mixtures. Only cell-state programs that remain robust
under stable-allograft, non-rejection injury and native-kidney reference definitions should enter a
drug-reversal query.
