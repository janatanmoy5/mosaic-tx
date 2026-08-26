# Scientific protocol

## Objective and falsifiable hypothesis

Develop and externally validate a model that ranks interventions for molecularly defined kidney
allograft injury endotypes. The primary hypothesis is that cell-state-aware, safety-constrained
ranking predicts ex vivo molecular rescue better than bulk signature reversal and target-only
baselines.

The project fails its central claim if it does not improve the preregistered primary metric in a
locked external cohort and a blinded perturbation experiment.

## Intended population and scope

The first study is restricted to adult kidney recipients with biopsy adjudication under the
contemporaneous Banff classification. ABMR, TCMR, mixed rejection, BK virus nephropathy, other
injury, and stable-allograft controls are kept distinct. Pediatric and other-organ generalization
are later studies, not assumed transportable.

## Work packages

### 1. Cohort and reference construction

- Freeze eligibility, index time, treatment exposure windows, and outcome definitions before analysis.
- Build stable, rejection, infection, and nonimmune-injury references matched on organ, time after
  transplant, tissue handling, and platform where possible.
- Track center, ancestry, sex, age, donor type, ischemia, immunosuppression, viral status, DSA,
  histology, kidney function, and pre-biopsy treatment as confounders or audit strata.
- Split by patient and center; all longitudinal samples from one recipient stay in one partition.

### 2. Multimodal state estimation

- Derive donor-aware cell states from single-cell/single-nucleus RNA, spatial assays, and validated
  protein panels; quantify ambient RNA, doublets, batch effects, and cell-composition uncertainty.
- Learn pathway-level states rather than relying on unstable individual-gene lists.
- Anchor biopsy states to less-invasive blood, urine, antibody, and donor-derived cell-free DNA
  measurements without claiming interchangeability.
- Use spatial neighborhoods to weight interactions involving endothelium, Fc-receptor-positive
  innate cells, plasma/B cells, tubules, and stromal injury programs.

### 3. Counterfactual evidence layer

- Curate perturbational evidence with cell context, dose, duration, direction, provenance, and
  direct-versus-inferred status.
- Map targets to cell-state pathways and ligand-receptor circuits; keep causal, associative, and
  predicted edges visibly separate.
- Estimate intervention effects with ensemble models and out-of-distribution checks. Do not treat
  a transcriptomic anticorrelation as proof of therapeutic benefit.

### 4. Safety and feasibility constraints

- Predefine hard exclusions and soft penalties for infection, marrow toxicity, nephrotoxicity,
  hepatic risk, malignancy, QT risk, pregnancy, and interactions with immunosuppression.
- Keep these rules versioned and pharmacist/clinician reviewed. Never infer a dose from this model.
- Report molecular score, evidence coverage, uncertainty, and each penalty separately.

### 5. Validation ladder

1. **Retrospective computational:** center-held-out discrimination, calibration, sensitivity to
   reference definition, negative controls, and comparison with bulk reversal, target proximity,
   and clinical-risk-only baselines.
2. **Orthogonal biological:** concordance across RNA, protein, spatial localization, and independent
   cohorts; leave-one-modality-out ablations.
3. **Blinded experimental:** patient-derived cells, precision-cut tissue, or organoid/co-culture
   systems with prespecified molecular rescue, viability, and functional readouts.
4. **Prospective observational:** lock the pipeline and assess feasibility and temporal prediction
   without changing care.
5. **Interventional:** only after regulatory, ethics, pharmacology, and clinical-trial review.

## Endpoints and statistics

Primary experimental endpoint: difference between ranked and control interventions in a blinded,
prespecified composite rescue score, measured in an external biological system. Secondary endpoints
include cell-state-specific rescue, viability, endothelial injury, complement/Fc programs, and
between-patient heterogeneity.

Use nested, grouped resampling; center-level external validation; bootstrap confidence intervals;
multiple-testing control; calibration curves; and decision-curve analysis only for a separately
defined clinical-use question. Report missingness and exclusions. Sample size and the clinically
meaningful effect must be set with a biostatistician before labels are unlocked.

## Stop/go gates

- **Go to model development:** adequate case/control diversity and complete provenance.
- **Go to wet lab:** candidate stable across reasonable references, models, and perturbation sources.
- **Go to prospective study:** externally replicated effect and acceptable safety/pharmacy review.
- **Stop:** leakage, unresolvable treatment confounding, poor reference transportability, inadequate
  overlap, or failure of blinded experimental superiority.

