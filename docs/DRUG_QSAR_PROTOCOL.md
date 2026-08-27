# Drug-reversal and QSAR protocol

## Two different questions

Drug-signature reversal asks whether a perturbation opposes a disease expression program. QSAR asks
whether chemical structure predicts activity against one defined molecular target. They are not one
model and must not share an outcome label.

## Drug-reversal inputs and output

The input signature is created only from genes that replicate directionally in bulk and patient-level
pseudobulk analysis. Separate up/down gene lists, effect sizes, FDR values, tissue, phenotype and cell
state are versioned. LINCS/CMap matches are summarized across cell line, dose and time. The output is
a ranked hypothesis table with compound, perturbation ID, connectivity, supported target, direction,
approval status, transplant-relevant safety concerns and evidence tier.

Compounds are removed when reversal depends on one context, known pharmacology would reinforce the
disease direction, or the molecule is incompatible with immunosuppressed kidney recipients. This is
prioritization for experimental follow-up, not treatment recommendation.

## Target selection

A target advances only when it is expressed in the implicated cell state, directionally associated
with disease in at least two cohorts, altered consistently under treatment or perturbation, and has a
plausible causal position rather than being solely a terminal injury marker. Genetic or orthogonal
protein evidence raises confidence but is not fabricated when absent.

## QSAR data contract

- one target and one mechanism per model;
- human target preferred;
- exact structures, standardized salts/charges/tautomers and duplicate handling recorded;
- comparable assay type and units;
- censoring (`>`, `<`) retained and handled explicitly;
- publication/assay provenance retained;
- no random split as the primary estimate—use Bemis–Murcko scaffold splits;
- duplicates and close analogues cannot cross train/test partitions.

Minimum sample size is determined after chemical-series and endpoint audit. If the activity landscape
cannot support a meaningful scaffold-held-out test, the honest output is “QSAR not feasible.”

## Model evaluation

Start with fingerprint similarity, regularized regression and random forest. Compare more complex
models only if they improve repeated scaffold-held-out performance. Report MAE/RMSE and R² for
continuous activity or AUROC/AUPRC/MCC for a justified activity threshold, plus applicability domain,
Y-randomization and uncertainty. The final test set remains untouched until model selection ends.

Docking is optional structural interpretation. It requires a suitable structure, prepared receptor,
known ligand, redocking benchmark and sensitivity to protonation/site definition. Molecular dynamics
is reserved for a small number of already-supported complexes and does not rescue a weak QSAR model.
