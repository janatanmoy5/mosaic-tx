# Evidence-gated single-cell-to-therapeutics extension

## Purpose

This extension asks whether the replicated ABMR immune-effector program can be assigned to specific
human cell states, converted into a patient-level biomarker, and connected to a tractable intervention.
It is not a linear promise that transcriptomics will yield a drug. Every transition has a stop rule.

## Cohort and platform rules

- The new single-cell discovery branch accepts **human Illumina NovaSeq 6000** studies only.
- Platform is verified from repository metadata, not inferred from filenames.
- Stable kidney allograft, non-rejection acute kidney injury, and native healthy kidney are separate
  reference domains. Native healthy tissue is never pooled with stable allograft as one negative class.
- ABMR, mixed ABMR, TCMR, DSA-positive/no rejection, and non-rejection are retained as distinct labels.
- The patient or biopsy—not the cell—is the biological replicate and the unit of data splitting.
- Existing microarray cohorts remain essential external validation and are not excluded by the
  NovaSeq rule, which applies only to the new single-cell branch.

## Stage 1 — Build auditable single-cell objects

The canonical representation is `AnnData` (`.h5ad`) created with Scanpy. A MuData (`.h5mu`) wrapper
is emitted when `muon` is installed so additional modalities can be added later without changing the
RNA object. Required metadata include `dataset_id`, `sample_id`, `patient_id`, `phenotype`, `tissue`,
`species`, `platform`, and `reference_domain`.

Raw counts are preserved in `layers['counts']`. QC, normalization, highly variable genes, PCA,
neighbors, UMAP and Leiden clustering are recorded with versions and parameters. Doublet detection
and ambient-RNA correction are dataset-specific sensitivity analyses, not silent defaults.

Cell labels require concordance of:

1. canonical positive and negative markers;
2. reference mapping or label transfer;
3. cluster-level differential expression;
4. manual review of ambiguous NK/T and monocyte/NK boundaries.

KIR expression alone is not an NK-cell label because KIR-positive CD8 T cells exist. The primary NK
program requires a coherent combination of `NKG7`, `KLRD1`, `NCR1`, `GNLY`, `PRF1`, `FCGR3A` and
absence of a dominant `CD3D`/`CD3E`/`TRAC` program.

## Stage 2 — Cell-state and pathway inference

Differential expression is performed on patient-level pseudobulk counts within prespecified cell
classes. Per-cell tests are descriptive only. Primary cell classes are FCGR3A-positive NK cells,
FCGR3A-positive monocytes, CD8 T cells, glomerular endothelial cells and peritubular endothelial
cells. Models include dataset and relevant clinical covariates when available.

Preranked GSEA uses Hallmark, Reactome and GO Biological Process collections. The prespecified
mechanisms are Fc receptor/ADCC, NK cytotoxicity, IFN-gamma response, endothelial activation,
microvascular injury and tissue-repair programs. Pathway claims require FDR control and directional
replication in an independent cohort.

## Stage 3 — Integrate bulk and single-cell evidence

Single-cell data are used to define cellular origin and coherent modules. Bulk cohorts provide the
patient numbers for replication. Integration uses gene-level effect-size concordance, module scoring,
and pseudobulk-to-bulk transfer; raw single cells and bulk samples are never concatenated into one
expression matrix.

A candidate biomarker must satisfy all of the following:

- directionally concordant in at least two independent patient cohorts;
- associated with ABMR versus non-rejection and evaluated versus TCMR;
- localized to a biologically plausible cell state;
- robust to leave-one-dataset-out analysis;
- not driven only by CMV/BK infection, library depth, platform, or one patient;
- measurable in the intended specimen type.

## Stage 4 — Machine learning, after feature freezing

The prediction unit is the patient/biopsy. Candidate features are frozen before model comparison and
limited to replicated genes, pathway scores and cell-state proportions. Models are intentionally
simple: penalized logistic regression is primary; linear SVM and tree boosting are secondary.

Evaluation uses nested, group-aware cross-validation with all samples from a patient in one fold.
Where several cohorts exist, leave-one-study-out validation is preferred. Feature selection,
normalization and tuning occur inside each training fold. Report AUROC, AUPRC, calibration slope,
Brier score, sensitivity and specificity with bootstrap confidence intervals. A model is rejected if
performance collapses in an untouched study, calibration is poor, or selected features are unstable.
Deep learning is not justified by the available number of independent patients.

## Stage 5 — Drug-reversal analysis

Only replicated, directionally interpretable disease signatures enter perturbational matching.
Up- and down-regulated genes are queried separately against public perturbation resources such as
LINCS/CMap. Ranking must account for cell line, dose, exposure time, mechanism, renal/transplant
relevance, known immunosuppressive toxicity and whether reversal is observed across multiple
contexts. A negative connectivity score is hypothesis generation, not evidence of efficacy.

Candidates advance only if their targets lie in the supported NK–Fc–endothelial network and the
direction of pharmacology is compatible with the observed disease state. Existing approved drugs
and investigational compounds are reported separately.

## Stage 6 — Target and QSAR gate

QSAR is appropriate only when all gates pass:

1. a single molecular target is prioritized by replicated human evidence;
2. activation or inhibition direction is known;
3. ChEMBL or an equivalent source contains a sufficiently large, comparable quantitative activity
   series for that target and assay type;
4. chemical structures and activities can be standardized without excessive censoring;
5. scaffold-aware train/test splitting is possible;
6. applicability domain and prospective value can be assessed.

The QSAR endpoint is one target-specific activity—not ABMR diagnosis and not a mixture of IC50,
EC50 and binding assays. Baselines include Morgan-fingerprint random forest and regularized models.
Evaluation uses scaffold splits, repeated nested validation, Y-randomization, applicability-domain
analysis and an untouched external/test set. Docking or molecular dynamics may explain a QSAR lead
only after redocking/reference-ligand validation; they cannot validate therapeutic efficacy.

## Evidence ladder and stop rules

| Transition | Evidence required | Stop condition |
|---|---|---|
| Cells to signature | donor-level pseudobulk replication | per-cell-only significance |
| Signature to biomarker | independent study and specificity | batch-driven or TCMR-nonspecific |
| Biomarker to model | grouped nested CV and calibration | unstable features or no external performance |
| Signature to drug | robust reversal plus target plausibility | cell-line-only reversal or wrong direction |
| Target to QSAR | adequate standardized activity series | too few compounds or mixed endpoints |

No stage is described as clinical validation. Wet-lab, pharmacology, toxicity and prospective
clinical testing remain necessary after this computational work.
