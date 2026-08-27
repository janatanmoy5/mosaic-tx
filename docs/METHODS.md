# Methods

## Analysis populations

GSE50084 biopsy and blood compartments are analyzed separately. Primary discovery uses biopsies.
GSE36059 uses ABMR, TCMR and non-rejecting allografts; mixed and nephrectomy groups are excluded from
primary contrasts but retained for sensitivity analyses. GSE275824 retains every available randomized
biopsy and blocks repeated observations by public patient identifier.

## Statistical rules

The prespecified mechanistic chain is DSA → endothelial engagement → FCGR3A-positive NK/monocyte
activation → cytotoxic and IFNG/chemokine response → microvascular endothelial injury → longitudinal
persistence or response. Operational definitions and stop rules are frozen in `MECHANISM_MAP.md` and
`configs/mechanism_modules.json`.


- Expression: platform-specific gene mapping followed by limma moderated models.
- Multiplicity: Benjamini–Hochberg FDR within each declared contrast family.
- Enrichment: preranked Hallmark GSEA; report NES, FDR and leading edge.
- Longitudinal treatment: `(follow-up − baseline) felzartamab − (follow-up − baseline) placebo`, with
  patient blocking. Published scores use mixed-effects models with patient random intercepts.
- Replication: same-direction effect plus validation-cohort FDR; pathway replication is reported
  separately from gene replication.
- Biological replicate: sample/patient. Cells, spots and ROIs are never substituted for patients.

## Arrow-level tests

- DSA is a clinical/serological exposure from metadata, never an expression-derived label.
- Antibody engagement is not directly observed; it requires DSA context plus endothelial response and is labeled indirect.
- FCGR3A signaling is scored separately in NK and monocyte patient pseudobulk.
- Cytotoxic/IFNG programs must be tested against TCMR and, where annotated, CMV/BK infection.
- Microvascular injury integrates endothelial pseudobulk with Banff MVI rather than expression alone.
- Longitudinal persistence uses treatment-by-time estimates and confidence intervals; nonsignificance alone is insufficient.
- Each arrow is classified as supported, unsupported, or not testable.

## Prohibited analytical shortcuts

No pooling of blood and biopsy, ABMR and TCMR, mixed and pure rejection, stable allograft and
nephrectomy, or cells across patients as independent replicates. No feature selection before data
splitting. No interpretation of GSEA as causal evidence. No clinical prediction claim from public
retrospective expression alone.

## Data processing and executable workflow

The pipeline compiles public GEO matrices and metadata with dataset-specific labeling rules. Platform-specific Bioconductor packages map probes to gene symbols; unmapped probes are removed and multiple probes per gene are averaged. Final gene matrices contain 19,998 genes for GSE50084, 21,358 for GSE36059, and 19,434 for GSE275824.

Cross-sectional differential expression uses limma moderated models. Longitudinal GSE275824 analysis uses patient blocking with duplicateCorrelation and treatment-by-time difference-in-differences. Preranked Hallmark GSEA uses fgsea. Deposited MMDx variables are analyzed by patient-random-intercept mixed models and are not reconstructed or renamed.

Mechanism modules are declared in configs/mechanism_modules.json. Genes are z-standardized within each eligible cohort, and each sample score is the mean of observed module genes subject to minimum coverage. Scores are dimensionless research quantities without clinical cutoffs.

PCA is performed independently by platform on 1,000 most-variable genes. It is descriptive and does not enter hypothesis testing or classification.

The complete workflow runs with:

```bash
bash scripts/run_three_dataset_pipeline.sh
```

The runner compiles data, annotates genes, fits all declared models, generates figures and tables, runs unit and syntax checks, validates output dimensions, and records SHA-256 checksums and a complete execution log.

## Outputs and auditability

Every included public profile and its metadata-derived analysis label are recorded in results/tables/Table_S1_sample_labels.csv. Primary effects and confidence intervals are in Table_2_primary_estimands.csv. Module definitions and observed coverage are in Table_S2_module_definitions_and_coverage.csv. Exact plotting values and PCA coordinates are retained under results/figures/.

## Ethics and interpretation

Only public, de-identified data are used. This secondary analysis involves no participant contact or protected health information. Bulk expression cannot establish antibody binding, cell origin, target causality, drug efficacy, clinical prediction, or treatment guidance.
